#include "summary.hpp"

po::variables_map parse_summary_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("summary options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("output-directory,d", po::value<std::string>()->default_value("./"),
        "Write output files to the target directory. Default is current directory.")
        ("samples,s", po::value<std::string>()->default_value(""),
        "Write a tsv listing all samples in the tree and their parsimony score (terminal branch length).")
        ("clades,c", po::value<std::string>()->default_value(""),
        "Write a tsv listing all clades and the count of associated samples.")
        ("sample-clades,C", po::value<std::string>()->default_value(""),
        "Write a tsv of all samples and their associated clade values")
        ("mutations,m", po::value<std::string>()->default_value(""),
        "Write a tsv listing all mutations in the tree and their occurrence count.")
        ("aberrant,a", po::value<std::string>()->default_value(""),
        "Write a tsv listing duplicate samples and internal nodes with no mutations and/or branch length 0.")
        ("calculate-roho,R", po::value<std::string>()->default_value(""),
        "Write a tsv containing the distribution of ROHO values calculated for all homoplasic mutations.")
        ("get-all,A", po::bool_switch(),
        "Use default filenames (samples.txt, clades.txt, etc) and save all summary tables to the output directory.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void write_sample_table(MAT::Tree& T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing samples to output %s\n", filename.c_str());
    std::ofstream samplefile;
    samplefile.open(filename);
    //print a quick column header (makes this specific file auspice compatible also!)
    samplefile << "sample\tparsimony\n";

    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        if (s->is_leaf()) {
            //leaves are samples (on the uncondensed tree)
            samplefile << s->identifier << "\t" << s->mutations.size() << "\n";
        }        
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_clade_table(MAT::Tree& T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing clades to output %s\n", filename.c_str());
    std::ofstream cladefile;
    cladefile.open(filename);
    cladefile << "clade\tcount\n";
    //clades will be a map object.
    std::map<std::string, size_t> cladecounts;
    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        std::vector<std::string> canns = s->clade_annotations;
        if (canns.size() > 0) {
            if (canns.size() > 1 || canns[0] != "") {
                //the empty string is the default clade identifier attribute
                //skip entries which are annotated with 1 clade but that clade is empty
                //but don't skip entries which are annotated with 1 clade and its not empty
                //get the set of samples descended from this clade root
                for (auto c: canns) {
                    //the empty string is a default clade identifier
                    //make sure not to include it.
                    //the first time you encounter it should be the root of the clade
                    //then if you see if after that you can ignore it. probably?
                    if (cladecounts.find(c) == cladecounts.end() && c != "") {
                        std::vector<std::string> sids = T.get_leaves_ids(s->identifier);
                        cladecounts[c] = sids.size();
                    }
                }        
            }
        }
    }
    //write the contents of map to the file.
    for (auto const &clade : cladecounts) {
        cladefile << clade.first << "\t" << clade.second << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_mutation_table(MAT::Tree& T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing mutations to output %s\n", filename.c_str());
    std::ofstream mutfile;
    mutfile.open(filename);
    mutfile << "ID\toccurrence\n";
    //mutations will be a map object, similar to clades, since we're looking for counts
    std::map<std::string, int> mutcounts;

    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        //occurrence here is the number of times a mutation occurred 
        //during the history of the pandemic
        //while the number of samples involved is interesting, sample
        //bias towards specific geographic regions (cough the UK cough)
        //means its not very useful without normalization; occurrence
        //across the tree will give a better idea about mutations which
        //happen again and again that may be positively selected.
        //want to include internal nodes this time.
        for (auto m: s->mutations) {
            std::string mname = m.get_string();
            if (mname != "MASKED") {
                if (mutcounts.find(mname) != mutcounts.end()) {
                    mutcounts[mname]++ ;
                } else {
                    mutcounts[mname] = 1;
                }
            }
        }
    }
    //write the contents of map to the file
    for (auto const &mut : mutcounts) {
        mutfile << mut.first << "\t" << mut.second << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_aberrant_table(MAT::Tree& T, std::string filename) {
    /*
    This function identifies nodes which have missing information or are otherwise potentially problematic.
    These nodes do not necessarily break matUtils or Usher in general, but they may need to be accounted for
    in downstream analysis. Nodes like these can be created by the masking out of specific samples or by not 
    collapsing trees constructed from bifurcated/fully resolved or other types of tree.
    */
    timer.Start();
    fprintf(stderr, "Writing bad nodes to output %s\n", filename.c_str());
    std::ofstream badfile;
    badfile.open(filename);
    badfile << "NodeID\tIssue\n";
    size_t num_annotations = T.get_num_annotations();
    std::set<std::string> dup_tracker;
    auto dfs = T.depth_first_expansion();
    for (auto n: dfs) {
        if (dup_tracker.find(n->identifier) == dup_tracker.end()) {
            dup_tracker.insert(n->identifier);
        } else {
            badfile << n->identifier << "\tduplicate-node-id\n";
        }
        if (n->mutations.size() == 0 && !n->is_leaf() && !n->is_root()) {
            badfile << n->identifier << "\tinternal-no-mutations\n";
        }
        if (num_annotations != n->clade_annotations.size()) {
            badfile << n->identifier << "\tclade-annotations (" << n->clade_annotations.size() <<
              " not " << num_annotations << ")\n";
        }
    }
}

void write_sample_clades_table (MAT::Tree& T, std::string sample_clades) {
    timer.Start();
    fprintf(stderr, "Writing clade associations for all samples to output %s\n", sample_clades.c_str());
    std::ofstream scfile;
    scfile.open(sample_clades);
    scfile << "sample\tannotation_1\tannotation_2\n";
    for (auto n: T.get_leaves()) {
        std::string ann_1 = "None";
        std::string ann_2 = "None";
        for (auto a: T.rsearch(n->identifier, false)) {
            std::vector<std::string> canns = a->clade_annotations;
            for (size_t i = 0; i < canns.size(); i++) {
                if ((i == 0) && (canns[i] != "") && (ann_1 == "None")) {
                    ann_1 = canns[i];
                }
                if ((i == 1) && (canns[i] != "") && (ann_2 == "None")) {
                    ann_2 = canns[i];
                }
                if ((ann_1 != "None") && (ann_2 != "None")) {
                    break;
                }
            }
        }
        scfile << n->identifier << "\t" << ann_1 << "\t" << ann_2 << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    scfile.close();
}

void write_roho_table(MAT::Tree& T, std::string roho_file) {
    /*
    RoHo is a metric used for SARS-CoV-2 genetic analysis described in van Dorp et al 2021 (doi: 10.1038/s41467-020-19818-2)
    Essentially it is the proportion of offspring of a parent node which have a mutation vs do not have a mutation of interest
    It is 1 when the mutation is irrelevant to growth, with a significant amount of stochasticity to it
    meaning homoplasy is required to replicate and give statistical significance.

    This function calculates these ROHO values for every mutation which occurs across the entire tree, including single occurrences, where 
    possible, to allow for downstream analysis.
    */
    timer.Start();
    fprintf(stderr, "Calculating and writing RoHo values to output %s\n", roho_file.c_str());
    std::ofstream rhfile;
    rhfile.open(roho_file);
    // rhfile << "mutation\tparent_node\tchild_count\toccurrence_node\toffspring_with\tmedian_offspring_without\tmean_offspring_without\tbranch_length_with\tbranch_length_without\n";
    rhfile << "mutation\tparent_node\tchild_count\toccurrence_node\toffspring_with\tmedian_offspring_without\tsingle_roho\n";
    for (auto n: T.depth_first_expansion()) {
        //candidate mutations maps each mutation to its child where it occurred
        //child counter maps the number of offspring of a given child
        //together, they store the results we need.
        std::map<std::string,std::string> candidate_mutations;
        std::map<std::string,size_t> child_counter;
        // //temporary check, let's try skipping nodes with >2 children and see what that does.
        // //there are four alternatives 1. drop these occurrences (lose a lot of tree) 2. pick one at random to compare 3. pick biggest to compare 4. compare to the mean of the others
        // if (n->children.size() > 2) {
        //     continue;
        // }

        //step one: collect all potential candidate nodes.
        size_t ccheck = 0;
        std::map<std::string,size_t> blentracker;

        for (auto c: n->children) {
            if (!c->is_leaf()) {
                //mutations occurring on any non-leaf child are potentially valid RoHo targets for this node.
                ccheck++;
                blentracker[c->identifier] = c->mutations.size();
                for (auto m: c->mutations) {
                    //candidate_mutations.insert(m.get_string());
                    //assumption each mutation only occurs for one child. This is a good assumption, as if its broken, the tree is malformed.
                    candidate_mutations[m.get_string()] = c->identifier;
                }
            }
        }
        if (candidate_mutations.size() == 0) {
            continue;
        }
        //step two: remove candidates who do not 1. have at least two offspring associated with them 2. do not recur secondarily at any point after this parent
        for (auto c: n->children) {
            if (c->is_leaf()) {
                continue;
            }
            size_t cc = 0;
            for (auto dn: T.depth_first_expansion(c)) {
                if (dn->identifier == c->identifier) {
                    continue;
                }
                if (dn->is_leaf()) {
                    cc++;
                }
                for (auto m: dn->mutations) {
                    candidate_mutations.erase(m.get_string());
                    // auto find = candidate_mutations.get(m.get_string());
                    // if (find != candidate_mutations.end()) {
                    //     //remove mutations which recur secondarily at any point in the descendents from consideration.
                    //     candidate_mutations.erase(find);
                    // }
                }
            }
            //don't bother recording values of 0 or 1.
            if (cc > 1) {
                child_counter[c->identifier] = cc;
            }
        }
        //we need at least one valid candidate and at least two non-leaf children to continue.
        if ((candidate_mutations.size() == 0) || (child_counter.size() <= 1)) {
            continue;
        }
        //step 3: actually record the results.
        for (auto ms: candidate_mutations) {
            //size_t non_c = 0;
            //size_t sum_non = 0;
            std::vector<size_t> all_non;
            size_t sum_wit = 0;
            for (auto cs: child_counter) {
                if (cs.first != ms.second) {
                    all_non.push_back(cs.second);
                    //sum_non += cs.second;
                    //non_c++;
                } else {
                    sum_wit += cs.second;
                }
            }
            float med_non;
            std::sort(all_non.begin(), all_non.end());
            if (all_non.size() %2 == 0){
                med_non = (all_non[all_non.size()/2-1] + all_non[all_non.size()/2]) / 2;
            } else {
                med_non = all_non[all_non.size()/2];
            }
            // float mean_non = 0;
            // for (auto nv: all_non) {
                // mean_non += nv;
            // }
            // mean_non = mean_non / all_non.size();
            // size_t non_blc = 0;
            // size_t sum_nonblc = 0;
            // size_t sum_witblc = 0;
            // for (auto bs: blentracker) {
            //     if (bs.first != ms.second) {
            //         sum_nonblc += bs.second;
            //         non_blc++;
            //     } else {
            //         sum_witblc += bs.second;
            //     }
            // }
            //rhfile << ms.first << "\t" << n->identifier << "\t" << ccheck << "\t" << ms.second << "\t" << sum_wit << "\t" << med_non << "\t" << mean_non << "\t" << sum_witblc << "\t" << sum_nonblc/non_blc << "\n";
            rhfile << ms.first << "\t" << n->identifier << "\t" << ccheck << "\t" << ms.second << "\t" << sum_wit << "\t" << med_non << "\t" << std::log10(sum_wit/med_non) << "\n";
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    rhfile.close();
}

void summary_main(po::parsed_options parsed) {
    po::variables_map vm = parse_summary_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";

    std::string samples = dir_prefix + vm["samples"].as<std::string>();
    std::string clades = dir_prefix + vm["clades"].as<std::string>();
    std::string sample_clades = dir_prefix + vm["sample-clades"].as<std::string>();
    std::string mutations = dir_prefix + vm["mutations"].as<std::string>();
    std::string aberrant = dir_prefix + vm["aberrant"].as<std::string>();
    std::string roho = dir_prefix + vm["calculate-roho"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    bool get_all = vm["get-all"].as<bool>();
    if (get_all) {
        //if this is set, overwrite with default names for all output
        //and since these names are non-empty, all will be generated
        samples = dir_prefix + "samples.tsv";
        clades = dir_prefix + "clades.tsv";
        mutations = dir_prefix + "mutations.tsv";
        aberrant = dir_prefix + "aberrant.tsv";
        sample_clades = dir_prefix + "sample-clades.tsv";
        roho = dir_prefix + "roho-values.tsv";
    }
    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    // record the number of condensed leaves in case its needed for printing later.
    size_t num_condensed_leaves = T.condensed_leaves.size();
    size_t num_condensed_nodes = T.condensed_nodes.size();
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    bool no_print = true;
    if (clades != dir_prefix ) {
        write_clade_table(T, clades);
        no_print = false;
    }
    if (samples != dir_prefix) {
        write_sample_table(T, samples);
        no_print = false;
    }
    if (mutations != dir_prefix) {
        write_mutation_table(T, mutations);
        no_print = false;
    }  
    if (aberrant != dir_prefix) {
        write_aberrant_table(T, aberrant);
        no_print = false;
    }
    if (sample_clades != dir_prefix) {
        write_sample_clades_table(T, sample_clades);
        no_print = false;
    }
    if (roho != dir_prefix) {
        write_roho_table(T, roho);
        no_print = false;
    }
        
    if (no_print) {
        //just count the number of nodes in the tree and the number of leaves (samples)
        //additional basic statistics can also go in this code block (total tree depth?)
        timer.Start();
        fprintf(stderr, "No arguments set; getting basic statistics...\n");
        int nodecount = 0;
        int samplecount = 0;
        auto dfs = T.depth_first_expansion();
        for (auto s: dfs) {
            nodecount++;
            if (s->is_leaf()){
                samplecount++;
            }
        }
        fprintf(stdout, "Total Nodes in Tree: %d\n", nodecount);
        fprintf(stdout, "Total Samples in Tree: %d\n", samplecount);
        fprintf(stdout, "Total Condensed Nodes in Tree: %ld\n", num_condensed_nodes);
        fprintf(stdout, "Total Samples in Condensed Nodes: %ld\n", num_condensed_leaves);
        fprintf(stdout, "Total Tree Parsimony: %ld\n", T.get_parsimony_score());
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}
