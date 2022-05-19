#include "summary.hpp"
#include "introduce.hpp" //for date parsing functions.
#include "translate.hpp"
#include <boost/date_time/gregorian/gregorian.hpp>

po::variables_map parse_summary_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("summary options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
    ("input-gtf,g", po::value<std::string>()->default_value(""),
     "Input GTF annotation file. Required for --translate / -t")
    ("input-fasta,f", po::value<std::string>()->default_value(""),
     "Input FASTA reference sequence. Required for --translate / -t")
    ("output-directory,d", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("samples,s", po::value<std::string>()->default_value(""),
     "Write a tsv listing all samples in the tree, their parsimony score (terminal branch length), and the ID of their parent node.")
    ("clades,c", po::value<std::string>()->default_value(""),
     "Write a tsv listing all clades and the (inclusive and exclusive of nested clades) count of associated samples.")
    ("sample-clades,C", po::value<std::string>()->default_value(""),
     "Write a tsv of all samples and their associated clade values")
    ("mutations,m", po::value<std::string>()->default_value(""),
     "Write a tsv listing all mutations in the tree and their occurrence count.")
    ("translate,t", po::value<std::string>()->default_value(""),
     "Write a tsv listing the amino acid and nucleotide mutations at each node.")
    ("aberrant,a", po::value<std::string>()->default_value(""),
     "Write a tsv listing duplicate sample identifiers and internal nodes with no mutations and/or branch length 0.")
    ("haplotype,H", po::value<std::string>()->default_value(""),
     "Write a tsv listing haplotypes represented by comma-delimited sets of mutations and their total frequency across the tree.")
    ("calculate-roho,R", po::value<std::string>()->default_value(""),
     "Write a tsv containing the distribution of RoHO values calculated for all homoplasic mutations.")
    ("get-all-basic,A", po::bool_switch(),
     "Use default filenames and save samples, clades, mutations, and aberrant summary tables to the output directory.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("expanded-roho,E", po::bool_switch(),
     "Use to include date and other contextual information in RoHO table output. Significantly slows calculation time.")
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void write_sample_table(MAT::Tree* T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing samples to output %s\n", filename.c_str());
    std::ofstream samplefile;
    samplefile.open(filename);
    //print a quick column header (makes this specific file auspice compatible also!)
    samplefile << "sample\tparsimony\tparent_id\n";

    auto dfs = T->depth_first_expansion();
    for (auto s: dfs) {
        if (s->is_leaf()) {
            //leaves are samples (on the uncondensed tree)
            samplefile << s->identifier << "\t" << s->mutations.size() << "\t" << s->parent->identifier << "\n";
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_clade_table(MAT::Tree* T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing clades to output %s\n", filename.c_str());
    std::ofstream cladefile;
    cladefile.open(filename);
    cladefile << "clade\tinclusive_count\texclusive_count\n";
    //clades will be a map object.
    std::map<std::string, size_t> incl_cladecounts;
    std::map<std::string, size_t> excl_cladecounts;
    for (auto s: T->get_leaves()) {
        bool first_encountered_a1 = true;
        bool first_encountered_a2 = true;
        for (auto a: T->rsearch(s->identifier, false)) {
            std::vector<std::string> canns = a->clade_annotations;
            if (canns.size() >= 1) {
                if (canns[0] != "" && incl_cladecounts.find(canns[0]) == incl_cladecounts.end()) {
                    //these should always be matched in keys.
                    //if a sample is present in one, its present in the other, and vice versa
                    incl_cladecounts[canns[0]] = 0;
                    excl_cladecounts[canns[0]] = 0;
                }
                if (canns[0] != "") {
                    incl_cladecounts[canns[0]]++;
                    if (first_encountered_a1) {
                        excl_cladecounts[canns[0]]++;
                        first_encountered_a1 = false;
                    }
                }
                if (canns.size() >= 2) {
                    if (canns[1] != "" && incl_cladecounts.find(canns[1]) == incl_cladecounts.end()) {
                        incl_cladecounts[canns[1]] = 0;
                        excl_cladecounts[canns[1]] = 0;
                    }
                    if (canns[1] != "") {
                        incl_cladecounts[canns[1]]++;
                        if (first_encountered_a2) {
                            excl_cladecounts[canns[1]]++;
                            first_encountered_a2 = false;
                        }
                    }
                }
            }
        }
    }
    //write the contents of map to the file.
    for (auto const &clade : incl_cladecounts) {
        cladefile << clade.first << "\t" << clade.second << "\t" << excl_cladecounts[clade.first] << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_mutation_table(MAT::Tree* T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing mutations to output %s\n", filename.c_str());
    std::ofstream mutfile;
    mutfile.open(filename);
    mutfile << "ID\toccurrence\n";
    //mutations will be a map object, similar to clades, since we're looking for counts
    std::map<std::string, int> mutcounts;

    auto dfs = T->depth_first_expansion();
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

void write_translate_table(MAT::Tree* T, std::string output_filename, std::string gtf_filename, std::string fasta_filename) {
    fprintf(stderr, "Writing aa and nt mutations per node to output %s\n", output_filename.c_str());
    translate_main(T, output_filename, gtf_filename, fasta_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

std::map<std::set<std::string>,size_t> count_haplotypes(MAT::Tree* T) {
    std::map<std::set<std::string>,size_t> hapmap;
    //naive method. for each sample, rsearch along the history to collect each of its mutations
    //and add those to the set. At the end, add that set to the hapcount keys if its not there, and increment.
    for (auto s: T->get_leaves()) {
        std::set<std::string> mset;
        for (auto a: T->rsearch(s->identifier, true)) {
            for (auto m: a->mutations) {
                mset.insert(m.get_string());
            }
        }
        if (hapmap.find(mset) == hapmap.end()) {
            hapmap[mset] = 0;
        }
        hapmap[mset]++;
    }
    return hapmap;
}

void write_haplotype_table(MAT::Tree* T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing haplotype frequencies to output %s\n", filename.c_str());
    std::ofstream hapfile;
    hapfile.open(filename);
    hapfile << "mutation_set\tsample_count\n";
    auto hapmap = count_haplotypes(T);
    for (auto const &hapc : hapmap) {
        std::ostringstream msetstr;
        for (auto m: hapc.first) {
            msetstr << m << ",";
        }
        std::string final_str = msetstr.str();
        final_str.pop_back();
        if (final_str.c_str()[0] == '\0') {
            final_str = "reference";
        }
        hapfile << final_str.c_str() << "\t" << hapc.second << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_aberrant_table(MAT::Tree* T, std::string filename) {
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
    size_t num_annotations = T->get_num_annotations();
    std::set<std::string> dup_tracker;
    auto dfs = T->depth_first_expansion();
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

void write_sample_clades_table (MAT::Tree* T, std::string sample_clades) {
    timer.Start();
    fprintf(stderr, "Writing clade associations for all samples to output %s\n", sample_clades.c_str());
    std::ofstream scfile;
    scfile.open(sample_clades);
    size_t num_annotations = T->get_num_annotations();
    scfile << "sample";
    for (size_t i = 1;  i <= num_annotations;  i++) {
        scfile << "\tannotation_" << i;
    }
    scfile << "\n";
    for (auto n: T->get_leaves()) {
        std::vector<std::string> annotations_found (num_annotations, "None");
        for (auto a: T->rsearch(n->identifier, false)) {
            std::vector<std::string> canns = a->clade_annotations;
            for (size_t i = 0; i < num_annotations; i++) {
                if (canns[i] != "" && annotations_found[i] == "None") {
                    annotations_found[i] = canns[i];
                }
            }
            bool all_found = true;
            for (size_t i = 0;  i < num_annotations;  i++) {
                if (annotations_found[i] == "None") {
                    all_found = false;
                    break;
                }
            }
            if (all_found) {
                break;
            }
        }
        scfile << n->identifier;
        for (size_t i = 0;  i < num_annotations;  i++) {
            scfile << "\t" << annotations_found[i];
        }
        scfile << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    scfile.close();
}

void write_roho_table(MAT::Tree* T, std::string roho_file, bool get_dates) {
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
    rhfile << "mutation\tparent_node\tchild_count\toccurrence_node\toffspring_with\tmedian_offspring_without\tsingle_roho";
    //the above are the basic columns, the following are additional columns requested to help sort out roho output.
    if (get_dates) {
        rhfile << "\tsister_clade_offspring_counts\tidentical_sample_sibling_count\tearliest_date\tlatest_date\tearliest_identical_sibling\tlatest_identical_sibling\tearliest_clade_sibling_dates\tlatest_clade_sibling_dates\n";
    } else {
        rhfile << "\n";
    }
    for (auto n: T->depth_first_expansion()) {
        //candidate mutations maps each mutation to its child where it occurred
        //child counter maps the number of offspring of a given child
        //together, they store the results we need.
        std::map<std::string,std::string> candidate_mutations;
        std::map<std::string,std::set<std::string>> child_counter;
        //using a separate size tracker to speed things when extra columns are not requested and therefore
        //set membership is irrelevant and wasteful to compute/store
        std::map<std::string,size_t> child_increment;
        std::set<std::string> parent_identical_samples;

        //step one: collect all potential candidate nodes.
        std::vector<std::string> ccheck;
        std::map<std::string,size_t> blentracker;

        for (auto c: n->children) {
            if (!c->is_leaf()) {
                //mutations occurring on any non-leaf child are potentially valid RoHo targets for this node.
                ccheck.push_back(c->identifier);
                blentracker[c->identifier] = c->mutations.size();
                for (auto m: c->mutations) {
                    //assumption each mutation only occurs for one child. This is a good assumption, as if its broken, the tree is malformed.
                    candidate_mutations[m.get_string()] = c->identifier;
                }
            } else if (c->mutations.size() == 0) {
                //leaf children with zero mutations are of interest for date tracking downstream.
                parent_identical_samples.insert(c->identifier);
            }
        }
        if (candidate_mutations.size() == 0) {
            continue;
        }
        //step two: remove candidates who do not 1. have at least two offspring associated with them 2. do not recur secondarily at any point after this parent
        for (auto c: n->children) {
            std::set<std::string> child_samples;
            size_t ccount = 0;
            if (c->is_leaf()) {
                continue;
            }
            for (auto dn: T->depth_first_expansion(c)) {
                if (dn->identifier == c->identifier) {
                    continue;
                }
                if (dn->is_leaf()) {
                    ccount++;
                    if (get_dates) {
                        child_samples.insert(dn->identifier);
                    }
                }
                for (auto m: dn->mutations) {
                    candidate_mutations.erase(m.get_string());
                }
            }
            //don't bother recording values of 0 or 1.
            if (ccount > 1) {
                child_increment[c->identifier] = ccount;
                if (get_dates) {
                    child_counter[c->identifier] = std::move(child_samples);
                }
            }
        }
        //we need at least one valid candidate and at least two non-leaf children to continue.
        if ((candidate_mutations.size() == 0) || (child_increment.size() <= 1)) {
            continue;
        }
        //step 2.5; prerecord date information for each child of node N if requested.
        std::map<std::string, std::pair<boost::gregorian::date,boost::gregorian::date>> datemap;
        std::pair<boost::gregorian::date,boost::gregorian::date> parent_identical_dates;
        if (get_dates) {
            for (auto cc: child_counter) {
                auto cdates = get_nearest_date(T, n, &cc.second);
                datemap[cc.first] = std::move(cdates);
            }
            parent_identical_dates = get_nearest_date(T, n, &parent_identical_samples);
        }

        //step 3: actually record the results.
        for (auto ms: candidate_mutations) {
            //size_t non_c = 0;
            //size_t sum_non = 0;
            std::vector<size_t> all_non;
            size_t sum_wit = 0;
            for (auto cs: child_increment) {
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
            if (all_non.size() %2 == 0) {
                med_non = (all_non[all_non.size()/2-1] + all_non[all_non.size()/2]) / 2;
            } else {
                med_non = all_non[all_non.size()/2];
            }
            std::stringstream nonstrs;
            std::stringstream nonearlydates;
            std::stringstream nonlatedates;
            if (get_dates) {
                for (auto cc: child_counter) {
                    //build all at once to assert same order
                    //auto dates = get_nearest_date(&T, n, &cc.second);
                    if (cc.first != ms.second) {
                        nonstrs << cc.second.size() << ",";
                        auto dates = datemap[cc.first];
                        nonearlydates << dates.first << ",";
                        nonlatedates << dates.second << ",";
                    }
                }
            }

            rhfile << ms.first << "\t" << n->identifier << "\t" << ccheck.size() << "\t" << ms.second << "\t" << sum_wit << "\t" << med_non << "\t" << std::log10(sum_wit/med_non) << "\t";
            if (get_dates) {
                std::string ns = nonstrs.str();
                ns.pop_back();
                rhfile << ns << "\t" << parent_identical_samples.size() << "\t";
                auto descendent_dates = datemap[ms.second];
                rhfile << descendent_dates.first << "\t" << descendent_dates.second << "\t";
                if (parent_identical_samples.size() > 0) {
                    rhfile << parent_identical_dates.first << "\t" << parent_identical_dates.second << "\t";
                } else {
                    rhfile << "None\tNone\t";
                }
                std::string ned = nonearlydates.str();
                std::string nld = nonlatedates.str();
                ned.pop_back();
                nld.pop_back();
                rhfile << ned << "\t" << nld << "\n";
            } else {
                rhfile << "\n";
            }
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    rhfile.close();
}

void summary_main(po::parsed_options parsed) {
    po::variables_map vm = parse_summary_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    bool get_dates = vm["expanded-roho"].as<bool>();

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
    std::string translate = dir_prefix + vm["translate"].as<std::string>();
    std::string gtf = dir_prefix + vm["input-gtf"].as<std::string>();
    std::string fasta = dir_prefix + vm["input-fasta"].as<std::string>();
    std::string aberrant = dir_prefix + vm["aberrant"].as<std::string>();
    std::string roho = dir_prefix + vm["calculate-roho"].as<std::string>();
    std::string hapfile = dir_prefix + vm["haplotype"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    bool get_all = vm["get-all-basic"].as<bool>();
    if (get_all) {
        //if this is set, overwrite with default names for all output
        //and since these names are non-empty, all will be generated
        samples = dir_prefix + "samples.tsv";
        clades = dir_prefix + "clades.tsv";
        mutations = dir_prefix + "mutations.tsv";
        aberrant = dir_prefix + "aberrant.tsv";
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
        write_clade_table(&T, clades);
        no_print = false;
    }
    if (samples != dir_prefix) {
        write_sample_table(&T, samples);
        no_print = false;
    }
    if (mutations != dir_prefix) {
        write_mutation_table(&T, mutations);
        no_print = false;
    }
    if (translate != dir_prefix) {
        bool quit = false;
        if (gtf == dir_prefix) {
            fprintf(stderr, "ERROR: You must specify a GTF file with -g\n");
            quit = true;
        }
        if (fasta == dir_prefix) {
            fprintf(stderr, "ERROR: You must specify a FASTA reference file with -f\n");
            quit = true;
        }
        if (quit) {
            exit(1);
        }
        write_translate_table(&T, translate, gtf, fasta);
        no_print = false;
    }
    if (aberrant != dir_prefix) {
        write_aberrant_table(&T, aberrant);
        no_print = false;
    }
    if (sample_clades != dir_prefix) {
        write_sample_clades_table(&T, sample_clades);
        no_print = false;
    }
    if (roho != dir_prefix) {
        if (get_dates) {
            fprintf(stderr, "RoHO contextual columns requested; collecting...\n");
        }
        write_roho_table(&T, roho, get_dates);
        no_print = false;
    }
    if (hapfile != dir_prefix) {
        write_haplotype_table(&T, hapfile);
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
        size_t slevel = 0;
        size_t mlevel = 0;
        for (auto s: dfs) {
            nodecount++;
            if (s->is_leaf()) {
                slevel += s->level;
                if (s->level > mlevel) {
                    mlevel = s->level;
                }
                samplecount++;
            }
        }
        fprintf(stdout, "Total Nodes in Tree: %d\n", nodecount);
        fprintf(stdout, "Total Samples in Tree: %d\n", samplecount);
        fprintf(stdout, "Total Condensed Nodes in Tree: %ld\n", num_condensed_nodes);
        fprintf(stdout, "Total Samples in Condensed Nodes: %ld\n", num_condensed_leaves);
        fprintf(stdout, "Total Tree Parsimony: %ld\n", T.get_parsimony_score());
        fprintf(stdout, "Number of Clade Annotations: %ld\n", T.get_num_annotations());
        fprintf(stdout, "Max Tree Depth: %ld\n", mlevel);
        fprintf(stdout, "Mean Tree Depth: %f\n", static_cast<float>(slevel) / samplecount);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}
