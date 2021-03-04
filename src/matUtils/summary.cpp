#include "summary.hpp"

po::variables_map parse_summary_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("summary options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("samples,s", po::value<std::string>()->default_value(""),
        "Write a tsv listing all samples in the tree and their parsimony score (terminal branch length).")
        ("clades,c", po::value<std::string>()->default_value(""),
        "Write a tsv listing all clades and the count of associated samples.")
        ("mutations,m", po::value<std::string>()->default_value(""),
        "Write a tsv listing all mutations in the tree and their occurrence count.")
        ("aberrant,a", po::value<std::string>()->default_value(""),
        "Write a tsv listing duplicate samples and internal nodes with no mutations and/or branch length 0.")
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
            samplefile << s->identifier << "\t" << s->branch_length << "\n";
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
    std::set<std::string> dup_tracker;
    auto dfs = T.depth_first_expansion();
    for (auto n: dfs) {
        if (dup_tracker.find(n->identifier) == dup_tracker.end()) {
            dup_tracker.insert(n->identifier);
        } else {
            badfile << n->identifier << "\tduplicate-node-id\n";
        }
        if (n->branch_length == 0 && !n->is_leaf() && !n->is_root()) {
            badfile << n->identifier << "\tinternal-branch-length-0\n";
        }
        if (n->mutations.size() == 0 && !n->is_leaf() && !n->is_root()) {
            badfile << n->identifier << "\tinternal-no-mutations\n";
        }
    }
}

void summary_main(po::parsed_options parsed) {
    po::variables_map vm = parse_summary_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples = vm["samples"].as<std::string>();
    std::string clades = vm["clades"].as<std::string>();
    std::string mutations = vm["mutations"].as<std::string>();
    std::string aberrant = vm["aberrant"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    bool no_print = true;
    if (clades != "" ) {
        write_clade_table(T, clades);
        no_print = false;
    }
    if (samples != "") {
        write_sample_table(T, samples);
        no_print = false;
    }
    if (mutations != "") {
        write_mutation_table(T, mutations);
        no_print = false;
    }  
    if (aberrant != "") {
        write_aberrant_table(T, aberrant);
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
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}