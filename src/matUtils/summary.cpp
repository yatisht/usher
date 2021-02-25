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

void write_sample_table(const& MAT::Tree T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing samples to output %s", filename.c_str());
    std::ofstream samplefile;
    samplefile.open(filename);
    //print a quick column header (makes this specific file auspice compatible also!)
    samplefile << "sample\tparsimony\n";

    dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        if (s->is_leaf()) {
            //leaves are samples (on the uncondensed tree)
            samplefile << s->identifier << "\t" << s->branch_length << "\n";
        }        
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void write_clade_table(const& MAT::Tree T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing clades to output %s", filename.c_str());
    std::ofstream cladefile;
    cladefile.open(filename);
    cladefile << "clade\tcount\n";
    //clades will be a map object.
    std::map<std::string, int> cladecounts;
    
    dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        if (s->is_leaf()) {
            for (auto c: s->clade_annotations) {
                if (cladecounts.find(c) != cladecounts.end()) {
                    cladecounts[c]++; //increment should work?
                } else {
                    cladecounts[c] = 1;
                }
            }
        }
    }
    //write the contents of map to the file.
    for (auto const &clade : cladecounts) {
        cladefile << clade.first << "\t" << clade.second << "\n";s
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
void write_mutation_table(const& MAT::Tree T, std::string filename) {
    timer.Start();
    fprintf(stderr, "Writing mutations to output %s", filename.cstr());
    std::ofstream mutfile;
    mutfile.open(filename);
    mutfile << "ID\toccurrence\n";
    //mutations will be a map object, similar to clades, since we're looking for counts
    std::map<std::string, int> mutcounts;

    dfs = T.depth_first_expansion();
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
            if (mutcounts.find(c) != mutcounts.end()) {
                mutcounts[c]++ ;
            } else {
                mutcounts[c] = 1;
            }
        }
    }
    //write the contents of map to the file
    for (auto const &mut : mutcounts) {
        mutfile << mut.first << "\t" << mut.second << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void summary_main(po::parsed_options parsed) {
    po::variables_map vm = parse_summary_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples = vm["samples"].as<std::string>();
    std::string clades = vm["clades"].as<std::string>();
    std::string mutations = vm["mutations"].as<std::string>();

    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    if (clades != "" ) {
        write_clade_table(T, clades);
    }
    if (samples != "") {
        write_sample_table(T, samples);
    }
    if (mutations != "") {
        write_mutation_table(T, mutations);
    }
}