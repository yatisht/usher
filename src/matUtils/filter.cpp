#include "filter.hpp"

po::variables_map parse_filter_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("filter options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->required(),
         "Output mutation-annotated tree file [REQUIRED]")
        ("samples,s", po::value<std::string>()->required(),
        "File containing names of samples (one per line).")
        ("prune,p", po::bool_switch(),
        "Use to remove the indicated samples from the tree. Default is to return a subtree MAT containing only the indicated samples.")
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

void read_sample_names (std::string sample_filename, std::vector<std::string>& sample_names) {
    std::ifstream infile(sample_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", sample_filename.c_str());
        exit(1);
    }    
    std::string line;

    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() != 1) {
            fprintf(stderr, "ERROR: Incorrect format for input file: %s!\n", sample_filename.c_str());
            exit(1);
        }
        else {
            sample_names.push_back(std::move(words[0]));
        }
    }
    infile.close();
}

void filter_main(po::parsed_options parsed) {
    po::variables_map vm = parse_filter_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string sample_filename = vm["samples"].as<std::string>();
    bool prune = vm["prune"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    std::vector<std::string> sample_names;
    timer.Start();
    fprintf(stderr, "Loading sample names from %s.\n", sample_filename.c_str()); 
    read_sample_names(sample_filename, sample_names);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    MAT::Tree subtree;
    if (prune) {
        subtree = prune_leaves(T, sample_names);
    } else if (sample_names.size() < 10000) {
        //for retaining only a subtree, get_subtree is the most effective method
        //for subtress up to about 10 thousand samples in size; after that, pruning
        //all other nodes becomes more efficient, because get_subtree scales with 
        //the size of the input sample set, while prune takes a similar time for 
        //any sample set size, while scaling on total tree size
        subtree = get_sample_subtree(T, sample_names);
    } else {
        subtree = get_sample_prune(T, sample_names);
    }
    
    timer.Start();
    fprintf(stderr, "Saving output MAT file %s.\n", output_mat_filename.c_str()); 
    // Save output MAT
    subtree.condense_leaves();
    MAT::save_mutation_annotated_tree(subtree, output_mat_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

MAT::Tree prune_leaves (const MAT::Tree& T, std::vector<std::string> sample_names) {

    timer.Start();
    fprintf(stderr, "Pruning %zu samples.\n", sample_names.size());
    auto subtree = MAT::get_tree_copy(T);
    for (auto s: sample_names) {
        if (subtree.get_node(s) == NULL) {
            fprintf(stderr, "ERROR: Sample %s not found in the tree!\n", s.c_str()); 
        }
        else {
            subtree.remove_node(s, true);
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_subtree (const MAT::Tree& T, std::vector<std::string> sample_names) {

    timer.Start();
    fprintf(stderr, "Extracting subtree of %zu samples.\n", sample_names.size());
    auto subtree = MAT::get_subtree(T, sample_names);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_prune (const MAT::Tree& T, std::vector<std::string> sample_names) {
    
    timer.Start();
    fprintf(stderr, "Large sample input; building subtree by pruning all but %zu samples.\n", sample_names.size());
    //its going to be much quicker to check membership in a set then a vector O(logN)
    //which is important when we're looking at large sample inputs
    //so convert the sample names vector to an unordered set
    std::unordered_set<std::string> setnames(sample_names.begin(),sample_names.end());

    auto subtree = MAT::get_tree_copy(T);
    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        //only call the remover on leaf nodes (can't be deleting the root...)
        if (s->is_leaf()) {
            //if the node is NOT in the set, remove it
            if (setnames.find(s->identifier) == setnames.end()) {
                subtree.remove_node(s->identifier, true);
            }
        }
    }    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    return subtree;
}