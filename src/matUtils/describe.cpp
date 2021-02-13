#include "describe.hpp"

po::variables_map parse_describe_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("describe options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("mutation-paths,m", po::value<std::string>()->required(),
         "File containing sample names for which mutation paths should be displayed [REQUIRED]")
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

void describe_main(po::parsed_options parsed) {
    po::variables_map vm = parse_describe_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string mutation_paths_filename = vm["mutation-paths"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Displaying mutation paths for the samples in the input file %s.\n", mutation_paths_filename.c_str()); 
    mutation_paths(T, mutation_paths_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void mutation_paths(const MAT::Tree& T, std::string mutation_paths_filename) {
    std::ifstream infile(mutation_paths_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", mutation_paths_filename.c_str());
        exit(1);
    }    
    std::string line;

    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() != 1) {
            fprintf(stderr, "ERROR: Incorrect format for input file: %s!\n", mutation_paths_filename.c_str());
            exit(1);
        }
        else {
            auto sample = words[0];
            if (T.get_node(sample) == NULL) {
                fprintf(stderr, "ERROR: Input sample %s not found in the tree!\n", sample.c_str());
                exit(1);
            }
            else {
                auto root_to_sample = T.rsearch(sample, true);
                std::reverse(root_to_sample.begin(), root_to_sample.end());
                fprintf(stdout, "%s: ", sample.c_str());
                for (auto n: root_to_sample) {
                    for (size_t i=0; i<n->mutations.size(); i++) {
                        fprintf(stdout, "%s", n->mutations[i].get_string().c_str());
                        if (i+1 < n->mutations.size()) {
                            fprintf(stdout, ",");
                        }
                    }
                    if (n != root_to_sample.back()) {
                        fprintf(stdout, " > ");
                    }
                }
                fprintf(stdout, "\n");
            }
        }
    }
    infile.close();
}

