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
    std::vector<std::string> samples = read_sample_names(mutation_paths_filename);
    std::vector<std::string> pathstring = mutation_paths(T, samples);
    for (auto s: pathstring) {
        std::cout << s << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

std::vector<std::string> mutation_paths(const MAT::Tree& T, std::vector<std::string> samples) {
    std::vector<std::string> mpaths;
    for (auto sample: samples) {
        std::string cpath = sample + ": ";
        auto root_to_sample = T.rsearch(sample, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        //fprintf(stdout, "%s: ", sample.c_str());
        for (auto n: root_to_sample) {
            for (size_t i=0; i<n->mutations.size(); i++) {
                cpath += n->mutations[i].get_string();
                //fprintf(stdout, "%s", n->mutations[i].get_string().c_str());
                if (i+1 < n->mutations.size()) {
                    cpath += ",";
                    //fprintf(stdout, ",");
                }
            }
            if (n != root_to_sample.back()) {
                cpath += " > ";
                //fprintf(stdout, " > ");
            }
        }
        mpaths.push_back(cpath);
        //fprintf(stdout, "\n");
    }
    return mpaths;
}

std::vector<std::string> clade_paths(MAT::Tree T, std::vector<std::string> clades) {
    //get the set of clade path strings for printing
    //similar to the above, construct a series of strings to be printed or redirected later on
    std::vector<std::string> clpaths;
    clpaths.push_back("clade\tspecific\tfrom_root");
    //do a breadth-first search
    //the first time a clade that is in clades is encountered, that's the root;
    //get the path of mutations to that root (rsearch), save the unique mutations + that path
    //unique mutations being ones that occurred in the clade root, and the path being all mutations from that root back to the tree root
    //then continue. if a clade has already been encountered in the breadth first, its
    //not clade root, and should be skipped.
    std::unordered_set<std::string> clades_seen;

    auto dfs = T.breadth_first_expansion();
    for (auto n: dfs) {
        std::string curpath = n->identifier + "\t";
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                //if its one of our target clades and it hasn't been seen...
                if (std::find(clades.begin(), clades.end(), ann) != clades.end() && clades_seen.find(ann) != clades_seen.end()){
                    //get its own mutations, if there are any
                    std::string unique = "";
                    for (size_t i=0; i<n->mutations.size(); i++) {
                        unique += n->mutations[i].get_string();
                        if (i+1 < n->mutations.size()) {
                            unique += ",";
                        }
                    }
                    curpath += unique + "\t";
                    //get the ancestral mutations back to the root
                    std::string root = "";
                    auto root_to_sample = T.rsearch(n->identifier, true);
                    for (auto n: root_to_sample) {
                        for (size_t i=0; i<n->mutations.size(); i++) {
                            root += n->mutations[i].get_string();
                            if (i+1 < n->mutations.size()) {
                                root += ",";
                            }
                        }
                        if (n != root_to_sample.back()) {
                            root += " > ";
                        }
                    }
                    //save values to the string, mark the clade as seen, and save the string
                    curpath += unique;
                    curpath += root;
                    clades_seen.insert(ann);
                    clpaths.push_back(curpath);
                }
            }
        }
    }
    return clpaths;
}