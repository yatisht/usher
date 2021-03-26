#include "introduce.hpp"

po::variables_map parse_introduce_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description filt_desc("mask options");
    filt_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("population-samples,s", po::value<std::string>()->required(), 
         "Names of samples from the population of interest [REQUIRED].") 
        ("max-mutations,m", po::value<int>()->default_value(3),
        "Maximum number of mutations allowed per branch on path to putative introduction")
        ("output,o", po::value<std::string>()->required(),
        "Name of the file to save the introduction information to.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << filt_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

std::vector<std::string> find_introductions(MAT::Tree* T, std::vector<std::string> samples, int k) {
    /*
    This function applies a five-step heuristic process to identify putative introductions to a population defined by the samples vector
    1. The root is initialized as OUT.
    2. Traversal of nodes in DFS order. Labels of leaf nodes don't change. For internal nodes, if there is a direct leaf descendant with 0 mutations to the internal node labelled IN, the internal node is labelled IN.
    3. If there is no direct leaf descendant labelled as IN, the internal node is labelled IN only if its parent is labelled IN and the closest IN leaf node is at most k mutations from this internal node.
    4. If 2 and 3 are not satisfied, the node is labelled as OUT.
    5. Now all nodes are labelled IN or OUT. The number of introductions is simply the number of OUT->IN transitions in parent-child pairs which can be computed by traversing the tree again in a DFS order.
    */
    std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
    std::map<std::string, int> assignments;
    //the first real step is traversing nodes in DFS order.
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        if (n->is_root()) {
            //rule 1
            assignments[n->identifier] = 0;
        } else if (n->is_leaf()) {
            //rule 2. Leaves are in or out based on membership in samples only.
            if (sample_set.find(n->identifier) != sample_set.end()) {
                assignments[n->identifier] = 1;
            } else {
                assignments[n->identifier] = 0;
            }
        } else {
            //check rule 2
            bool r2 = false;
            for (auto child: n->children) {
                if (child->is_leaf()) {
                    if (sample_set.find(child->identifier) != sample_set.end() & (child->mutations.size() == 0)) {
                        assignments[n->identifier] = 1;
                        r2 = true;
                    }
                }
            }
            //if not assigned via rule 2, check rule 3
            bool r3 = false;
            if (!r2) {
                auto sv = assignments.find(n->parent->identifier);
                if (sv != assignments.end()) {
                    if (sv->second == 1) {
                        r3 = true;
                        assignments[n->identifier] = 1;
                    } else {
                        //now we get the most complex bit of logic
                        //get terminals descended from this node, and for each one, if theyre members of population
                        //then check to see if there are less than k total mutations on the rsearch
                        auto leaves = T->get_leaves(n->identifier);
                        for (auto l: leaves) {
                            if (sample_set.find(l->identifier) != sample_set.end() & !r3) {
                                auto path = T->rsearch(l->identifier, true);
                                size_t mcount = 0;
                                for (auto p: path) {
                                    if (p->identifier == n->identifier & (mcount <= k)) {
                                        assignments[n->identifier] = 1;
                                        r3 = true;
                                        break;
                                    } else if (mcount > k) {
                                        break;
                                    }
                                    mcount += p->mutations.size();
                                }
                            }
                        }
                    }
                }
            }
            if (!r3) {
                assignments[n->identifier] = 0;
            }
        }
    }
    //now that we have assignments in a map objects of node identifiers to states
    //we can traverse again in dfs order and check where it changes from 0 to 1
    //rule 4
    std::vector<std::string> outstrs;
    int last = 0;
    for (auto n: dfs) {
        auto search = assignments.find(n->identifier);
        //every single node should be in assignments.
        assert (search != assignments.end());
        if (last == 0 & search->second == 1) {
            //this is an introduction!
            //collect some info.
            //column 1 is the internal node identifier, column 2 is a comma-delineated list of terminals from
            //this internal node which are IN, and column 3 is a comma-delineated list of terminals which are OUT
            std::string ostr = n->identifier + "\t";
            std::string interms = "";
            std::string oterms = "";
            auto terminals = T->get_leaves_ids(n->identifier);
            for (auto t: terminals) {
                if (assignments.find(t)->second == 1) {
                    interms += "," + t;
                } else {
                    oterms += "," + t;
                }
            }
            ostr += interms + "\t";
            ostr += oterms + "\n";
            outstrs.push_back(ostr);
        }
        last = search->second;
    }
    return outstrs;
}

void introduce_main(po::parsed_options parsed) {
    po::variables_map vm = parse_introduce_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples_filename = vm["population-samples"].as<std::string>();
    int max_mutations = vm["max-mutations"].as<int>();
    std::string output_file = vm["output"].as<std::string>();
    //uint32_t num_threads = vm["threads"].as<uint32_t>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    auto popsamples = read_sample_names(samples_filename);
    auto outstrings = find_introductions(&T, popsamples, max_mutations);

    std::ofstream of;
    of.open(output_file);
    for (std::string o: outstrings) {
        of << o;
    }
    of.close();
}
