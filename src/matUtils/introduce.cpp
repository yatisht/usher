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
std::vector<std::string> find_introductions(MAT::Tree* T, std::vector<std::string> samples) {
    //k is going to be unused for the moment while I test out this new algo
    //if I like it, I'll cut the k argument
    /*
    This function applies a heuristic series of steps to label internal nodes as in or out of a geographic area
    based on their relationship to the samples in the input list. The rules are:
    1. If the node is a leaf, it is IN if it is in samples, otherwise OUT
    2. If all descendents of the node are IN, its IN
    3. If all descendents of the node are OUT, its OUT
    4. Nodes are assigned the state which yields the minimum distance to the nearest leaf of that type divided by the total number of leaves of that type descended from it
    This last step is the most complex one, and is a weighted minimum which takes into account both the number of descendents and the total distance to those descendents
    It's worth noting that if a sample is an identical child to the internal node, that state will always win
    5. On a tie, the node is assigned to the state of its parent

    Introductions are identified as locations where assignments shift from OUT to IN on a depth first search.
    */
    std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
    std::map<std::string, int> assignments;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        if (n->is_leaf()) {
            //rule 1
            if (sample_set.find(n->identifier) != sample_set.end()) {
                assignments[n->identifier] = 1;
            } else {
                assignments[n->identifier] = 0;
            }
        } else {
            auto leaves = T->get_leaves_ids(n->identifier);
            //to apply rules 2-3, we need to check the state of each leaf
            std::unordered_set<std::string> in_leaves;
            std::unordered_set<std::string> out_leaves;
            for (auto l: leaves) {
                if (sample_set.find(l) != sample_set.end()) {
                    in_leaves.insert(l);
                } else {
                    out_leaves.insert(l);
                }
            }
            if (out_leaves.size() == 0) {
                //rule 2
                assignments[n->identifier] = 1;
            } else if (in_leaves.size() == 0) {
                //rule 3
                assignments[n->identifier] = 0;
            } else {
                //rule 4...
                //the best way to do this is to keep in mind that the nearest descendent leaf 
                //is going to be the next leaf encountered in DFS order
                //so we're going to iterate over DFS until we encounter a leaf of each type and record those distances
                size_t min_to_in = 0;
                size_t min_to_out = 0;
                for (auto d: T->depth_first_expansion(n)) {
                    if ((min_to_in > 0) & (min_to_out > 0)){
                        //if both of these are assigned, we're done. break out of this loop and move on
                        break;
                    }
                    if (d->is_leaf()){
                        if ((sample_set.find(d->identifier) != sample_set.end()) & (min_to_in == 0)) {
                            //found the nearest IN.
                            //rsearch back from it so that we're getting the direct path to the original node
                            //maybe there's a better way to do this, but I want to make sure I'm not including side branches
                            //on my iteration
                            size_t total_traveled = 0;
                            for (auto a: T->rsearch(d->identifier, true)) {
                                total_traveled += a->mutations.size();
                                if (a->identifier == n->identifier) {
                                    //back to the original query, break here
                                    min_to_in = total_traveled;
                                    break;
                                }
                            }
                        } else if ((sample_set.find(d->identifier) == sample_set.end()) & (min_to_out == 0)) {
                            //found the nearest out
                            //same routine as above
                            size_t total_traveled = 0;
                            for (auto a: T->rsearch(d->identifier, true)) {
                                total_traveled += a->mutations.size();
                                if (a->identifier == n->identifier) {
                                    //back to the original query, break here
                                    min_to_out = total_traveled;
                                    break;
                                }
                            }
                        } 
                    }
                }
                //now we complete rule 4 by checking the balance.
                if ((min_to_in / in_leaves.size()) < (min_to_out / out_leaves.size())) {
                    //IN wins
                    assignments[n->identifier] = 1;
                } else if ((min_to_in / in_leaves.size()) > (min_to_out / out_leaves.size())) {
                    //OUT wins
                    assignments[n->identifier] = 0;
                } else {
                    //if we end up with a tie, the parent state is tiebreaker.
                    //check if this is the root (this should rarely be the case that its perfectly balanced across the tree)
                    //(but if it is, we're going to default to OUT)
                    if (n->is_root()) {
                        assignments[n->identifier] = 0;
                    } else {
                        auto search = assignments.find(n->parent->identifier);
                        assignments[n->identifier] = search->second;
                    }
                }
            }
        }
    }
    //after the wild ride that was that set of heuristics
    //its time to identify introductions. we're going to do this by iterating
    //over all of the leaves. For each sample, rsearch back until it hits a 0 assignment
    //then record the last encountered 1 assignment as the point of introduction
    std::vector<std::string> outstrs;
    for (auto s: samples) {
        //everything in this vector is going to be 1 (IN)
        std::string last_encountered = s;
        size_t traversed = 0;
        for (auto a: T->rsearch(s,true)) {
            if (a->is_root()) {
                //if we get back to the root, the root is necessarily the point of introduction for this sample
                last_encountered = a->identifier;
                std::stringstream ostr;
                ostr << s << "\t" << last_encountered << "\t" << traversed << "\n";
                outstrs.push_back(ostr.str());
                break;
            }
            //every node should be in assignments at this point.
            int anc_state = assignments.find(a->identifier)->second;
            if (anc_state == 0) {
                //encountered a 0, record the last encountered.
                std::stringstream ostr;
                ostr << s << "\t" << last_encountered << "\t" << traversed << "\n";
                outstrs.push_back(ostr.str());
                break;
            } else {
                last_encountered = a->identifier;
                traversed += a->mutations.size();
            }
        }
    }
    return outstrs;
}

//old code below
// std::vector<std::string> find_introductions(MAT::Tree* T, std::vector<std::string> samples, int k) {
//     /*
//     This function applies a five-step heuristic process to identify putative introductions to a population defined by the samples vector
//     1. The root is initialized as OUT.
//     2. Traversal of nodes in DFS order. Labels of leaf nodes don't change. For internal nodes, if there is a direct leaf descendant with 0 mutations to the internal node labelled IN, the internal node is labelled IN.
//     3. If there is no direct leaf descendant labelled as IN, the internal node is labelled IN only if its parent is labelled IN and the closest IN leaf node is at most k mutations from this internal node.
//     4. If 2 and 3 are not satisfied, the node is labelled as OUT.
//     5. Now all nodes are labelled IN or OUT. The number of introductions is simply the number of OUT->IN transitions in parent-child pairs which can be computed by traversing the tree again in a DFS order.
//     */
//     std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
//     std::map<std::string, int> assignments;
//     //the first real step is traversing nodes in DFS order.
//     auto dfs = T->depth_first_expansion();
//     for (auto n: dfs) {
//         if (n->is_root()) {
//             //rule 1
//             assignments[n->identifier] = 0;
//         } else if (n->is_leaf()) {
//             //rule 2. Leaves are in or out based on membership in samples only.
//             if (sample_set.find(n->identifier) != sample_set.end()) {
//                 assignments[n->identifier] = 1;
//             } else {
//                 assignments[n->identifier] = 0;
//             }
//         } else {
//             //check rule 2
//             bool r2 = false;
//             for (auto child: n->children) {
//                 if (child->is_leaf()) {
//                     if (sample_set.find(child->identifier) != sample_set.end() & (child->mutations.size() == 0)) {
//                         assignments[n->identifier] = 1;
//                         r2 = true;
//                     }
//                 }
//             }
//             //if not assigned via rule 2, check rule 3
//             bool r3 = false;
//             if (!r2) {
//                 auto sv = assignments.find(n->parent->identifier);
//                 if (sv != assignments.end()) {
//                     if (sv->second == 1) {
//                         r3 = true;
//                         assignments[n->identifier] = 1;
//                     } else {
//                         //now we get the most complex bit of logic
//                         //get terminals descended from this node, and for each one, if theyre members of population
//                         //then check to see if there are less than k total mutations on the rsearch
//                         auto leaves = T->get_leaves(n->identifier);
//                         for (auto l: leaves) {
//                             if (sample_set.find(l->identifier) != sample_set.end() & !r3) {
//                                 auto path = T->rsearch(l->identifier, true);
//                                 size_t mcount = 0;
//                                 for (auto p: path) {
//                                     if (p->identifier == n->identifier & (mcount <= static_cast<size_t>(k))) {
//                                         assignments[n->identifier] = 1;
//                                         r3 = true;
//                                         break;
//                                     } else if (mcount > k) {
//                                         break;
//                                     }
//                                     mcount += p->mutations.size();
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//             if (!r3 & !r2) {
//                 assignments[n->identifier] = 0;
//             }
//         }
//     }
//     //now that we have assignments in a map objects of node identifiers to states
//     //we can traverse again in dfs order and check where it changes from 0 to 1
//     //rule 4
//     std::vector<std::string> outstrs;
//     for (auto n: dfs) {
//         auto search = assignments.find(n->identifier);
//         assert (search != assignments.end());
//         bool has_ancestral = false;
//         if (search->second == 1) {
//             //check whether any ancestor of this node is IN
//             //we're only reporting the furthest back introductions on the tree.
//             size_t total_traversed = 0;
//             for (auto a: T->rsearch(n->identifier, true)) {
//                 if (a->identifier != n->identifier) {
//                     auto ancsearch = assignments.find(a->identifier);
//                     if ((ancsearch->second == 1) & (total_traversed <= static_cast<size_t>(k))) {
//                         //this had some IN in its history that is a relatively short mutational distance away
//                         //don't treat it as a novel event
//                         has_ancestral = true;
//                         break;
//                     }
//                 }
//                 total_traversed += a->mutations.size();
//             }
//             if (!has_ancestral) {
//                 //novel introduction
//                 //collect some info.
//                 //column 1 is the internal node identifier, column 2 is a comma-delineated list of terminals from
//                 //this internal node which are IN, and column 3 is a comma-delineated list of terminals which are OUT
//                 auto terminals = T->get_leaves(n->identifier);
//                 for (auto t: terminals) {
//                     if (assignments.find(t->identifier)->second == 1) {
//                         std::string ostr = t->identifier + "\t" + n->identifier + "\n";
//                         outstrs.push_back(ostr);
//                     } 
//                 }
//             }
//         }
//     }
//     return outstrs;
// }

void introduce_main(po::parsed_options parsed) {
    po::variables_map vm = parse_introduce_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples_filename = vm["population-samples"].as<std::string>();
    std::string output_file = vm["output"].as<std::string>();
    //uint32_t num_threads = vm["threads"].as<uint32_t>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    auto popsamples = read_sample_names(samples_filename);
    auto outstrings = find_introductions(&T, popsamples);

    std::ofstream of;
    of.open(output_file);
    of << "sample\tintroduction_node\tdistance\n";
    for (std::string o: outstrings) {
        of << o;
    }
    of.close();
}
