#include "optimize.hpp"
#include <time.h>  
#include <set>  

Timer timer;

po::variables_map check_options(int argc, char** argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::options_description desc("optimize options");
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file to optimize [REQUIRED].")
        ("output-mat,o", po::value<std::string>()->required(),
         "Output optimized mutation-annotated tree file [REQUIRED].") 
        ("sample-names,f", po::value<std::string>()->default_value(""),
         "File containing sample names whose ancestral nodes would be considered for optimization.") 
        ("radius,r", po::value<uint32_t>()->default_value(10), \
         "Radius in which to restrict the SPR moves.")
        ("optimization-seconds,s", po::value<uint32_t>()->default_value(3600), \
         "Approximate number of seconds to run the tree optimization stage. The stage terminates as soon as the elapsed time exceeds this value.")
        ("save-every-seconds,S", po::value<uint32_t>()->default_value(300), \
         "Periodically save the optimized tree after every specified number of seconds have elapsed since the last save.") 
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(p)
                .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless
        // the user specifies help
            if (vm.count("help"))
                exit(0);
            else
                exit(1);
    }
    return vm;
}

struct Pruned_Sample {
    std::string sample_name;
    std::vector<MAT::Mutation> sample_mutations;
    std::unordered_set<uint32_t> positions;

    // Assumes mutations are added in reverse chrono order
    void add_mutation (MAT::Mutation mut) {
        // If not reversal to reference allele 
        if ((mut.ref_nuc != mut.mut_nuc) && (positions.find(mut.position) == positions.end())) {
            auto iter = std::lower_bound(sample_mutations.begin(), sample_mutations.end(), mut);
            mut.par_nuc = mut.ref_nuc;
            sample_mutations.insert(iter, mut);
        }
        positions.insert(mut.position);
    }

    Pruned_Sample (std::string name) {
        sample_name = name;
        sample_mutations.clear();
        positions.clear();
    }
};


inline std::string get_reversal_mutation(const std::string& mut_string) {
    std::vector<std::string> words;
    MAT::string_split(mut_string, words);
    assert(words.size() == 4);
    return words[0] + "\t" + words[3] + "\t" + words[2] + "\t" + words[1]; 
}

size_t get_node_distance (const MAT::Tree& T, MAT::Node* source, MAT::Node* dest) {
    auto s_anc = T.rsearch(source->identifier, true);
    auto d_anc = T.rsearch(dest->identifier, true);

    size_t s_off=0, d_off=0, num_iter=0;
    if (s_anc.size() <= d_anc.size()) {
        d_off = d_anc.size() - s_anc.size();
        num_iter = s_anc.size();
    }
    else {
        s_off = s_anc.size() - d_anc.size();
        num_iter = d_anc.size();
    }
    
    auto distance = s_off+d_off;
    for (size_t i=0; i<num_iter; i++) {
        if (s_anc[s_off+i] == d_anc[d_off+i]) {
            break;
        }
        else {
            distance+=2;
        }
    }

    return distance;
}

/*
void optimize_main_old(po::parsed_options parsed) {
    po::variables_map vm = parse_optimize_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    uint32_t max_seconds = vm["optimization-seconds"].as<uint32_t>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    tbb::task_scheduler_init init(num_threads);
    srand (time(NULL));
    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "The parsimony score for this MAT is %zu\n", T.get_parsimony_score());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    timer.Start();
    fprintf(stderr, "Starting tree optimization.\n\n"); 
    
    fprintf(stderr, "Finding the nodes with recurrent or reversal mutations.\n");
    auto dfs = T.depth_first_expansion();
    
    std::map<std::string, int> mutation_counts;
    tbb::mutex tbb_lock;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
            [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
            auto n = dfs[i];
            n->clear_annotations();
            
            for (auto m: n->mutations) {
               std::string mut_string = m.chrom + "\t" + std::to_string(m.par_nuc) + "\t" + 
                   std::to_string(m.position) + "\t" + std::to_string(m.mut_nuc);
               tbb_lock.lock();
               if (mutation_counts.find(mut_string) == mutation_counts.end()) {
                   mutation_counts[mut_string] = 1;
               }
               else {
                   mutation_counts[mut_string] += 1;
               }
               tbb_lock.unlock();
           }
        }
    }, ap);
    std::vector<std::string> nodes_to_prune;
    for (auto n: dfs) {
        if ((n == T.root) || (T.get_num_leaves(n) > 1000)) {
            continue;
        }
        for (auto m: n->mutations) {
            std::string mut_string = m.chrom + "\t" + std::to_string(m.par_nuc) + "\t" + 
                std::to_string(m.position) + "\t" + std::to_string(m.mut_nuc);
            // Contains recurrent mutation
            if (mutation_counts[mut_string] > 1) {
                nodes_to_prune.emplace_back(n->identifier);
                break;
            }
            auto reversal_string = get_reversal_mutation(mut_string);
            // reversal mutation found
            if (mutation_counts.find(reversal_string) != mutation_counts.end()) {
                nodes_to_prune.emplace_back(n->identifier);
                break;
            }
        }
    }
    // Shuffle the nodes to be pruned
    std::random_shuffle(nodes_to_prune.begin(), nodes_to_prune.end());
    fprintf(stderr, "%zu nodes found with recurrent or reversal mutations.\n\n", nodes_to_prune.size());
    auto best_parsimony_score = T.get_parsimony_score();
    for (auto node_to_prune: nodes_to_prune) {
        if (T.get_node(node_to_prune) == NULL) {
            continue;
        }
        auto copy = MAT::get_tree_copy(T);
        fprintf(stderr, "Pruning subtree at node %s.\n", node_to_prune.c_str());
        auto leaves_to_prune = T.get_leaves(node_to_prune);
        std::vector<Pruned_Sample> pruned_samples;
        for (auto l: leaves_to_prune) {
            Pruned_Sample to_prune(l->identifier);
            auto root_to_node = T.rsearch(l->identifier, true); 
            std::reverse(root_to_node.begin(), root_to_node.end());
            for (auto curr: root_to_node) {
                for (auto m: curr->mutations) {
                    to_prune.add_mutation(m);
                }
            }
            // move_level set to false as moving it would be an unnecessary overhead
            T.remove_node(l->identifier, false);
            pruned_samples.emplace_back(to_prune);
        }
        fprintf(stderr, "Pruned %zu samples at node %s.\n", pruned_samples.size(), node_to_prune.c_str());
        fprintf(stderr, "Placing %zu samples at node %s in random order.\n", pruned_samples.size(), node_to_prune.c_str());
        std::random_shuffle(pruned_samples.begin(), pruned_samples.end());
        for (auto s: pruned_samples) {
            dfs = T.depth_first_expansion();
            size_t total_nodes = dfs.size();
            // Stores the excess mutations to place the sample at each
            // node of the tree in DFS order. When placement is as a
            // child, it only contains parsimony-increasing mutations in
            // the sample. When placement is as a sibling, it contains 
            // parsimony-increasing mutations as well as the mutations
            // on the placed node in common with the new sample. Note
            // guaranteed to be corrrect only for optimal nodes since
            // the mapper can terminate the search early for non-optimal
            // nodes
            std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
            // Stores the imputed mutations for ambiguous bases in the
            // sampled in order to place the sample at each node of the 
            // tree in DFS order. Again, guaranteed to be corrrect only 
            // for pasrimony-optimal nodes 
            std::vector<std::vector<MAT::Mutation>> node_imputed_mutations(total_nodes);
            std::vector<int> node_set_difference;
            size_t best_node_num_leaves = 0;
            // The maximum number of mutations is bound by the number
            // of mutations in the missing sample (place at root)
            //int best_set_difference = 1e9;
            // TODO: currently number of root mutations is also added to
            // this value since it forces placement as child but this
            // could be changed later 
            int best_set_difference = s.sample_mutations.size() + T.root->mutations.size() + 1;
            size_t best_j = 0;
            bool best_node_has_unique = false;
            std::vector<bool> node_has_unique(total_nodes, false);
            std::vector<size_t> best_j_vec;
            size_t num_best = 1;
            MAT::Node* best_node = T.root;
            best_j_vec.emplace_back(0);
            // Parallel for loop to search for most parsimonious
            // placements. Real action happens within mapper2_body
            static tbb::affinity_partitioner ap;
            tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k){
                    mapper2_input inp;
                    inp.T = &T;
                    inp.node = dfs[k];
                    inp.missing_sample_mutations = &s.sample_mutations;
                    inp.excess_mutations = &node_excess_mutations[k];
                    inp.imputed_mutations = &node_imputed_mutations[k];
                    inp.best_node_num_leaves = &best_node_num_leaves;
                    inp.best_set_difference = &best_set_difference;
                    inp.best_node = &best_node;
                    inp.best_j =  &best_j;
                    inp.num_best = &num_best;
                    inp.j = k;
                    inp.has_unique = &best_node_has_unique;
                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);
                    mapper2_body(inp, false);
                    }       
                    }, ap); 
            // Ensure sample not already in the tree
            if (T.get_node(s.sample_name) == NULL) {
                // Is placement as sibling
                if (best_node->is_leaf() || best_node_has_unique) {
                    std::string nid = std::to_string(++T.curr_internal_node);
                    T.create_node(nid, best_node->parent->identifier);
                    T.create_node(s.sample_name, nid);
                    T.move_node(best_node->identifier, nid);
                    // common_mut stores mutations common to the
                    // best node branch and the sample, l1_mut
                    // stores mutations unique to best node branch
                    // and l2_mut stores mutations unique to the
                    // sample not in best node branch
                    std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                    std::vector<MAT::Mutation> curr_l1_mut;
                    // Compute current best node branch mutations
                    for (auto m1: best_node->mutations) {
                        MAT::Mutation m = m1.copy();
                        curr_l1_mut.emplace_back(m);
                    }
                    // Clear mutations on the best node branch which
                    // will be later replaced by l1_mut
                    best_node->clear_mutations();
                    // Compute l1_mut
                    for (auto m1: curr_l1_mut) {
                        bool found = false;
                        for (auto m2: node_excess_mutations[best_j]) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            MAT::Mutation m = m1.copy();
                            l1_mut.emplace_back(m);
                        }
                    }
                    // Compute l2_mut
                    for (auto m1: node_excess_mutations[best_j]) {
                        bool found = false;
                        for (auto m2: curr_l1_mut) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    MAT::Mutation m = m1.copy();
                                    common_mut.emplace_back(m);
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            MAT::Mutation m = m1.copy();
                            l2_mut.emplace_back(m);
                        }
                    }
                    // Add mutations to new node using common_mut
                    for (auto m: common_mut) {
                        T.get_node(nid)->add_mutation(m);
                    }
                    // Add mutations to best node using l1_mut
                    for (auto m: l1_mut) {
                        T.get_node(best_node->identifier)->add_mutation(m);
                    }
                    // Add new sample mutations using l2_mut
                    for (auto m: l2_mut) {
                        T.get_node(s.sample_name)->add_mutation(m);
                    }
                }
                // Else placement as child
                else {
                    MAT::Node* node = T.create_node(s.sample_name, best_node->identifier);
                    std::vector<MAT::Mutation> node_mut;
                    std::vector<MAT::Mutation> curr_l1_mut;
                    for (auto m1: best_node->mutations) {
                        MAT::Mutation m = m1.copy();
                        curr_l1_mut.emplace_back(m);
                    }
                    for (auto m1: node_excess_mutations[best_j]) {
                        bool found = false;
                        for (auto m2: curr_l1_mut) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            MAT::Mutation m = m1.copy();
                            node_mut.emplace_back(m);
                        }
                    }
                    for (auto m: node_mut) {
                        node->add_mutation(m);
                    }
                }
            }
        }
        fprintf(stderr, "Done pruning and placing samples at node %s\n", node_to_prune.c_str());
        auto new_parsimony_score = T.get_parsimony_score();
        if (new_parsimony_score >= best_parsimony_score) {
            fprintf(stderr, "Placement resulted in a parsimony score of %zu. Retaining old tree.\n", new_parsimony_score);
            MAT::clear_tree(T);
            T = copy;
        }
        else {
            MAT::clear_tree(copy);
            best_parsimony_score = new_parsimony_score; 
            fprintf(stderr, "Placement lowered parsimony score to %zu!\n", best_parsimony_score);
        }
        auto elapsed_time = timer.Stop()/1000;
        fprintf(stderr, "Elapsed time (tree optimization stage): %ld seconds.\n\n", elapsed_time);
        if (elapsed_time >= max_seconds) {
            break;
        }
    }
    
    fprintf(stderr, "Optimization complete. Saving the final MAT with a parsimony score of %zu\n", T.get_parsimony_score());
    T.collapse_tree();
    T.condense_leaves();
    MAT::save_mutation_annotated_tree(T, output_mat_filename);
    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
*/

int main(int argc, char** argv) {
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string sample_name_filename = vm["sample-names"].as<std::string>();
    uint32_t max_seconds = vm["optimization-seconds"].as<uint32_t>();
    uint32_t save_every = vm["save-every-seconds"].as<uint32_t>();
    uint32_t radius = vm["radius"].as<uint32_t>();

    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);
    srand (time(NULL));

    static tbb::affinity_partitioner ap;

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    // Collapse tree for optimal performance and results
    T.collapse_tree();
    
    fprintf(stderr, "The parsimony score for this MAT is %zu\n", T.get_parsimony_score());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    float last_update = 0;
    
    timer.Start();
    fprintf(stderr, "Finding the nodes to prune.\n");
    auto bfs = T.breadth_first_expansion();
    
    std::set<std::string> nodes_to_prune;
    
    if (sample_name_filename != "") {
        std::ifstream infile(sample_name_filename);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the sample name file: %s!\n", sample_name_filename.c_str());
            exit(1);
        }    
        std::string line;

        fprintf(stderr, "Reading sample name file.\n"); 
        timer.Start();
        while (std::getline(infile, line)) {
            std::vector<std::string> words;
            MAT::string_split(line, words);
            if (words.size() != 1) {
                fprintf(stderr, "ERROR: Incorrect format for sample name file: %s!\n", sample_name_filename.c_str());
                exit(1);
            }
            auto n = T.get_node(words[0]);
            if (n == NULL) {
                fprintf(stderr, "WARNING: Node id %s not found in the tree!\n", words[0].c_str());
                exit(1);
            }
            else {
                for (auto anc: T.rsearch(words[0], true)) {
                    if (anc->level >= 4) { 
                        nodes_to_prune.insert(anc->identifier);
                    }
                }
            }
        }
        infile.close();
    }

    else {
        for (auto n: bfs) {
            if ((n == T.root) || (n->level < 4)) {
                continue;
            }
            size_t count=0;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, bfs.size()),
                    [&](const tbb::blocked_range<size_t> r) {
                    for (size_t i=r.begin(); i<r.end(); ++i){
                    auto n2 = bfs[i];
                    if (n == n2) {
                    continue;
                    }
                    for (auto m1: n->mutations) {
                    for (auto m2: n2->mutations) {
                    if  (m1.position == m2.position)  {
                    if (get_node_distance(T, n, n2) <= radius) {
                    __sync_fetch_and_add(&count, 1);
                    }
                    }
                    }
                    }
                    }
                    }, ap);

            if (count > 0) {
                nodes_to_prune.insert(n->identifier);
            }
        }
    }

    fprintf(stderr, "%zu nodes found to prune.\n", nodes_to_prune.size());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    
    timer.Start();
    fprintf(stderr, "Starting tree optimization.\n\n"); 

    auto best_parsimony_score = T.get_parsimony_score();
    for (auto nid_to_prune: nodes_to_prune) {
        auto node_to_prune = T.get_node(nid_to_prune);
        if (node_to_prune == NULL) {
            continue;
        }
        
        fprintf(stderr, "Pruning subtree at node %s.\n", nid_to_prune.c_str());

        // Find mutations on the node to prune
        Pruned_Sample pruned_sample(nid_to_prune);

        auto node_to_root = T.rsearch(nid_to_prune, true); 
        //std::reverse(root_to_node.begin(), root_to_node.end());
        for (auto curr: node_to_root) {
            for (auto m: curr->mutations) {
                pruned_sample.add_mutation(m);
            }
        }

        // Now prune the node_to_prune from the tree
        auto curr_parent = node_to_prune->parent; 
        auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), node_to_prune);
        assert (iter != curr_parent->children.end());
        curr_parent->children.erase(iter);

        // Set source to current parent
        auto source = curr_parent;

        // Remove curr_parent if it has no children
        if (curr_parent->children.size() == 0) {
            auto ancestors = T.rsearch(curr_parent->identifier, true);
            T.remove_node(curr_parent->identifier, true);

            // Since remove_node can remove multiple levels of ancestors,
            // reassign source to nearest ancestor currently in the tree
            for (auto anc: ancestors) {
                if (T.get_node(anc->identifier) != NULL) {
                    source = anc;
                    break;
                }
            }
        }

        // Move the remaining child one level up if it is the only child of its parent 
        if (curr_parent->children.size() == 1) {
            auto child = curr_parent->children[0];
            if (curr_parent->parent != NULL) {
                child->parent = curr_parent->parent;
                child->level = curr_parent->parent->level + 1;

                std::vector<MAT::Mutation> tmp;
                for (auto m: child->mutations) {
                    tmp.emplace_back(m.copy());
                }

                //Clear and add back mutations in chrono order
                child->clear_mutations();
                for (auto m: curr_parent->mutations) {
                    child->add_mutation(m.copy());
                }
                for (auto m: tmp) {
                    child->add_mutation(m.copy());
                }

                curr_parent->parent->children.push_back(child);

                curr_parent->children.clear();
                T.remove_node(curr_parent->identifier, false);

                source = child;
            }
        }

        bfs = T.breadth_first_expansion();
        size_t total_nodes = bfs.size();

        std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
        std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);

        std::vector<int> node_set_difference;

        size_t best_node_num_leaves = 0;
        int best_set_difference = 1e9;

        std::vector<bool> node_has_unique(total_nodes);
        size_t best_j = 0;
        bool best_node_has_unique = false;

        size_t best_distance = 1e9;

        std::vector<size_t> best_j_vec;

        size_t num_best = 1;
        MAT::Node* best_node = T.root;
        best_j_vec.emplace_back(0);

        std::vector<size_t> node_distance(total_nodes);

        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                [&](tbb::blocked_range<size_t> r) {
                for (size_t k=r.begin(); k<r.end(); ++k){
                node_distance[k] = get_node_distance(T, source, bfs[k]);
                if (node_distance[k] > radius) {
                   continue;
                }
                
                node_has_unique[k] = false;

                mapper2_input inp;
                inp.T = &T;
                inp.node = bfs[k];
                inp.missing_sample_mutations = &pruned_sample.sample_mutations;
                inp.excess_mutations = &node_excess_mutations[k];
                inp.imputed_mutations = &imputed_mutations[k];
                inp.best_node_num_leaves = &best_node_num_leaves;
                inp.best_set_difference = &best_set_difference;
                inp.best_node = &best_node;
                inp.best_j =  &best_j;
                inp.num_best = &num_best;
                inp.j = k;
                inp.has_unique = &best_node_has_unique;

                inp.distance = node_distance[k];
                inp.best_distance = &best_distance;

                inp.best_j_vec = &best_j_vec;
                inp.node_has_unique = &(node_has_unique);

                mapper2_body(inp, false);
                }       
                }, ap); 

        auto distance = node_distance[best_j]; 

        // Clear current mutations of the node to prune
        node_to_prune->clear_mutations();

        // Is placement as sibling
        if (best_node->is_leaf() || best_node_has_unique) {
            std::string nid = std::to_string(++T.curr_internal_node);
            auto new_node = T.create_node(nid, best_node->parent->identifier);
            node_to_prune->parent = new_node;
            new_node->children.emplace_back(node_to_prune);
            T.move_node(best_node->identifier, nid);
            
            std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
            std::vector<MAT::Mutation> curr_l1_mut;

            // Compute current best node branch mutations
            for (auto m1: best_node->mutations) {
                MAT::Mutation m = m1.copy();
                curr_l1_mut.emplace_back(m);
            }
            // Clear mutations on the best node branch which
            // will be later replaced by l1_mut
            best_node->clear_mutations();

            // Compute l1_mut
            for (auto m1: curr_l1_mut) {
                bool found = false;
                for (auto m2: node_excess_mutations[best_j]) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if (m1.position == m2.position) {
                        if (m1.mut_nuc == m2.mut_nuc) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    l1_mut.emplace_back(m);
                }
            }
            // Compute l2_mut
            for (auto m1: node_excess_mutations[best_j]) {
                bool found = false;
                for (auto m2: curr_l1_mut) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if (m1.position == m2.position) {
                        if (m1.mut_nuc == m2.mut_nuc) {
                            found = true;
                            MAT::Mutation m = m1.copy();
                            common_mut.emplace_back(m);
                            break;
                        }
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    l2_mut.emplace_back(m);
                }
            }

            // Add mutations to new node using common_mut
            for (auto m: common_mut) {
                new_node->add_mutation(m);
            }
            // Add mutations to best node using l1_mut
            for (auto m: l1_mut) {
                best_node->add_mutation(m);
            }
            // Add new sample mutations using l2_mut
            for (auto m: l2_mut) {
                node_to_prune->add_mutation(m);
            }
        }
        // Else placement as child
        else {
            best_node->children.emplace_back(node_to_prune);
            node_to_prune->parent = best_node;

            std::vector<MAT::Mutation> node_mut;

            std::vector<MAT::Mutation> curr_l1_mut;

            for (auto m1: best_node->mutations) {
                MAT::Mutation m = m1.copy();
                curr_l1_mut.emplace_back(m);
            }

            for (auto m1: node_excess_mutations[best_j]) {
                bool found = false;
                for (auto m2: curr_l1_mut) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if (m1.position == m2.position) {
                        if (m1.mut_nuc == m2.mut_nuc) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    node_mut.emplace_back(m);
                }
            }
            for (auto m: node_mut) {
                node_to_prune->add_mutation(m);
            }
        }
        
        fprintf(stderr, "Done pruning and placing samples at node %s\n", best_node->identifier.c_str());
        fprintf(stderr, "Distance from original placement: %zu\n", distance);

        auto new_parsimony_score = T.get_parsimony_score();
        assert(new_parsimony_score <= best_parsimony_score);
        if (new_parsimony_score == best_parsimony_score) {
            fprintf(stderr, "Parsimony score of %zu unchanged after placement.\n", new_parsimony_score);
        }
        else {
            best_parsimony_score = new_parsimony_score; 
            fprintf(stderr, "Placement lowered parsimony score to %zu!\n", best_parsimony_score);
        }
        
        float elapsed_sec = static_cast<float>(timer.Stop())/1000;
        fprintf(stderr, "Elapsed time (tree optimization stage): %.3f sec\n\n", elapsed_sec);

        if (elapsed_sec >= max_seconds) {
            break;
        }
        
        if (elapsed_sec-last_update >= save_every) {
            // Levels need to be adjusted after placement 
            bfs = T.depth_first_expansion();
            for (auto n: bfs) {
                if (n == T.root) {
                    n->level = 1;
                }
                else {
                    n->level = n->parent->level+1;
                }
            }

            fprintf(stderr, "Saving the current MAT with a parsimony score of %zu\n\n", T.get_parsimony_score());
            T.condense_leaves();
            MAT::save_mutation_annotated_tree(T, output_mat_filename);
            last_update = elapsed_sec;

            // Uncondense for the next iteration
            T.uncondense_leaves();
        }
    }
    
    // Levels need to be adjusted after placement 
    bfs = T.breadth_first_expansion();
    for (auto n: bfs) {
        if (n == T.root) {
            n->level = 1;
        }
        else {
            n->level = n->parent->level+1;
        }
    }

    fprintf(stderr, "Optimization complete. Saving the final MAT with a parsimony score of %zu\n", T.get_parsimony_score());
    T.condense_leaves();
    MAT::save_mutation_annotated_tree(T, output_mat_filename);
    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

