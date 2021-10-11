#include "merge.hpp"

po::variables_map parse_merge_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description merge_desc("merge options");
    merge_desc.add_options()("input-mat-1", po::value<std::string>()->required(),
                             "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
    ("input-mat-2", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
    ("output-mat,o", po::value<std::string>()->required(),
     "Write output files to the target directory. Default is current directory.")
    ("max-depth,d", po::value<uint32_t>()->default_value(20),
     "Max depth to consider in the subtree rooted at the consistent node.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str());

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(merge_desc)
                  .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << merge_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

std::string get_first_leaf(MAT::Node* n) {
    auto ret = n;
    while (ret->children.size() > 0) {
        ret = ret->children[0];
    }
    return ret->identifier;
}

/**
 * Checks for consistency between both MAT files to ensure that they are able to merge
 **/
bool consistent(MAT::Tree A, MAT::Tree B, concurMap& consistNodes) {
    //vectors of all leaves in both input trees
    std::vector<std::string> A_leaves = A.get_leaves_ids();
    std::vector<std::string> B_leaves = B.get_leaves_ids();
    tbb::parallel_sort(A_leaves.begin(), A_leaves.end());
    tbb::parallel_sort(B_leaves.begin(), B_leaves.end());

    //creates a vector of common_leaves between two input trees
    std::vector<std::string> common_leaves;
    set_intersection(B_leaves.begin(), B_leaves.end(), A_leaves.begin(), A_leaves.end(), std::back_inserter(common_leaves));
    tbb::parallel_sort(common_leaves.begin(), common_leaves.end());
    fprintf(stderr, "%zu common leaves.\n", common_leaves.size());


    if (common_leaves.size() == 0) {
        return true;
    }

    //creates two subtrees using the common_leaves
    std::vector<std::string> A_to_remove, B_to_remove;

    auto Asub = get_tree_copy(A);
    std::set_difference(A_leaves.begin(), A_leaves.end(), common_leaves.begin(), common_leaves.end(), std::back_inserter(A_to_remove));

    for (auto l : A_to_remove) {
        Asub.remove_node(l, false);
    }
    Asub.remove_single_child_nodes();

//    auto Bsub = get_tree_copy(B);
//    std::set_difference(B_leaves.begin(), B_leaves.end(), common_leaves.begin(), common_leaves.end(), std::back_inserter(B_to_remove));
//
//    for (auto l : B_to_remove) {
//        Bsub.remove_node(l, false);
//    }
//    Bsub.remove_single_child_nodes();


    auto Adfs = Asub.depth_first_expansion();
//    auto Bdfs = Bsub.depth_first_expansion();
//    if (Adfs.size() != Bdfs.size()) {
//        fprintf(stderr, "ERROR: Different DFS sizes for the common subtrees (%zu and %zu). MATs may not be consistent.",
//                Adfs.size(), Bdfs.size());
//    }
//    Asub.rotate_for_consistency();
//    Bsub.rotate_for_consistency();

    bool ret = true;
    static tbb::affinity_partitioner ap;
    /**
     * Parallel loop that parses through the depth first expansion of
     * both subtrees and compares each node's mutations
     **/
    tbb::mutex m;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, Adfs.size()),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t k = r.begin(); k < r.end(); ++k) {

            auto n = Adfs[k];
            if (n->children.size() > 1) {
                auto c1 = n->children[0];
                auto c2 = n->children[1];
                auto l1 = get_first_leaf(c1);
                auto l2 = get_first_leaf(c2);
                auto lca1 = MAT::LCA(A, l1, l2);
                auto lca2 = MAT::LCA(B, l1, l2);

                if ((lca1 != NULL) && (lca2!= NULL)) {
                    consistNodes.emplace(std::pair<std::string, std::string> (lca2->identifier, lca1->identifier));
                } else {
                    fprintf(stderr, "NULL found!\n");
                }
            }

            if (n->children.size() == 0) {
                consistNodes.emplace(std::pair<std::string, std::string> (n->identifier, n->identifier));
            }
        }
    },
    ap);

    if (consistNodes.size() != Adfs.size()) {
        fprintf (stderr, "WARNING: MATs not completely consistent!\n");
    }
    fprintf (stderr, "%zu of %zu nodes consistent.\n", consistNodes.size(), Adfs.size());
    return ret;
}

void merge_main(po::parsed_options parsed) {
    timer.Start();
    //Takes user input and loads the specified MAT files
    po::variables_map vm = parse_merge_command(parsed);
    std::string mat1_filename = vm["input-mat-1"].as<std::string>();
    std::string mat2_filename = vm["input-mat-2"].as<std::string>();
    std::string output_filename = vm["output-mat"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    uint32_t max_levels = vm["max-depth"].as<uint32_t>();


    fprintf(stderr, "Initializing %u worker threads.\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);

    fprintf(stderr, "Loading first input MAT. Existing clade annotations will be cleared\n");
    MAT::Tree mat1 = MAT::load_mutation_annotated_tree(mat1_filename);
    for (auto n: mat1.depth_first_expansion()) {
        n->clear_annotations();
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Loading second input MAT. Existing clade annotations will be cleared\n");
    MAT::Tree mat2 = MAT::load_mutation_annotated_tree(mat2_filename);
    for (auto n: mat2.depth_first_expansion()) {
        n->clear_annotations();
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Checking MAT consistency.\n");
    MAT::Tree baseMat;
    MAT::Tree otherMat;
    concurMap consistNodes;

    //uncondenses nodes of both MAT files
    if (mat1.condensed_nodes.size() > 0) {
        mat1.uncondense_leaves();
    }
    if (mat2.condensed_nodes.size() > 0) {
        mat2.uncondense_leaves();
    }
    //Assigns largest MAT to baseMat and smaller one to otherMat
    if (mat1.get_num_leaves() > mat2.get_num_leaves()) {
        baseMat = mat1;
        otherMat = mat2;
    } else {
        baseMat = mat2;
        otherMat = mat1;
    }

    //Checks for consistency in mutation paths between the two trees
    //if (consistent(baseMat, otherMat) == false) {
    //fprintf(stderr, "WARNING: MAT files are not consistent!\n");
    //exit(1);
    //}
    consistent(baseMat, otherMat, consistNodes);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Finding new samples to merge\n");
    //Makes a copy of the baseMat files to add the samples to later
    MAT::Tree finalMat = get_tree_copy(baseMat);
    std::vector<std::string> samples;

    auto otherLeaves = otherMat.get_leaves_ids();
    auto baseLeaves = baseMat.get_leaves_ids();

    tbb::parallel_sort(otherLeaves.begin(), otherLeaves.end());
    tbb::parallel_sort(baseLeaves.begin(), baseLeaves.end());

    //creates vector of new samples to be added
    std::set_difference(otherLeaves.begin(), otherLeaves.end(), baseLeaves.begin(), baseLeaves.end(), std::back_inserter(samples));
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();

    static tbb::affinity_partitioner ap;
    fprintf(stderr, "Merging %lu new samples\n", samples.size());

    tbb::mutex m;
    size_t num_mutex = 8*num_threads;
    std::vector<tbb::mutex> m_vec(num_mutex);

    std::map<std::string, size_t> node_to_first_sample;

    //parallel loop that concurrently parses through all the new samples
    int i = 0;
    for (size_t k = 0; k < samples.size(); k++) {

        std::string x = samples[k];
        auto ancestors = otherMat.rsearch(x, true);
        std::string curr = finalMat.root->identifier;

        /**
         * Parses through the the ancestors of each sample on otherMat and finds
         * Closest consistent node on baseMat
         **/
        for (auto anc: ancestors) {
            if (consistNodes.find(anc->identifier) != consistNodes.end()) {
                curr = consistNodes[anc->identifier];
                break;
            }
        }

        MAT::Node s(x, -1);
        s.mutations.clear();
        for (int y = ancestors.size()-1; y >= 0; y--) {
            for (auto m: ancestors[y]->mutations) {
                s.add_mutation(m);
            }
        }
        std::sort(s.mutations.begin(), s.mutations.end());

        //Roots bfs at closest consistent node identified in previous loop
        //Restricts tree search to a smaller subtree
        auto bfs = finalMat.breadth_first_expansion(curr);

        std::vector<std::vector<MAT::Mutation> > node_excess_mutations(bfs.size());
        std::vector<std::vector<MAT::Mutation> > node_imputed_mutations(bfs.size());
        size_t best_node_num_leaves = 0;
        int best_set_difference = s.mutations.size() + bfs[0]->mutations.size() + 1;
        size_t best_j = 0;
        size_t num_best = 1;
        bool best_node_has_unique = false;
        MAT::Node *best_node = bfs[0];
        std::vector<bool> node_has_unique(bfs.size(), false);
        std::vector<size_t> best_j_vec;
        best_j_vec.emplace_back(0);


        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, bfs.size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t k = r.begin(); k < r.end(); k++) {
                if (bfs[k]->level - bfs[0]->level > max_levels) {
                    continue;
                }

                mapper2_input inp;

                inp.T = &finalMat;
                inp.node = bfs[k];
                inp.missing_sample_mutations = &s.mutations;
                inp.excess_mutations = &node_excess_mutations[k];
                inp.imputed_mutations = &node_imputed_mutations[k];

                inp.best_node_num_leaves = &best_node_num_leaves;
                inp.best_set_difference = &best_set_difference;
                inp.best_node = &best_node;
                inp.best_j = &best_j;
                inp.num_best = &num_best;
                inp.j = k;
                inp.has_unique = &best_node_has_unique;
                inp.best_j_vec = &best_j_vec;
                inp.node_has_unique = &(node_has_unique);
                mapper2_body(inp, false);
            }
        }, ap);

        if (finalMat.get_node(x) == NULL) {
            // Is placement as sibling
            if (best_node->is_leaf() || best_node_has_unique) {
                m.lock();
                std::string nid = finalMat.new_internal_node_id();
                finalMat.create_node(nid, best_node->parent->identifier);
                finalMat.create_node(x, nid);
                finalMat.move_node(best_node->identifier, nid);
                m.unlock();
                // common_mut stores mutations common to the
                // best node branch and the sample, l1_mut
                // stores mutations unique to best node branch
                // and l2_mut stores mutations unique to the
                // sample not in best node branch
                std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                std::vector<MAT::Mutation> curr_l1_mut;

                // Compute current best node branch mutations
                for (auto m1 : best_node->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }
                // Clear mutations on the best node branch which
                // will be later replaced by l1_mut
                best_node->clear_mutations();

                // Compute l1_mut
                for (auto m1 : curr_l1_mut) {
                    bool found = false;
                    for (auto m2 : node_excess_mutations[best_j]) {
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
                for (auto m1 : node_excess_mutations[best_j]) {
                    bool found = false;
                    for (auto m2 : curr_l1_mut) {
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
                for (auto m : common_mut) {
                    finalMat.get_node(nid)->add_mutation(m);
                }
                // Add mutations to best node using l1_mut
                for (auto m : l1_mut) {
                    finalMat.get_node(best_node->identifier)->add_mutation(m);
                }
                // Add new sample mutations using l2_mut

                for (auto m : l2_mut) {
                    finalMat.get_node(x)->add_mutation(m);
                }
            }
            // Else placement as child
            else {
                m.lock();
                finalMat.create_node(x, best_node->identifier);
                MAT::Node *node = finalMat.get_node(x);
                m.unlock();
                std::vector<MAT::Mutation> node_mut;

                std::vector<MAT::Mutation> curr_l1_mut;

                for (auto m1 : best_node->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }

                for (auto m1 : node_excess_mutations[best_j]) {
                    bool found = false;
                    for (auto m2 : curr_l1_mut) {
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

                for (auto m : node_mut) {

                    node->add_mutation(m);
                }
            }
        }

        fprintf(stderr, "\rAdded %d of %zu samples", ++i, samples.size());
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    //Save final MAT file
    timer.Start();
    fprintf(stderr, "Condensing and saving final MAT\n");
    finalMat.condense_leaves();
    MAT::save_mutation_annotated_tree(finalMat, output_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

