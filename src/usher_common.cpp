#include "usher_common.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

//returns exit code
int usher_common(std::string dout_filename, std::string outdir, uint32_t max_trees,
                 uint32_t max_uncertainty, uint32_t max_parsimony, bool sort_before_placement_1, bool sort_before_placement_2, bool sort_before_placement_3,
                 bool reverse_sort, bool collapse_tree, bool collapse_output_tree, bool print_uncondensed_tree, bool print_parsimony_scores,
                 bool retain_original_branch_len, bool no_add, bool detailed_clades, size_t print_subtrees_size, size_t print_subtrees_single,
                 std::vector<Missing_Sample>& missing_samples, std::vector<std::string>& low_confidence_samples, MAT::Tree* loaded_MAT) {


    if (print_subtrees_size == 1) {
        std::cerr << "ERROR: print-subtrees-size should be larger than 1\n";
        return 1;
    }

    if (sort_before_placement_1 + sort_before_placement_2 + sort_before_placement_3 > 1) {
        std::cerr << "ERROR: Can't use two or more of sort-before-placement-1, sort-before-placement-2 and sort-before-placement-3 simultaneously. Please specify only one.\n";
        return 1;
    }

    if (sort_before_placement_1 || sort_before_placement_2 || sort_before_placement_3) {
        std::cerr << "WARNING: Using experimental option ";
        if (sort_before_placement_1) {
            std::cerr << "--sort-before-placement-1 (-s)\n";
        }
        if (sort_before_placement_2) {
            std::cerr << "--sort-before-placement-2 (-S)\n";
        }
        if (sort_before_placement_3) {
            std::cerr << "--sort-before-placement-3 (-A)\n";
        }
    } else if (reverse_sort) {
        std::cerr << "ERROR: Can't use reverse-sort without sorting options (sort-before-placement-1 or sort-before-placement-2 or sort-before-placement-3)\n";
        return 1;
    }

    if (print_parsimony_scores) {
        if (max_trees > 1) {
            std::cerr << "ERROR: cannot use --multiple-placements (-M) and --print_parsimony_scores (-p) options simulaneously.\n";
            return 1;
        }

        if (sort_before_placement_1 || sort_before_placement_2 || sort_before_placement_3 ||
                collapse_tree || collapse_output_tree || print_uncondensed_tree || (print_subtrees_size > 0) || (dout_filename != "")) {
            fprintf (stderr, "WARNING: --print-parsimony-scores-per-node is set. Will terminate without modifying the original tree.\n");
        }
    }

    if (max_trees == 0) {
        std::cerr << "ERROR: Number of trees specified by --multiple-placements (-M) should be >= 1\n";
        return 1;
    }
    /*
    if (max_trees > 1) {
        if (max_trees < 256) {
            std::cerr << "WARNING: Using experimental option --multiple-placements (-M)\n";
        }
        else {
            std::cerr << "ERROR: Number of trees specified by --multiple-placements (-M) should be <= 255\n";
            return 1;
        }
    }
    */

    if (no_add && (print_subtrees_size > 0 || print_subtrees_single)) {
        std::cerr << "ERROR: Sorry, cannot output subtrees when -n/--no-add is specified.\n";
        return 1;
    }

    if (retain_original_branch_len) {
        fprintf(stderr, "Output newick files will retain branch lengths from the input tree (unspecified at branches modified during the placement).\n\n");
    } else {
        fprintf(stderr, "Output newick files will have branch lengths equal to the number of mutations of that branch.\n\n");
    }

    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();

    // timer object to be used to measure runtimes of individual stages
    Timer timer;

//    fprintf(stderr, "Initializing %u worker threads.\n\n", num_threads);
//    tbb::task_scheduler_init init(num_threads);

#if SAVE_PROFILE == 1
    Instrumentor::Get().BeginSession("test-main", "p1.json");
#endif

    // Vector to store multiple trees, each corresponding to a different
    // possibility of a  parsimony-optimal placement, when --multiple-placements
    // is used. Otherwise, this vector maintains a single tree througout the
    // execution in which a tie-breaking strategy defined in usher_mapper is
    // used for multiple parsimony-optimal placements.
    std::vector<MAT::Tree> optimal_trees;



    // Tree pointer to point to some element in optimal_trees that would be
    // updated several times during the execution

    optimal_trees.emplace_back(std::move(*loaded_MAT));
    MAT::Tree* T = &optimal_trees[0];
    // Since --multiple-placements can result in trees with different parsimony
    // scores, the vector below will be used to maintain the final parsimony
    // score of each tree
    std::vector<size_t> tree_parsimony_scores;

    auto num_trees = optimal_trees.size();

    // Collapses the tree nodes not carrying a mutation and also condenses
    // identical sequences into a single node.
    if (collapse_tree) {
        timer.Start();

        fprintf(stderr, "Collapsing input tree.\n");

        assert (optimal_trees.size() == 1);

        T->collapse_tree();

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        fprintf(stderr, "Condensing identical sequences. \n");

        T->condense_leaves();

        auto condensed_tree_filename = outdir + "/condensed-tree.nh";
        fprintf(stderr, "Writing condensed input tree to file %s\n", condensed_tree_filename.c_str());

        FILE* condensed_tree_file = fopen(condensed_tree_filename.c_str(), "w");
        fprintf(condensed_tree_file, "%s\n", MAT::get_newick_string(*T, T->root, true, true, retain_original_branch_len).c_str());
        fclose(condensed_tree_file);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    fprintf(stderr, "Found %zu missing samples.\n\n", missing_samples.size());

    FILE* parsimony_scores_file = NULL;

    //Sort samples based on number of ambiguous bases if specified
    if (sort_before_placement_3) {
        timer.Start();
        fprintf(stderr, "Sorting missing samples based on the number of ambiguous bases \n");
        std::stable_sort(missing_samples.begin(), missing_samples.end());
        // Reverse sorted order if specified
        if (reverse_sort) {
            std::reverse(missing_samples.begin(), missing_samples.end());
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    // If samples found in VCF that are missing from the input tree, they are
    // now placed using maximum parsimony or if print_parsimony_score is set,
    // the parsimony scores for placning the sample is printed (without the
    // actual placement).
    if (missing_samples.size() > 0) {

        // indexes stores the order in which the missing samples should be
        // placed sequentially on the tree. It is initialized with sequentially
        // increasing values 0,1,2,.. which will be modified if one of the
        // sorting options is used
        std::vector<size_t> indexes(missing_samples.size());
        std::iota(indexes.begin(), indexes.end(), 0);

        // Write current tree with internal nodes labelled. Parsimony score of
        // placement at each node (including internal nodes) will be printed later.
        if (print_parsimony_scores) {
            timer.Start();
            auto current_tree_filename = outdir + "/current-tree.nh";

            fprintf(stderr, "Writing current tree with internal nodes labelled to file %s \n", current_tree_filename.c_str());
            FILE* current_tree_file = fopen(current_tree_filename.c_str(), "w");
            fprintf(current_tree_file, "%s\n", MAT::get_newick_string(*T, true, true, retain_original_branch_len).c_str());
            fclose(current_tree_file);

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        } else {
            if ((sort_before_placement_1 || sort_before_placement_2) && (missing_samples.size() > 1)) {
                timer.Start();
                fprintf(stderr, "Computing parsimony scores and number of parsimony-optimal placements for new samples and using them to sort the samples.\n");
                if (max_trees > 1) {
                    fprintf(stderr, "WARNING: --multiple-placements option is used but note that the samples will be sorted only once using the parsimony scores on input tree (without actual placements).\n");
                }

                // vectors to store the best parsimony scores and the number of
                // parsimony-optimal placements for each of the new samples to
                // be placed on the tree
                std::vector<int> best_parsimony_scores;
                std::vector<size_t> num_best_placements;

                for (size_t s=0; s<missing_samples.size(); s++) {

                    //Sort the missing sample mutations by position
                    std::sort(missing_samples[s].mutations.begin(), missing_samples[s].mutations.end());

                    auto bfs = T->breadth_first_expansion();
                    size_t total_nodes = bfs.size();

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

                    // Stores the parsimony score to place the sample at each
                    // node of the tree in DFS order.
                    std::vector<int> node_set_difference;
                    if (print_parsimony_scores) {
                        node_set_difference.resize(total_nodes);
                    }

                    size_t best_node_num_leaves = 0;
                    // The maximum number of mutations is bound by the number
                    // of mutations in the missing sample (place at root)
                    //int best_set_difference = 1e9;
                    // TODO: currently number of root mutations is also added to
                    // this value since it forces placement as child but this
                    // could be changed later
                    int best_set_difference = missing_samples[s].mutations.size() + T->root->mutations.size() + 1;

                    size_t best_j = 0;
                    size_t num_best = 1;
                    bool best_node_has_unique = false;
                    MAT::Node* best_node = T->root;

                    std::vector<bool> node_has_unique(total_nodes, false);
                    std::vector<size_t> best_j_vec;
                    best_j_vec.emplace_back(0);

                    // Parallel for loop to search for most parsimonious
                    // placements. Real action happens within mapper2_body
                    static tbb::affinity_partitioner ap;
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {
                        for (size_t k=r.begin(); k<r.end(); ++k) {
                            mapper2_input inp;
                            inp.T = T;
                            inp.node = bfs[k];
                            inp.missing_sample_mutations = &missing_samples[s].mutations;
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

                    best_parsimony_scores.emplace_back(best_set_difference);
                    num_best_placements.emplace_back(num_best);
                }

                // Sort samples order in indexes based on parsimony scores
                // and number of parsimony-optimal placements
                if (sort_before_placement_1) {
                    std::stable_sort(indexes.begin(), indexes.end(),
                    [&num_best_placements, &best_parsimony_scores](size_t i1, size_t i2) {
                        return ((best_parsimony_scores[i1] < best_parsimony_scores[i2]) || \
                                ((best_parsimony_scores[i1] == best_parsimony_scores[i2]) && (num_best_placements[i1] < num_best_placements[i2])));
                    });
                } else if (sort_before_placement_2) {
                    std::stable_sort(indexes.begin(), indexes.end(),
                    [&num_best_placements, &best_parsimony_scores](size_t i1, size_t i2) {
                        return ((num_best_placements[i1] < num_best_placements[i2]) || \
                                ((num_best_placements[i1] == num_best_placements[i2]) && (best_parsimony_scores[i1] < best_parsimony_scores[i2])));
                    });
                }

                // Reverse sorted order if specified
                if (reverse_sort) {
                    std::reverse(indexes.begin(), indexes.end());
                }

                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }

            fprintf(stderr, "Adding missing samples to the tree.\n");
        }

        std::string placement_stats_filename = outdir + "/placement_stats.tsv";
        FILE *placement_stats_file = fopen(placement_stats_filename.c_str(), "w");

        // Traverse in sorted sample order
        for (size_t idx=0; idx<indexes.size(); idx++) {

            num_trees = optimal_trees.size();

            for (size_t t_idx=0; t_idx < num_trees; t_idx++) {
                timer.Start();

                T = &optimal_trees[t_idx];

                if (num_trees > 1) {
                    fprintf(stderr, "==Tree %zu=== \n", t_idx+1);
                }

                size_t s = indexes[idx];
                auto sample = missing_samples[s].name;

                if (T->get_node(sample) != NULL) {
                    fprintf(stderr, "WARNING: Sample %s already in the tree! Ignoring.\n\n", sample.c_str());
                    continue;
                }

                if (print_parsimony_scores) {
                    auto parsimony_scores_filename = outdir + "/parsimony-scores.tsv";
                    if (s==0) {
                        fprintf(stderr, "\nNow computing branch parsimony scores for adding the missing samples at each of the %zu nodes in the existing tree without modifying the tree.\n", T->breadth_first_expansion().size());
                        fprintf(stderr, "The branch parsimony scores will be written to file %s\n\n", parsimony_scores_filename.c_str());

                        parsimony_scores_file = fopen(parsimony_scores_filename.c_str(), "w");
                        fprintf (parsimony_scores_file, "#Sample\tTree node\tParsimony score\tOptimal (y/n)\tParsimony-increasing mutations (for optimal nodes)\n");
                    }
                }

                auto bfs = T->breadth_first_expansion();
                size_t total_nodes = bfs.size();

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

                if (print_parsimony_scores) {
                    node_set_difference.resize(total_nodes);
                }

                size_t best_node_num_leaves = 0;
                // The maximum number of mutations is bound by the number
                // of mutations in the missing sample (place at root)
                //int best_set_difference = 1e9;
                // TODO: currently number of root mutations is also added to
                // this value since it forces placement as child but this
                // could be changed later
                int best_set_difference = missing_samples[s].mutations.size() + T->root->mutations.size() + 1;

                size_t best_j = 0;
                bool best_node_has_unique = false;

                std::vector<bool> node_has_unique(total_nodes, false);
                std::vector<size_t> best_j_vec;
                best_j_vec.emplace_back(0);

                size_t num_best = 1;
                MAT::Node* best_node = T->root;

                // Parallel for loop to search for most parsimonious
                // placements. Real action happens within mapper2_body
                static tbb::affinity_partitioner ap;
                tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k) {
                        mapper2_input inp;
                        inp.T = T;
                        inp.node = bfs[k];
                        inp.missing_sample_mutations = &missing_samples[s].mutations;
                        inp.excess_mutations = &node_excess_mutations[k];
                        inp.imputed_mutations = &node_imputed_mutations[k];
                        inp.best_node_num_leaves = &best_node_num_leaves;
                        inp.best_set_difference = &best_set_difference;
                        inp.best_node = &best_node;
                        inp.best_j =  &best_j;
                        inp.num_best = &num_best;
                        inp.j = k;
                        inp.has_unique = &best_node_has_unique;

                        if (print_parsimony_scores) {
                            inp.set_difference = &node_set_difference[k];
                        }
                        inp.best_j_vec = &best_j_vec;
                        inp.node_has_unique = &(node_has_unique);

                        mapper2_body(inp, print_parsimony_scores, print_parsimony_scores);
                    }
                }, ap);

                if (!print_parsimony_scores) {
                    best_set_difference += 1;

                    auto tmp_vec = std::vector<size_t>(best_j_vec.begin(), best_j_vec.end());

                    num_best = 0;
                    best_j_vec.clear();

                    // Parallel for loop to search for most parsimonious
                    // placements. Real action happens within mapper2_body
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, tmp_vec.size()),
                    [&](tbb::blocked_range<size_t> r) {
                        for (size_t l=r.begin(); l<r.end(); ++l) {
                            auto k = tmp_vec[l];
                            mapper2_input inp;
                            inp.T = T;
                            inp.node = bfs[k];
                            inp.missing_sample_mutations = &missing_samples[s].mutations;
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

                    fprintf(stderr, "Current tree size (#nodes): %zu\tSample name: %s\tParsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", total_nodes, sample.c_str(), \
                            best_set_difference, num_best);
                    fprintf(placement_stats_file, "%s\t%d\t%zu\t", sample.c_str(), best_set_difference, num_best);
                    // Prints a warning message if 2 or more
                    // parsimony-optimal placements found
                    if (num_best > 1) {
                        if (max_trees == 1) {
                            low_confidence_samples.emplace_back(sample);
                        }
                        if (num_best > max_uncertainty) {
                            fprintf(stderr, "WARNING: Number of parsimony-optimal placements exceeds maximum allowed value (%u). Ignoring sample %s.\n", max_uncertainty, sample.c_str());
                        }  else if (best_set_difference <= max_parsimony) {
                            fprintf(stderr, "WARNING: Multiple parsimony-optimal placements found. Placement done without high confidence.\n");
                        }
                    }

                    if (best_set_difference > max_parsimony) {
                        fprintf(stderr, "WARNING: Parsimony score of the most parsimonious placement exceeds the maximum allowed value (%u). Ignoring sample %s.\n", max_parsimony, sample.c_str());
                    }
                }  else {
                    fprintf(stderr, "Missing sample: %s\t Best parsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", sample.c_str(), \
                            best_set_difference, num_best);

                }

                // Debugging information to be printed if -DDEBUG compile-time
                // flag is set. This includes sample mutations, details of the
                // best node and the list of mutations at the best node
#if DEBUG == 1
                fprintf (stderr, "Sample mutations:\t");
                if (missing_samples[s].mutations.size() > 0) {
                    for (auto m: missing_samples[s].mutations) {
                        if (m.is_missing) {
                            continue;
                        }
                        fprintf(stderr, "|%s", (MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc)).c_str());
                        fprintf(stderr, "| ");
                    }
                }
                fprintf (stderr, "\n");

                assert(num_best > 0);

                //best_node_vec.emplace_back(best_node);
                if ((num_best > 0) && (num_best <= max_uncertainty) && (best_set_difference <= max_parsimony)) {
                    for (auto j: best_j_vec) {
                        auto node = bfs[j];

                        std::vector<std::string> muts;

                        fprintf(stderr, "Best node ");
                        if (node->is_leaf() || node_has_unique[j]) {
                            fprintf(stderr, "(sibling)");
                        } else {
                            fprintf(stderr, "(child)");
                        }

                        if (node == best_node) {
                            fprintf(stderr, "*: %s\t", node->identifier.c_str());
                        } else {
                            fprintf(stderr, ": %s\t", node->identifier.c_str());
                        }

                        std::string s = "|";
                        for (auto m: node->mutations) {
                            s += MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + '|';
                        }
                        if (node->mutations.size() > 0) {
                            muts.emplace_back(std::move(s));
                        }

                        for (auto anc: T->rsearch(node->identifier)) {
                            s = "|";
                            for (auto m: anc->mutations) {
                                s += MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + '|';
                            }
                            if (anc->mutations.size() > 0) {
                                muts.emplace_back(std::move(s));
                            }
                        }


                        std::reverse(muts.begin(), muts.end());

                        fprintf(stderr, "Mutations: ");
                        for (size_t m = 0; m < muts.size(); m++) {
                            fprintf(stderr, "%s", muts[m].c_str());
                            if (m+1 < muts.size()) {
                                fprintf(stderr, " > ");
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                    fprintf(stderr, "\n");
                }

#endif

                // If number of parsimony-optimal trees is more than 1 and if
                // the number of trees has not already exceeded the maximum
                // limit, create a copy of the current tree in curr_tree
                MAT::Tree curr_tree;
                if ((max_trees > 1) && (num_best > 1) && (num_trees < max_trees)) {
                    curr_tree = MAT::get_tree_copy(*T);
                }

                if (print_parsimony_scores) {
                    for (size_t k = 0; k < total_nodes; k++) {
                        char is_optimal = (node_set_difference[k] == best_set_difference) ? 'y' : 'n';
                        fprintf (parsimony_scores_file, "%s\t%s\t%d\t\t%c\t", sample.c_str(), bfs[k]->identifier.c_str(), node_set_difference[k], is_optimal);
                        if (node_set_difference[k] == best_set_difference) {
                            if (node_set_difference[k] == 0) {
                                fprintf(parsimony_scores_file, "*");
                            }
                            for (size_t idx = 0; idx < static_cast<size_t>(node_set_difference[k]); idx++) {
                                auto m = node_excess_mutations[k][idx];
                                assert (m.is_masked() || ((m.mut_nuc & (m.mut_nuc-1)) == 0));
                                fprintf(parsimony_scores_file, "%s", m.get_string().c_str());
                                if (idx+1 < static_cast<size_t>(node_set_difference[k])) {
                                    fprintf(parsimony_scores_file, ",");
                                }
                            }
                        } else {
                            fprintf(parsimony_scores_file, "N/A");
                        }
                        fprintf(parsimony_scores_file, "\n");
                    }
                }
                // Do placement only if number of parsimony-optimal placements
                // does not exceed the maximum allowed value and the parsimony
                // score for the most parsimonious placement does not exceed
                // the maximum allowed value
                else if ((num_best <= max_uncertainty) && (best_set_difference <= max_parsimony)) {
                    if (num_best > 1) {
                        if (max_trees > 1) {
                            // Sorting by bfs order ensures reproducible results
                            // during multiple placements
                            std::sort(best_j_vec.begin(), best_j_vec.end());
                        }

                        // Update num_best so that the number of trees does
                        // not exceed maximum limit
                        if ((optimal_trees.size() <= max_trees) && (num_best + optimal_trees.size() > max_trees)) {
                            if ((num_best + optimal_trees.size() > max_trees+1) && (max_trees > 1))
                                fprintf (stderr, "%zu parsimony-optimal placements found but total trees has already exceed the max possible value (%i)!\n", num_best, max_trees);
                            num_best = 1 + max_trees - optimal_trees.size();
                        }
                    }

                    // Assign clades if maximum number of trees is 1
                    if (max_trees == 1) {
                        missing_samples[s].clade_assignments.clear();
                        missing_samples[s].clade_assignments.resize(T->get_num_annotations());
                        missing_samples[s].best_clade_assignment.clear();
                        missing_samples[s].best_clade_assignment.resize(T->get_num_annotations());
                        for (size_t c=0; c < T->get_num_annotations(); c++) {
                            missing_samples[s].clade_assignments[c].resize(best_j_vec.size());
                            //TODO: can be parallelized
                            for (size_t k=0; k < best_j_vec.size(); k++) {
                                bool include_self = !bfs[best_j_vec[k]]->is_leaf() && !node_has_unique[best_j_vec[k]];
                                auto clade_assignment = T->get_clade_assignment(bfs[best_j_vec[k]], c, include_self);
                                missing_samples[s].clade_assignments[c][k] = clade_assignment;
                                if (bfs[best_j_vec[k]]==best_node) {
                                    missing_samples[s].best_clade_assignment[c] = clade_assignment;
                                }
                            }
                            std::sort(missing_samples[s].clade_assignments[c].begin(), missing_samples[s].clade_assignments[c].end());
                        }
                    }

                    // Iterate over the number of parsimony-optimal placements
                    // for which a new tree will be created
                    for (size_t k = 0; k < num_best; k++) {

                        // best_j is updated using best_j_vec if multiple
                        // placements are allowed and the number of new trees
                        // for the given sample is greater than 1. If not, the
                        // default tie-breaking strategy used in mapper2_body has
                        // already chosen a single best_j
                        if ((max_trees > 1) && (num_best > 1)) {
                            if ((k==0) && (num_best > 1)) {
                                fprintf (stderr, "Creating %zu additional tree(s) for %zu parsimony-optimal placements.\n", num_best-1, num_best);
                            }
                            // If at second placement or higher, a new tree needs to
                            // be added to optimal_trees and T needs to point to its
                            // last element. If not, T is already pointing to the
                            // last element of optimal_trees on which placement will
                            // be carried out
                            if (k > 0) {
                                auto tmp_T = MAT::get_tree_copy(curr_tree);
                                optimal_trees.emplace_back(std::move(tmp_T));
                                T = &optimal_trees[optimal_trees.size()-1];
                                bfs = T->breadth_first_expansion();
                            }

                            best_j = best_j_vec[k];
                            best_node_has_unique = node_has_unique[k];
                            best_node = bfs[best_j];
                        }

                        // Add sample to tree unless --no-add or it is already in the tree
                        if (!no_add && T->get_node(sample) == NULL) {
                            // Is placement as sibling
                            if (best_node->is_leaf() || best_node_has_unique) {
                                std::string nid = T->new_internal_node_id();
                                T->create_node(nid, best_node->parent->identifier);
                                T->create_node(sample, nid);
                                T->move_node(best_node->identifier, nid);
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
                                    T->get_node(nid)->add_mutation(m);
                                }
                                // Add mutations to best node using l1_mut
                                for (auto m: l1_mut) {
                                    T->get_node(best_node->identifier)->add_mutation(m);
                                }
                                // Add new sample mutations using l2_mut
                                for (auto m: l2_mut) {
                                    T->get_node(sample)->add_mutation(m);
                                }
                            }
                            // Else placement as child
                            else {
                                T->create_node(sample, best_node->identifier);
                                MAT::Node* node = T->get_node(sample);
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

                        if (node_imputed_mutations[best_j].size() > 0) {
                            fprintf (stderr, "Imputed mutations:\t");
                            size_t tot = node_imputed_mutations[best_j].size();
                            for (size_t curr = 0; curr < tot; curr++) {
                                MAT::Mutation& mut = node_imputed_mutations[best_j][curr];
                                if (curr < tot-1) {
                                    fprintf (stderr, "%i:%c;", mut.position, MAT::get_nuc(mut.mut_nuc));
                                    fprintf (placement_stats_file, "%i:%c;", mut.position, MAT::get_nuc(mut.mut_nuc));
                                } else {
                                    fprintf (stderr, "%i:%c", mut.position, MAT::get_nuc(mut.mut_nuc));
                                    fprintf (placement_stats_file, "%i:%c", mut.position, MAT::get_nuc(mut.mut_nuc));
                                }
                            }
                            fprintf(stderr, "\n");
                        }

                        if (max_trees == 1) {
                            break;
                        }
                    }
                }
                fputc('\n', placement_stats_file);

                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }
        }
        fclose(placement_stats_file);
    }

    num_trees = optimal_trees.size();

    // If user specified print_parsimony_scores, close corresponding file and
    // terminate normally
    if (print_parsimony_scores) {
        if (parsimony_scores_file) {
            fclose(parsimony_scores_file);
        }
        return 0;
    }

    // Collapse the output tree
    if (collapse_output_tree) {
        timer.Start();

        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();

            T = &optimal_trees[t_idx];

            if (num_trees > 1) {
                fprintf(stderr, "Collapsing output tree %lu.\n", t_idx+1);
            } else {
                fprintf(stderr, "Collapsing output tree.\n");
            }

            T->collapse_tree();

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    }

    // If user need uncondensed tree output, write uncondensed tree(s) to
    // file(s)
    if (print_uncondensed_tree) {
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();

            T = &optimal_trees[t_idx];

            auto uncondensed_final_tree_filename = outdir + "/uncondensed-final-tree.nh";
            if (num_trees > 1) {
                uncondensed_final_tree_filename = outdir + "/uncondensed-final-tree-" + std::to_string(t_idx+1) + ".nh";
                fprintf(stderr, "Writing uncondensed final tree %zu to file %s \n", (t_idx+1), uncondensed_final_tree_filename.c_str());
            } else {
                fprintf(stderr, "Writing uncondensed final tree to file %s \n", uncondensed_final_tree_filename.c_str());
            }

            auto parsimony_score = T->get_parsimony_score();
            fprintf(stderr, "The parsimony score for this tree is: %zu \n", parsimony_score);
            std::ofstream uncondensed_final_tree_file(uncondensed_final_tree_filename.c_str(), std::ofstream::out);
            std::stringstream newick_ss;
            write_newick_string(newick_ss, *T, T->root, true, true, retain_original_branch_len, true);
            uncondensed_final_tree_file << newick_ss.rdbuf();
            uncondensed_final_tree_file.close();

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    } else {
        // Write final tree(s) to file(s)
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();

            T = &optimal_trees[t_idx];

            auto final_tree_filename = outdir + "/final-tree.nh";
            if (num_trees > 1) {
                final_tree_filename = outdir + "/final-tree-" + std::to_string(t_idx+1) + ".nh";
                fprintf(stderr, "Writing final tree %zu to file %s \n", t_idx+1, final_tree_filename.c_str());
            } else {
                fprintf(stderr, "Writing final tree to file %s \n", final_tree_filename.c_str());
            }
            auto parsimony_score = T->get_parsimony_score();
            fprintf(stderr, "The parsimony score for this tree is: %zu \n", parsimony_score);

            std::ofstream final_tree_file(final_tree_filename.c_str(), std::ofstream::out);
            std::stringstream newick_ss;
            write_newick_string(newick_ss, *T, T->root, true, true, retain_original_branch_len);
            final_tree_file << newick_ss.rdbuf();
            final_tree_file.close();

            tree_parsimony_scores.emplace_back(parsimony_score);

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    }

    if (missing_samples.size() > 0) {
        // For each final tree write the path of mutations from tree root to the
        // sample for each newly placed sample
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();

            T = &optimal_trees[t_idx];
            bool use_tree_idx = false;
            if (num_trees > 1) {
                use_tree_idx = true;
            }
            std::vector<std::string> targets;
            for (auto s: missing_samples) {
                targets.emplace_back(s.name);
            }
            auto mutation_paths_filename = outdir + "/mutation-paths.txt";
            if (use_tree_idx) {
                mutation_paths_filename = outdir + "/mutation-paths-" + std::to_string(t_idx+1) + ".txt";
                fprintf(stderr, "Writing mutation paths for tree %zu to file %s \n", t_idx+1, mutation_paths_filename.c_str());
            } else {
                fprintf(stderr, "Writing mutation paths to file %s \n", mutation_paths_filename.c_str());
            }
            MAT::get_sample_mutation_paths(T, targets, mutation_paths_filename);
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
        // For each final tree write the annotations for each sample
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            T = &optimal_trees[t_idx];

            size_t num_annotations = T->get_num_annotations();

            if (num_annotations > 0) {
                timer.Start();

                auto annotations_filename = outdir + "/clades.txt";
                if (num_trees > 1) {
                    annotations_filename = outdir + "/clades" + std::to_string(t_idx+1) + ".txt";
                    fprintf(stderr, "Writing clade annotations for tree %zu to file %s \n", (t_idx+1), annotations_filename.c_str());
                } else {
                    fprintf(stderr, "Writing clade annotations to file %s \n", annotations_filename.c_str());
                }

                FILE* annotations_file = fopen(annotations_filename.c_str(), "w");

                for (size_t s=0; s<missing_samples.size(); s++) {
                    if (missing_samples[s].best_clade_assignment.size() == 0) {
                        // Sample was not placed (e.g. exceeded max EPPs) so no clades assigned
                        continue;
                    }
                    auto sample = missing_samples[s].name;

                    fprintf(annotations_file, "%s\t", sample.c_str());
                    for (size_t k=0; k< num_annotations; k++) {
                        fprintf(annotations_file, "%s", missing_samples[s].best_clade_assignment[k].c_str());
                        //TODO
                        if (max_trees == 1 && detailed_clades) {
                            fprintf(annotations_file, "*|");
                            std::string curr_clade = "";
                            int curr_count = 0;
                            for (auto clade: missing_samples[s].clade_assignments[k]) {
                                if (clade == curr_clade) {
                                    curr_count++;
                                } else {
                                    if (curr_count > 0) {
                                        fprintf(annotations_file, "%s(%i/%zu),", curr_clade.c_str(), curr_count,
                                                missing_samples[s].clade_assignments[k].size());
                                    }
                                    curr_clade = clade;
                                    curr_count = 1;
                                }
                            }
                            if (curr_count > 0) {
                                fprintf(annotations_file, "%s(%i/%zu)", curr_clade.c_str(), curr_count,
                                        missing_samples[s].clade_assignments[k].size());
                            }
                        }
                        if (k+1 < num_annotations) {
                            fprintf(annotations_file, "\t");
                        }
                    }
                    fprintf(annotations_file, "\n");
                }

                fclose(annotations_file);

                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }
        }
    }

    if ((print_subtrees_single > 1) && (missing_samples.size() > 0)) {
        fprintf(stderr, "Computing the single subtree for added samples with %zu random leaves. \n\n", print_subtrees_single);
        timer.Start();
        // For each final tree, write a subtree of user-specified size around
        // each newly placed sample in newick format
        std::vector<std::string> targets;
        for (auto m: missing_samples) {
            targets.emplace_back(m.name);
        }
        bool use_tree_idx = false;
        if (num_trees > 1) {
            use_tree_idx = true;
        }
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            optimal_trees[t_idx].uncondense_leaves();
            MAT::get_random_single_subtree(&optimal_trees[t_idx], targets, outdir, print_subtrees_single, t_idx, use_tree_idx, retain_original_branch_len);
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    if ((print_subtrees_size > 1) && (missing_samples.size() > 0)) {
        fprintf(stderr, "Computing subtrees for added samples. \n\n");

        // For each final tree, write a subtree of user-specified size around
        // each newly placed sample in newick format
        timer.Start();
        std::vector<std::string> targets;
        for (auto m: missing_samples) {
            targets.emplace_back(m.name);
        }
        bool use_tree_idx = false;
        if (num_trees > 1) {
            use_tree_idx = true;
        }
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            optimal_trees[t_idx].uncondense_leaves();
            MAT::get_random_sample_subtrees(&optimal_trees[t_idx], targets, outdir, print_subtrees_size, t_idx, use_tree_idx, retain_original_branch_len);
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    // Print warning message with a list of all samples placed with low
    // confidence (>=2 parsimony-optimal placements)
    if (low_confidence_samples.size() > 0) {
        fprintf(stderr, "WARNING: Following samples had multiple possibilities of parsimony-optimal placements:\n");
        for (auto lcs: low_confidence_samples) {
            fprintf(stderr, "%s\n", lcs.c_str());
        }
    }

    // Store mutation-annotated tree to a protobuf file if user has asked for it
    if (dout_filename != "") {

        timer.Start();

        fprintf(stderr, "Saving mutation-annotated tree object to file (after condensing identical sequences) %s\n", dout_filename.c_str());
        if (num_trees > 1) {
            fprintf(stderr, "WARNING: --multiple-placements option was used but only the first mutation-annotated tree object will be saved to file.\n");
        }

        Parsimony::data data;

        T = &optimal_trees[0];
        // Recondense tree with new samples
        if (T->condensed_nodes.size() > 0) {
            T->uncondense_leaves();
        }
        T->condense_leaves();
        MAT::save_mutation_annotated_tree(*T, dout_filename);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    /*
    // When multiple placements was used, print the tree names with lowest total
    // parsimony scores
    if (max_trees > 1) {
        size_t min_parsimony_score = *std::min_element(tree_parsimony_scores.begin(), tree_parsimony_scores.end());
        if (optimal_trees.size() > 1) {
            fprintf(stderr, "Trees below have the lowest parsimony score of %zu among all trees: \n", min_parsimony_score);
            for (size_t t_idx = 0; t_idx < optimal_trees.size(); t_idx++) {
                if (tree_parsimony_scores[t_idx] == min_parsimony_score) {
                    fprintf(stderr, "tree-%zu ", t_idx+1);
                }
            }
            fprintf(stderr, "\n\n");
        }
        else {
            fprintf(stderr, "Single best tree with a parsimony score of %zu was found during multiple placements.\n\n", min_parsimony_score);;
        }
    }
    */

    google::protobuf::ShutdownProtobufLibrary();

#if SAVE_PROFILE == 1
    Instrumentor::Get().EndSession();
#endif

    return 0;
}
