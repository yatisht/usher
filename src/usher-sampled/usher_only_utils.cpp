#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "src/usher-sampled/usher.hpp"
#include <tbb/parallel_for.h>
extern int process_count;
void prep_tree(MAT::Tree &tree) {
    if (!tree.root->mutations.empty()) {
        auto node = tree.create_node();
        auto ori_root = tree.root;
        node->children.push_back(ori_root);
        ori_root->parent = node;
        tree.root = node;
    }
    assign_descendant_muts(tree);
    assign_levels(tree.root);
    set_descendant_count(tree.root);
}
bool sort_samples(const Leader_Thread_Options& options,std::vector<Sample_Muts>& samples_to_place, MAT::Tree& tree,size_t sample_start_idx) {
    bool reordered=false;
    if (options.sort_by_ambiguous_bases) {
        fprintf(stderr, "Sorting missing samples based on the number of ambiguous bases \n");
        tbb::parallel_for(tbb::blocked_range<size_t>(0,samples_to_place.size()),
        [&samples_to_place](tbb::blocked_range<size_t> range) {
            for (size_t idx=range.begin(); idx<range.end(); idx++) {
                int ambiguous_count=0;
                for(const auto& mut:samples_to_place[idx].muts) {
                    if (mut.mut_nuc==0xf) {
                        ambiguous_count+=mut.range;
                    } else if (__builtin_popcount(mut.mut_nuc)!=1) {
                        ambiguous_count++;
                    }
                }
                samples_to_place[idx].sorting_key1=ambiguous_count;
            }
        });
        if (options.reverse_sort) {
            fprintf(stderr, "Reverse sort \n");
            std::sort(samples_to_place.begin(), samples_to_place.end(),
            [](const auto& samp1,const auto& samp2) {
                return samp1.sorting_key1>samp2.sorting_key1;
            });
        } else {
            std::sort(samples_to_place.begin(), samples_to_place.end(),
            [](const auto& samp1,const auto& samp2) {
                return samp1.sorting_key1<samp2.sorting_key1;
            });

        }        // Reverse sorted order if specified
        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        reordered=true;
    }
    if (options.collapse_tree) {
        //timer.Start();

        fprintf(stderr, "Collapsing input tree.\n");

        auto condensed_tree_filename = options.out_options.outdir + "/condensed-tree.nh";
        fprintf(stderr, "Writing condensed input tree to file %s\n", condensed_tree_filename.c_str());

        FILE* condensed_tree_file = fopen(condensed_tree_filename.c_str(), "w");
        fprintf(condensed_tree_file, "%s\n",tree.get_newick_string(true,true,options.out_options.retain_original_branch_len).c_str());
        fclose(condensed_tree_file);

        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    std::vector<std::string> low_confidence_samples;
    std::vector<Clade_info> samples_clade;
    if (options.print_parsimony_scores) {
        // timer.Start();
        auto current_tree_filename = options.out_options.outdir + "/current-tree.nh";

        fprintf(
            stderr,
            "Writing current tree with internal nodes labelled to file %s \n",
            current_tree_filename.c_str());
        FILE *current_tree_file = fopen(current_tree_filename.c_str(), "w");
        fprintf(current_tree_file, "%s\n",
                tree.get_newick_string(true, true, options.out_options.retain_original_branch_len)
                .c_str());
        fclose(current_tree_file);

        // fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    } else {
        if ((options.sort_before_placement_1 || options.sort_before_placement_2) &&
                (samples_to_place.size() > 1)) {
            // timer.Start();
            fprintf(stderr, "Computing parsimony scores and number of "
                    "parsimony-optimal placements for new samples and "
                    "using them to sort the samples.\n");
            assign_descendant_muts(tree);
            assign_levels(tree.root);
            set_descendant_count(tree.root);
            if (process_count>1) {
                fprintf(stderr, "Main sending tree\n");
                tree.MPI_send_tree();
            }
            std::atomic_size_t curr_idx(0);
            place_sample_leader(samples_to_place, tree, 100, curr_idx, INT_MAX,
                                true, nullptr, options.max_parsimony,
                                options.max_uncertainty, low_confidence_samples,
                                samples_clade,sample_start_idx,nullptr);
            fputc('\n', stderr);
            //check_repeats(samples_to_place, sample_start_idx);
            // Sort samples order in indexes based on parsimony scores
            // and number of parsimony-optimal placements
            if (options.sort_before_placement_1) {
                std::sort(samples_to_place.begin(), samples_to_place.end(),
                [&options](const auto &samp1, const auto &samp2) {
                    if (samp1.sorting_key1 < samp2.sorting_key1) {
                        return options.reverse_sort;
                    }
                    if (samp1.sorting_key1 == samp2.sorting_key1 &&
                            samp1.sorting_key2 < samp2.sorting_key2) {
                        return options.reverse_sort;
                    }
                    return !options.reverse_sort;
                });
            } else if (options.sort_before_placement_2) {
                std::sort(samples_to_place.begin(), samples_to_place.end(),
                [&options](const auto &samp1, const auto &samp2) {
                    if (samp1.sorting_key2 < samp2.sorting_key2) {
                        return options.reverse_sort;
                    }
                    if (samp1.sorting_key2 == samp2.sorting_key2 &&
                            samp1.sorting_key1 < samp2.sorting_key1) {
                        return options.reverse_sort;
                    }
                    return !options.reverse_sort;
                });
            }
            // fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            reordered=true;
        }
        //check_repeats(samples_to_place, sample_start_idx);

        fprintf(stderr, "Adding missing samples to the tree.\n");
    }
    return reordered;
}
void clean_up_leaf(std::vector<MAT::Node*>& dfs) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dfs.size()),[&dfs](tbb::blocked_range<size_t> range) {
        for (size_t node_idx=range.begin(); node_idx<range.end(); node_idx++) {
            auto node=dfs[node_idx];
            if (node->is_leaf()) {
                auto& muts=node->mutations.mutations;
                muts.erase(std::remove_if(muts.begin(), muts.end(), [](const auto& mut) {
                    return mut.get_par_one_hot()&mut.get_mut_one_hot();
                }),muts.end());
                for (auto &mut : muts) {
                    /*if (mut.get_par_one_hot()&mut.get_mut_one_hot()) {
                        raise(SIGTRAP);
                    }*/
                    mut.set_mut_one_hot(1<<__builtin_ctz(mut.get_mut_one_hot()));
                }
            }
        }
    });
}
