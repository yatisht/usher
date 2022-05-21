#include "place_sample.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>
static void place_sample(MAT::Tree &main_tree, const Sample_Muts &samp,
               Main_Tree_Target &target) {
    MAT::Mutations_Collection sample_mutations;
    discretize_mutations(target.sample_mutations, target.shared_mutations,
                         target.parent_node, sample_mutations);
    auto out = update_main_tree(sample_mutations, target.splited_mutations,
                                target.shared_mutations, target.target_node,
                                samp.sample_idx, main_tree, 0,true);

    if (out.deleted_nodes) {
        delete out.deleted_nodes;
    }
}
//void check_leaves(const MAT::Tree& T);
void place_sample_multiple_tree(
    std::vector<Sample_Muts> &sample_to_place,
    std::vector<MAT::Tree>& trees,
    FILE *placement_stats_file, int max_trees) {
    fprintf(stderr, "Max tree size %d\n",max_trees);
    for (const auto &samp : sample_to_place) {
        std::vector<std::tuple<std::vector<Main_Tree_Target>, int>>
            placement_result(trees.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, trees.size()),
                          [&samp, &trees, &placement_result](
                              tbb::blocked_range<size_t> range) {
                              for (size_t idx = range.begin();
                                   idx < range.end(); idx++) {
                                  placement_result[idx] =
                                      place_main_tree(samp.muts, trees[idx]);
                              }
                          });
        int min_par_score = INT_MAX;
        const auto& samp_name= trees[0].get_node_name(samp.sample_idx);
        for (size_t idx = 0; idx < placement_result.size(); idx++) {
            auto best_score=std::get<1>(placement_result[idx]);
            auto num_best=std::get<0>(placement_result[idx]).size();
            min_par_score = std::min(min_par_score,best_score);
            fprintf(stderr,
                    "Tree: %zu\tSample name: %s\tParsimony score: %d\tNumber "
                    "of parsimony-optimal placements: %zu\n",
                    idx, samp_name.c_str(), best_score, num_best);
            fprintf(placement_stats_file, "%s\t%d\t%zu\t", samp_name.c_str(), best_score, num_best);
        }
        for (size_t idx = 0; idx < trees.size(); idx++) {
            if (std::get<1>(placement_result[idx]) > min_par_score) {
                trees.erase(trees.begin() + idx);
                placement_result.erase(placement_result.begin() + idx);
            }
        }
        for (size_t idx = 0; idx < placement_result.size(); idx++) {
            auto to_place = std::get<0>(placement_result[idx]);
            if (to_place.size()>1) {
            if ((size_t)max_trees<trees.size()+to_place.size()) {
                fprintf(
                    stderr,
                    "%zu parsimony-optimal placements found but total trees "
                    "has already exceed the max possible value (%i)!\n",
                    to_place.size(), max_trees);
            }else {
                fprintf(stderr,
                        "Creating %zu additional tree(s) for %zu "
                        "parsimony-optimal placements.\n",
                        to_place.size() - 1, to_place.size());
            }
            }
            for (size_t place_idx = 1; place_idx < to_place.size();
                 place_idx++) {
                if (trees.size() >= (size_t)max_trees) {
                    break;
                }
                trees.push_back(trees[idx].copy_tree());
                auto temp=to_place[place_idx];
                temp.target_node=trees.back().get_node(temp.target_node->node_id);
                temp.parent_node=trees.back().get_node(temp.parent_node->node_id);
                place_sample(trees.back(), samp, temp);
            }
            place_sample(trees[idx], samp, to_place[0]);
        }
    }
}