#include "mapper.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <tbb/parallel_for.h>
#include <utility>
#include <vector>
static void update_possible_descendant_alleles(
    std::vector<Sampled_Tree_Mutation> &mutations_to_set,
    Sampled_Tree_Node *node) {
    std::unordered_map<int, uint8_t> alleles;
    alleles.reserve(mutations_to_set.size());
    for (auto &mut : mutations_to_set) {
        alleles.emplace(mut.position, mut.mut_nuc);
    }
    while (!alleles.empty() && node) {
        for (auto &mut : node->mutations) {
            auto iter = alleles.find(mut.position);
            if (iter != alleles.end()) {
                if ((mut.descendent_possible_nuc & iter->second) ==
                    iter->second) {
                    alleles.erase(iter);
                } else {
                    mut.descendent_possible_nuc |= iter->second;
                }
            }
        }
        node = node->parent;
    }
}

static void update_sampled_tree(Sampled_Place_Target &target,
                                const MAT::Node *main_tree_node) {
    auto new_node = new Sampled_Tree_Node();
    new_node->corresponding_main_tree_node = main_tree_node;
    Sampled_Tree_Node *parent_node =
        const_cast<Sampled_Tree_Node *>(target.target_node);
    new_node->mutations = std::move(target.muts);
    new_node->parent = parent_node;
    parent_node->children.push_back(new_node);
    update_possible_descendant_alleles(new_node->mutations, new_node);
}

static const MAT::Node *update_main_tree(Main_Tree_Target &target,
                                         std::string &&sample_string) {
    // Split branch?
    MAT::Node *sample_node = new MAT::Node;
    sample_node->identifier = std::move(sample_string);
    sample_node->mutations = std::move(target.sample_mutations);
    int sample_node_mut_count=0;
    for (const auto & mut : sample_node->mutations) {
        if (!(mut.get_par_one_hot()&mut.get_mut_one_hot())) {
            sample_node_mut_count++;
        }
    }
    sample_node->branch_length=sample_node_mut_count;
    if (target.splited_mutations.empty()&&(!target.target_node->is_leaf())) {
        sample_node->parent = target.target_node;
        target.target_node->children.push_back(sample_node);
    } else if (target.shared_mutations.empty()&&(!target.target_node->is_leaf())) {
        sample_node->parent = target.parent_node;
        target.parent_node->children.push_back(sample_node);
    } else {
        MAT::Node *new_target_node = new MAT::Node;
        new_target_node->identifier=target.target_node->identifier;
        new_target_node->children = target.target_node->children;
        new_target_node->mutations = std::move(target.splited_mutations);
        int target_node_mut_count=0;
        for (const auto& mut : new_target_node->mutations) {
            if (!(mut.get_mut_one_hot()&mut.get_par_one_hot())) {
                target_node_mut_count++;
            }   
        }
        new_target_node->branch_length=target_node_mut_count;
        MAT::Node *split_node = new MAT::Node;
        new_target_node->parent=split_node;
        sample_node->parent=split_node;
        split_node->identifier = "";
        split_node->parent = target.parent_node;
        split_node->mutations = std::move(target.shared_mutations);
        split_node->children.push_back(new_target_node);
        split_node->children.push_back(sample_node);
        split_node->branch_length=split_node->mutations.size();
        auto iter=std::find(target.parent_node->children.begin(),target.parent_node->children.end(),target.target_node);
        *iter=split_node;
    }
    return sample_node;
}
void place_sample(Sample_Muts &&sample_to_place,
                  Sampled_Tree_Node *sampled_tree_root, MAT::Tree &main_tree,
                  int sampling_radius
#ifndef NDEBUG
                  ,
                  Original_State_t &ori_state
#endif
) {
    std::vector<Sampled_Tree_Mutation> &&sample_mutations =
        std::move(sample_to_place.muts);
    std::string &&sample_string = std::move(sample_to_place.sample_name);
#ifndef NDEBUG
    Mutation_Set new_set;
    new_set.reserve(sample_mutations.size());
    for (const auto & mut : sample_mutations) {
        new_set.insert(MAT::Mutation(mut.chrom_idx,mut.position,MAT::Mutation::refs[mut.position],mut.mut_nuc));
    }
    ori_state.emplace(sample_string, new_set);
#endif
    auto sampled_tree_start=std::chrono::steady_clock::now();
    int sampled_mutations = 0;
    auto sampled_output = place_on_sampled_tree(
        sampled_tree_root, std::move(sample_mutations), sampled_mutations
        #ifndef NDEBUG
                        ,
                        new_set
#endif
        );
    auto sampled_tree_duration=std::chrono::steady_clock::now()-sampled_tree_start;
    fprintf(stderr, "sampled tree took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(sampled_tree_duration).count());
    auto main_tree_start=std::chrono::steady_clock::now();
    auto main_tree_out =
        place_main_tree(sampled_output, main_tree, sampling_radius
#ifndef NDEBUG
                        ,
                        new_set
#endif
        );
    auto main_tree_duration=std::chrono::steady_clock::now()-main_tree_start;
    fprintf(stderr, "main tree took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(main_tree_duration).count());
    auto whole_tree_start=std::chrono::steady_clock::now();
#ifndef NDEBUG
    optimality_check(new_set, std::get<2>(main_tree_out),
                     main_tree.root, sampling_radius, sampled_tree_root,
                     sampled_output);
#endif
    auto whole_tree_duration=std::chrono::steady_clock::now()-whole_tree_start;
    fprintf(stderr, "whole tree took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(whole_tree_duration).count());
    auto dist=std::get<0>(main_tree_out).distance_left;
    auto main_tree_node =
        update_main_tree(std::get<0>(main_tree_out), std::move(sample_string));
    if (dist >= sampling_radius) {
        update_sampled_tree(sampled_output[std::get<1>(main_tree_out)],
                            main_tree_node);
    }
#ifndef NDEBUG
    std::vector<Sampled_Tree_Node *> output;
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(main_tree, output, sampling_radius);
    fprintf(stderr, "%zu samples \n",ori_state.size());
    check_samples(main_tree.root, ori_state, &main_tree);
#endif
}