#include "mapper.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <utility>
#include <vector>
static void update_possible_descendant_alleles(
    const MAT::Mutations_Collection &mutations_to_set,
    MAT::Node *node) {
    std::unordered_map<int, uint8_t> alleles;
    alleles.reserve(mutations_to_set.size());
    for (auto &mut : mutations_to_set) {
        alleles.emplace(mut.get_position(), mut.get_mut_one_hot());
    }
    while (!alleles.empty() && node) {
        for (auto &mut : node->mutations) {
            auto iter = alleles.find(mut.get_position());
            if (iter != alleles.end()) {
                if ((mut.get_descendant_mut() & iter->second) ==
                    iter->second) {
                    alleles.erase(iter);
                } else {
                    mut.set_descendant_mut(mut.get_descendant_mut()|iter->second);
                }
            }
        }
        node = node->parent;
    }
}
static void gather_par_mutation_step(std::unordered_map<int, int>& to_find,MAT::Mutations_Collection& upstream, MAT::Mutations_Collection& output){
    for (const auto& mut : upstream) {
        auto iter=to_find.find(mut.get_position());
        if (iter!=to_find.end()) {
            output[iter->second].set_par_one_hot(mut.get_mut_one_hot());
            to_find.erase(iter);
        }
    }
}
static void gather_par_mut(std::unordered_map<int, int>& to_find,MAT::Node*& node, MAT::Mutations_Collection& output){
    while (node&&(!to_find.empty())) {
        gather_par_mutation_step(to_find,node->mutations,output);
        node=node->parent;
    }
    for(auto& temp: to_find){
        output[temp.second].set_par_one_hot(output[temp.second].get_ref_one_hot());
    }
}

static void discretize_mutations(std::vector<To_Place_Sample_Mutation> &in,
                                 MAT::Mutations_Collection &shared_mutations,
                                 MAT::Node *parent_node,
                                 MAT::Mutations_Collection &out) {
    out.reserve(in.size());
    std::unordered_map<int, int> par_nuc_idx;
    assert(in.back().position==INT_MAX);
    for (size_t idx=0;idx<(in.size()-1) ; idx++) {
        const auto & mut=in[idx];
        if (mut.mut_nuc == 0xf) {
            for (int pos = mut.position; pos <= mut.get_end_range(); pos++) {
                par_nuc_idx.emplace(pos, out.size());
                out.push_back(MAT::Mutation(mut.chrom_idx, pos, 0, 0xf));
                out.back().set_descendant_mut(0xf);
            }
        } else {
            out.push_back(MAT::Mutation(mut.chrom_idx, mut.position,
                                        mut.par_nuc, mut.mut_nuc));
            out.back().set_descendant_mut(mut.mut_nuc);
        }
    }
    gather_par_mutation_step(par_nuc_idx, shared_mutations,out);
    gather_par_mut(par_nuc_idx, parent_node, out);
}
static const MAT::Node *update_main_tree(Main_Tree_Target &target,
                                         std::string &&sample_string) {
    // Split branch?
    MAT::Node *sample_node = new MAT::Node;
    sample_node->identifier = std::move(sample_string);
    discretize_mutations(target.sample_mutations, target.shared_mutations, target.parent_node, sample_node->mutations);
    int sample_node_mut_count = 0;
    for (const auto &mut : sample_node->mutations) {
        if (!(mut.get_par_one_hot() & mut.get_mut_one_hot())) {
            sample_node_mut_count++;
        }
        assert(mut.get_position());
    }
    sample_node->branch_length = sample_node_mut_count;
    if (target.splited_mutations.empty() && (!target.target_node->is_leaf())) {
        sample_node->parent = target.target_node;
        target.target_node->children.push_back(sample_node);
    } else if (target.shared_mutations.empty() &&
               (!target.target_node->is_leaf())) {
        sample_node->parent = target.parent_node;
        target.parent_node->children.push_back(sample_node);
    } else {
        target.target_node->mutations = std::move(target.splited_mutations);
        int target_node_mut_count = 0;
        for (const auto &mut : target.target_node->mutations) {
            if (!(mut.get_mut_one_hot() & mut.get_par_one_hot())) {
                target_node_mut_count++;
            }
            assert(mut.get_position());
        }
        target.target_node->branch_length = target_node_mut_count;
        MAT::Node *split_node = new MAT::Node;
        target.target_node->parent = split_node;
        sample_node->parent = split_node;
        split_node->identifier = "";
        split_node->parent = target.parent_node;
        split_node->mutations = std::move(target.shared_mutations);
        split_node->children.push_back(target.target_node);
        split_node->children.push_back(sample_node);
        split_node->branch_length = split_node->mutations.size();
        for (const auto& mut : split_node->mutations) {
            assert(mut.get_position());
        }
        auto iter =
            std::find(target.parent_node->children.begin(),
                      target.parent_node->children.end(), target.target_node);
        if (iter==target.parent_node->children.end()||*iter!=target.target_node) {
            std::raise(SIGTRAP);
        }
        *iter = split_node;
    }
    update_possible_descendant_alleles(sample_node->mutations, sample_node->parent);
    #ifndef NDEBUG
    check_descendant_nuc(sample_node);
    check_descendant_nuc(target.target_node);
    check_descendant_nuc(sample_node->parent);
    #endif
    return sample_node;
}
void place_sample(Sample_Muts &&sample_to_place, MAT::Tree &main_tree
#ifndef NDEBUG
                  ,
                  Original_State_t &ori_state
#endif
) {
    std::vector<MAT::Mutation> &&sample_mutations =
        std::move(sample_to_place.muts);
    std::string &&sample_string = std::move(sample_to_place.sample_name);
    /*if (sample_string=="s1433144s") {
        raise(SIGTRAP);
    }
    if (sample_string=="s2886812s") {
        raise(SIGTRAP);
    }
    if (sample_string=="s2749940s") {
        raise(SIGTRAP);
    }*/
    std::vector<To_Place_Sample_Mutation> condensed_muts;
    convert_mut_type(sample_mutations,condensed_muts);
#ifndef NDEBUG
    Mutation_Set new_set;
    std::vector<To_Place_Sample_Mutation> condensed_muts_copy(condensed_muts);
    new_set.reserve(sample_mutations.size());
    for (const auto &mut : sample_mutations) {
        new_set.insert(mut);
    }
    ori_state.emplace(sample_string, new_set);
#endif
    auto main_tree_start = std::chrono::steady_clock::now();
    auto main_tree_out =
        place_main_tree(condensed_muts, main_tree
#ifndef NDEBUG
                        ,
                        new_set
#endif
        );
    auto main_tree_duration =
        std::chrono::steady_clock::now() - main_tree_start;
    fprintf(stderr, "main tree took %ld msec\n",
            std::chrono::duration_cast<std::chrono::milliseconds>(
                main_tree_duration)
                .count());
    auto main_tree_node =
        update_main_tree(std::get<0>(main_tree_out), std::move(sample_string));
#ifndef NDEBUG
    /*std::vector<Sampled_Tree_Node *> output;
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(main_tree, output, sampling_radius);
    fprintf(stderr, "%zu samples \n", ori_state.size());*/
    //check_samples(main_tree.root, ori_state, &main_tree);
#endif
}