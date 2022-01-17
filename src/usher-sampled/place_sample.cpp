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

static void get_par_mutations(
    const std::vector<To_Place_Sample_Mutation> &downstream_mutations,
    const MAT::Mutations_Collection &backward,
    std::vector<Sampled_Tree_Mutation> &output,
    const Sampled_Tree_Node *attach_node) {
        output.reserve(backward.size()+downstream_mutations.size());
        auto iter=downstream_mutations.begin();
        std::unordered_map<int, int> position_map;
        position_map.reserve(backward.size()+downstream_mutations.size());
        for (const auto& mut : backward) {
            while (iter->get_end_range()<mut.get_position()) {
                if (iter->mut_nuc!=0xf&&!(iter->mut_nuc&iter->par_nuc)) {
                    assert(__builtin_popcount(iter->mut_nuc)==1);
                    output.push_back(Sampled_Tree_Mutation{iter->position,iter->chrom_idx,iter->mut_nuc,iter->mut_nuc,iter->par_nuc});
                }
                iter++;
            }
            assert(mut.get_position()<=iter->get_end_range());
            //coincide
            if (iter->get_end_range()==mut.get_position()) {
                //Different
                position_map.emplace(mut.get_position(),output.size());
                output.push_back(
                        Sampled_Tree_Mutation{mut.get_position(),mut.get_chromIdx()
                        ,mut.get_par_one_hot(),mut.get_par_one_hot(),0});
                iter++;    
            }else if (mut.get_position()>=iter->position) {
                position_map.emplace(mut.get_position(),output.size());
                output.push_back(
                        Sampled_Tree_Mutation{mut.get_position(),mut.get_chromIdx()
                        ,mut.get_par_one_hot(),mut.get_par_one_hot(),0});            
            }else{
                //back mutation
                if (mut.get_par_one_hot()&mut.get_mut_one_hot()) {
                    fprintf(stderr, "%d not back mutation, iter position %d",mut.get_position(),iter->position);
                    raise(SIGTRAP);

                }
                    position_map.emplace(mut.get_position(),output.size());
                    output.push_back(
                        Sampled_Tree_Mutation{mut.get_position(),mut.get_chromIdx()
                        ,mut.get_par_one_hot(),mut.get_par_one_hot(),0});
            }
        }
        while (iter<(downstream_mutations.end()-1)) {
            if (!(iter->mut_nuc&iter->par_nuc)) {
                assert(__builtin_popcount(iter->mut_nuc)==1);
                output.push_back(Sampled_Tree_Mutation{iter->position,iter->chrom_idx,iter->mut_nuc,iter->mut_nuc,iter->par_nuc});
            }
            iter++;
        }
        while (attach_node&&(!position_map.empty())) {
            for (const auto& mut :attach_node->mutations ) {
                auto iter=position_map.find(mut.position);
                if (iter!=position_map.end()) {
                    output[iter->second].par_nuc=mut.mut_nuc;
                    position_map.erase(iter);
                }
            }
            attach_node=attach_node->parent;
        }
        for (auto& temp : position_map) {
            output[temp.second].par_nuc=MAT::Mutation::refs[temp.first];
        }
        output.erase(std::remove_if(output.begin(), output.end(), [](const auto& mut){
            return mut.par_nuc&mut.mut_nuc;
        }),output.end());
    }
static void update_sampled_tree(Sampled_Place_Target &target,
                                const MAT::Node *main_tree_node) {
    if (main_tree_node->parent ==
        target.target_node->corresponding_main_tree_node) {
        return;
    }
    auto new_node = new Sampled_Tree_Node();
    new_node->corresponding_main_tree_node = main_tree_node->parent;
    std::vector<Sampled_Tree_Mutation> new_node_mutations;
    Sampled_Tree_Node *parent_node =
        const_cast<Sampled_Tree_Node *>(target.target_node);
    get_par_mutations(target.muts, main_tree_node->mutations, new_node->mutations, parent_node);
    new_node->parent = parent_node;
    parent_node->children.push_back(new_node);
    update_possible_descendant_alleles(new_node->mutations, new_node);
    #ifndef NDEBUG
    check_sampled_main_correspondence(new_node);
    #endif
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
            }
        } else {
            out.push_back(MAT::Mutation(mut.chrom_idx, mut.position,
                                        mut.par_nuc, mut.mut_nuc));
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
        new_set.insert(MAT::Mutation(mut.chrom_idx, mut.position,
                                     MAT::Mutation::refs[mut.position],
                                     mut.mut_nuc));
    }
    ori_state.emplace(sample_string, new_set);
#endif
    auto sampled_tree_start = std::chrono::steady_clock::now();
    int sampled_mutations = 0;
    auto sampled_output = place_on_sampled_tree(
        sampled_tree_root, std::move(condensed_muts), sampled_mutations
#ifndef NDEBUG
        ,
        new_set
#endif
    );
    auto sampled_tree_duration =
        std::chrono::steady_clock::now() - sampled_tree_start;
    fprintf(stderr, "sampled tree took %ld msec\n",
            std::chrono::duration_cast<std::chrono::milliseconds>(
                sampled_tree_duration)
                .count());
    auto main_tree_start = std::chrono::steady_clock::now();
    auto main_tree_out =
        place_main_tree(sampled_output, main_tree, sampling_radius
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
#ifndef NDEBUG
    auto whole_tree_start = std::chrono::steady_clock::now();
    optimality_check(new_set, std::get<2>(main_tree_out), main_tree,
                     sampling_radius, sampled_tree_root, sampled_output,condensed_muts_copy);
    auto whole_tree_duration =
        std::chrono::steady_clock::now() - whole_tree_start;
    fprintf(stderr, "whole tree took %ld msec\n",
            std::chrono::duration_cast<std::chrono::milliseconds>(
                whole_tree_duration)
                .count());
#endif
    auto dist_left = std::get<0>(main_tree_out).distance_left;
    auto main_tree_node =
        update_main_tree(std::get<0>(main_tree_out), std::move(sample_string));
    for (const auto & temp : sampled_output) {
        assert(temp.muts.back().position==INT_MAX);
    }
    if (dist_left <= 0) {
        update_sampled_tree(sampled_output[std::get<1>(main_tree_out)],
                            main_tree_node);
    }
#ifndef NDEBUG
    /*std::vector<Sampled_Tree_Node *> output;
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(main_tree, output, sampling_radius);
    fprintf(stderr, "%zu samples \n", ori_state.size());*/
    //check_samples(main_tree.root, ori_state, &main_tree);
#endif
}