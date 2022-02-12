#include "mapper.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cassert>
#include <climits>
#include <cstdio>
#include <signal.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <vector>
Mutation_Set get_mutations(const MAT::Node *main_tree_node) {
    Mutation_Set out;
    while (main_tree_node) {
        for (const auto &mut : main_tree_node->mutations) {
                out.insert(mut);        
        }
        main_tree_node = main_tree_node->parent;
    }
    auto iter=out.begin();
    while (iter!=out.end()) {
        auto last_iter=iter;
        iter++;
        if (last_iter->get_mut_one_hot()==last_iter->get_ref_one_hot()) {
            out.erase(last_iter);
        }
    }
    return out;
}
static void set_covered_main_tree(const MAT::Node *start,
                                  std::vector<char> &checked, int distance_left,
                                  const MAT::Node *exclude_node) {
    checked[start->dfs_index] = true;
    int par_dist_left = distance_left - start->branch_length;
    if (start->parent && (start->parent != exclude_node) && distance_left > 0) {
        set_covered_main_tree(start->parent, checked, par_dist_left, start);
    }
    for (auto child : start->children) {
        int child_dist_left = distance_left - child->branch_length;
        if (distance_left > 0 && child != exclude_node) {
            set_covered_main_tree(child, checked, child_dist_left, start);
        }
    }
}
static void
check_mutation_helper(Mutation_Set &ref,
                      const MAT::Mutations_Collection &follows,
                      std::unordered_map<int, uint8_t> &parent_allele_check,
                      bool check_non_ambi) {
    for (const auto &mut : follows) {
        if (mut.get_position()==INT_MAX) {
            continue;
        }
        if (check_non_ambi &&
            (__builtin_popcount(mut.get_mut_one_hot()) != 1)) {
            raise(SIGTRAP);
        }
        auto ref_iter = ref.find(mut);
        if (ref_iter != ref.end()) {
            if (ref_iter->get_mut_one_hot() != mut.get_mut_one_hot()) {
                raise(SIGTRAP);
            }
            ref.erase(ref_iter);
        }   
        auto ins_result = parent_allele_check.emplace(mut.get_position(),
                                                      mut.get_par_one_hot());
        if (!ins_result.second) {
            if (ins_result.first->second != mut.get_mut_one_hot()) {
                raise(SIGTRAP);
            }
            ins_result.first->second = mut.get_par_one_hot();
        }
    }
}
static void
check_mutation_helper(Mutation_Set &ref,
                      const std::vector<To_Place_Sample_Mutation> &follows,
                      std::unordered_map<int, uint8_t> &parent_allele_check) {
    for (const auto &mut : follows) {
        if (mut.position==INT_MAX) {
            continue;
        }
        if (mut.mut_nuc==0xf) {
            for(int pos=mut.position;pos<=mut.get_end_range();pos++){
                auto ref_iter = ref.find(MAT::Mutation(mut.chrom_idx,pos,0,0xf));
                if (ref_iter != ref.end()) {
                    if (ref_iter->get_mut_one_hot() != 0xf) {
                        raise(SIGTRAP);
                    }
                    ref.erase(ref_iter);
                }        
            }
            continue;
        }
        auto ref_iter = ref.find(MAT::Mutation(mut.chrom_idx,mut.position,mut.par_nuc,mut.mut_nuc));
        if (ref_iter != ref.end()) {
            if (ref_iter->get_mut_one_hot() != mut.mut_nuc) {
                raise(SIGTRAP);
            }
            ref.erase(ref_iter);
        }   
        auto ins_result = parent_allele_check.emplace(mut.position,
                                                      mut.par_nuc);
        if(!ins_result.second){
            raise(SIGTRAP);
        }
    }
}

void check_ancestor(const MAT::Node* parent_node,Mutation_Set &ref,std::unordered_map<int, uint8_t>& parent_allele_check){
    while (parent_node&&(!ref.empty())) {
        check_mutation_helper(ref, parent_node->mutations, parent_allele_check,
                              true);
        parent_node = parent_node->parent;
    }
    if (!ref.empty()) {
        bool met=false;
        for (const auto & mut : ref) {
            fprintf(stderr, "pos: %d, mut_nuc %d, par_nuc %d \n",mut.get_position(),(int)mut.get_mut_one_hot(),(int)mut.get_par_one_hot());
            met=true;
        }
        if (met) {
            raise(SIGTRAP);        
        }
    }
    /*for (const auto &temp : parent_allele_check) {
        if (MAT::Mutation::refs[temp.first] != temp.second) {
            raise(SIGTRAP);
        }
    }*/
}
void check_mutations(Mutation_Set ref,
                     const MAT::Mutations_Collection &sample_mutations,
                     const MAT::Mutations_Collection &shared_mutations,
                     const MAT::Node *parent_node) {
    std::unordered_map<int, uint8_t> parent_allele_check;
    parent_allele_check.reserve(MAT::Mutation::refs.size());
    check_mutation_helper(ref, sample_mutations, parent_allele_check, false);
    check_mutation_helper(ref, shared_mutations, parent_allele_check, true);
    check_ancestor(parent_node, ref, parent_allele_check);

}
void check_mutations(Mutation_Set ref,
                     const std::vector<To_Place_Sample_Mutation> &sample_mutations,
                     const MAT::Mutations_Collection &shared_mutations,
                     const MAT::Node *parent_node) {
    std::unordered_map<int, uint8_t> parent_allele_check;
    parent_allele_check.reserve(MAT::Mutation::refs.size());
    check_mutation_helper(ref, sample_mutations, parent_allele_check);
    check_mutation_helper(ref, shared_mutations, parent_allele_check, true);
    check_ancestor(parent_node, ref, parent_allele_check);

}

void check_continuation(const MAT::Node* parent_node,Mutation_Set ref,const std::vector<To_Place_Sample_Mutation> &decendent_mutations){
    assert(decendent_mutations.back().position==INT_MAX);
    std::unordered_map<int, uint8_t> parent_allele_check;
    parent_allele_check.reserve(MAT::Mutation::refs.size());
    check_mutation_helper(ref, decendent_mutations, parent_allele_check);
    check_ancestor(parent_node, ref, parent_allele_check);
}
void check_mutations(Mutation_Set ref,
                     const Main_Tree_Target &target_to_check) {
    assert(target_to_check.sample_mutations.back().position==INT_MAX);
    check_mutations(ref, target_to_check.sample_mutations,
                    target_to_check.shared_mutations,
                    target_to_check.parent_node);
    auto split_node_mutations = get_mutations(target_to_check.target_node);
    check_mutations(split_node_mutations, target_to_check.splited_mutations,
                    target_to_check.shared_mutations,
                    target_to_check.parent_node);
}
void check_descendant_nuc(const MAT::Node* node){
    std::unordered_map<int, uint8_t> muts;
    muts.reserve(node->mutations.size());
    for(const auto& mut:node->mutations){
        muts.emplace(mut.get_position(),mut.get_mut_one_hot());
    }
    while (node) {
        for(const auto& mut:node->mutations){
            auto iter=muts.find(mut.get_position());
            if (iter!=muts.end()) {
                if ((iter->second&mut.get_descendant_mut())!=iter->second) {
                    raise(SIGTRAP);
                }
            }
        }
        node=node->parent;
    }
}