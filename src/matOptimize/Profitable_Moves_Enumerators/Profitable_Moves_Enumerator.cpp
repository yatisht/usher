#include "process_each_node.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>
template<typename output_type>
static void
merge_mutation_LCA_to_dst(MAT::Node *child,
                          range<Mutation_Count_Change> mutations,
                          output_type &child_mutations) {
    //probably a gross overestimate, as we don't merge invalid mutations
    child_mutations.reserve(mutations.size()+child->mutations.size());
    auto child_mutation_iter = child->mutations.begin();
    auto child_mutation_end = child->mutations.end();
    for ( ; mutations; mutations++) {
        const Mutation_Count_Change &m=*mutations;
        //skip invalid mutations, src don't have these allele distribution among its children
        while (child_mutation_iter != child_mutation_end &&
                !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }
        while (child_mutation_iter != child_mutation_end &&
                *child_mutation_iter < m) {
            //newly encountered mutation to merge, add as back mutation
            child_mutations.emplace_back(*child_mutation_iter,0,child_mutation_iter->get_par_one_hot());
            //the par_nuc if src is moved under child
            child_mutations.back().set_par_nuc(
                child_mutation_iter->get_mut_one_hot());
            child_mutation_iter++;
        }
        if (child_mutation_iter != child_mutation_end &&
                child_mutation_iter->get_position() == m.get_position()) {
            //match, the state of src (incremented allele in m) dosen't change,
            //but need to update  par_nuc as moving under another node
            nuc_one_hot new_par_allele = child_mutation_iter->get_mut_one_hot();
            if (m.get_incremented()!=new_par_allele) {
                child_mutations.emplace_back(m);
                child_mutations.back().set_par_nuc(new_par_allele);
            }
            child_mutation_iter++;
        } else {
            //copy over other undisturbed mutations
            child_mutations.emplace_back(m);
        }
    }
    while (child_mutation_iter != child_mutation_end) {
        child_mutations.emplace_back(*child_mutation_iter,0,0);
        child_mutations.back().set_change(
            0, child_mutation_iter->get_par_one_hot());
        child_mutations.back().set_par_nuc(
            child_mutation_iter->get_mut_one_hot());
        child_mutation_iter++;
    }
}
static void
merge_mutation_src_to_LCA(const MAT::Node *root,
                          Mutation_Count_Change_Collection &mutations) {
    Mutation_Count_Change_Collection merged_mutations;
    merged_mutations.reserve(mutations.size() + root->mutations.size());
    auto iter = mutations.begin();
    auto end = mutations.end();
    for (const auto &m : root->mutations) {
        //similar idea as merge_mutation_LCA_to_dst, ignore invalid mutations
        if (!m.is_valid()) {
            continue;
        }
        while (iter != end && iter->get_position() < m.get_position()) {
            //copy over input mutations
            merged_mutations.push_back(*iter);
            iter++;
        }
        if (iter != end && iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_nuc = m.get_par_one_hot();
            //match, only append if it is actually a state change
            if (new_par_nuc != iter->get_incremented()) {
                merged_mutations.push_back(*iter);
                merged_mutations.back().set_par_nuc(m.get_par_one_hot());
            }
            iter++;
        } else {
            merged_mutations.emplace_back(m,0, m.get_mut_one_hot());
        }
    }
    while (iter!=end) {
        merged_mutations.push_back(*iter);
        iter++;
    }
    //Swap won't help, as the mutation vector will just grow as it moves up, so couldn't recycle
    mutations = std::move(merged_mutations);
}
//Need to preserve all auxilary information in mutation vector of src, to know the allele distribution of its children after the move
static void init_mutation_change(MAT::Node* src, Mutation_Count_Change_Collection& mutations) {
    mutations.reserve(src->mutations.size());
    for (const auto &m : src->mutations) {
        if (m.is_valid()||m.get_all_major_allele()!=m.get_mut_one_hot()) {
            mutations.emplace_back(m,0, m.get_all_major_allele());
        }
    }
}
struct Profitable_Move_Comparator {
    bool operator()(const Profitable_Moves_ptr_t& lhs,const Profitable_Moves_ptr_t& rhs)const {
        return lhs->radius_left>rhs->radius_left;
    }
};
int individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA,output_t& out,bool do_drift
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    ,MAT::Tree* tree
#endif
                   ) {
    MAT::Node *root = src->parent;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection root_mutations_altered;
    Mutation_Count_Change_Collection new_alter_mutations;
    int parsimony_score_change = do_drift?-1:0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<Mutation_Count_Change_Collection> debug;
#endif
    MAT::Node * last_src_branch_node=src;
    if(src->parent->children.size()<=1) {
        return 0;
    }
    if(root!=LCA) {
        get_parent_altered_remove(root_mutations_altered, src, parsimony_score_change);
        last_src_branch_node=src->parent;
        merge_mutation_src_to_LCA(root, mutations);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        debug.push_back(root_mutations_altered);
#endif
        root=root->parent;
    }


    while (root&&root!=LCA) {
        get_intermediate_nodes_mutations(root,root_mutations_altered, new_alter_mutations, parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        debug.push_back(new_alter_mutations);
#endif
        last_src_branch_node=root;
        merge_mutation_src_to_LCA(root, mutations);
        root_mutations_altered = std::move(new_alter_mutations);
        // unresolved_ties_from_LCA = std::move(new_tie);
        root = root->parent;
    }
    std::vector<MAT::Node*> dst_path;
    root=dst;
    while (root!=LCA) {
        dst_path.push_back(root);
        root=root->parent;
    }
    for (int i=dst_path.size()-1; i>0; i--) {
        Mutation_Count_Change_Collection child_mutations;
        merge_mutation_LCA_to_dst(dst_path[i], mutations, child_mutations);
        mutations=std::move(child_mutations);
    }
    if (dst==LCA) {
        return check_move_profitable_LCA(src, LCA, mutations, root_mutations_altered,
                                         parsimony_score_change,last_src_branch_node,out,1
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                         ,debug,tree
#endif
                                        );;
    } else {
        return check_move_profitable_dst_not_LCA(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change,out,1
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                , debug,tree
#endif
                                                );
    }
}