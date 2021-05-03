#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>
int get_parsimmony_score_only(MAT::Node* src, MAT::Node* dst);


/*
static int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    // const New_Tie_Collection_t &unresolved_ties_from_LCA,
    int parsimony_score_change,output_t& output,
    const std::vector<MAT::Node *> node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> debug_from_src
#endif
) {
    int parsimony_score_change__=check_move_profitable_internal(src, dst, LCA,mutations,
    root_mutations_altered,parsimony_score_change,output,node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,debug_from_src
#endif
    );
    int ref=individual_move(src,dst,LCA);
    assert(parsimony_score_change__==ref);
    return parsimony_score_change__;
}*/
int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    // const New_Tie_Collection_t &unresolved_ties_from_LCA,
    int parsimony_score_change, output_t &output,
    const std::vector<MAT::Node *> node_stack_from_src, bool have_shared,
    const Mutation_Count_Change_Collection &dst_unique_mutations
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> debug_from_src
#endif
);
static void
merge_mutation_LCA_to_dst(MAT::Node *child,
                          const Mutation_Count_Change_Collection &mutations,
                          Mutation_Count_Change_Collection &child_mutations) {
    auto child_mutation_iter = child->mutations.begin();
    auto child_mutation_end = child->mutations.end();
    for (const Mutation_Count_Change &m : mutations) {
        while (child_mutation_iter != child_mutation_end &&
               !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }
        while (child_mutation_iter != child_mutation_end &&
               *child_mutation_iter < m) {
            child_mutations.push_back(*child_mutation_iter);
            child_mutations.back().set_change(
                0, child_mutation_iter->get_par_one_hot(),child_mutation_iter->get_par_one_hot());
            child_mutations.back().set_par_nuc(
                child_mutation_iter->get_mut_one_hot());
            child_mutation_iter++;
        }
        if (child_mutation_iter != child_mutation_end &&
            child_mutation_iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_allele = child_mutation_iter->get_mut_one_hot();
                child_mutations.emplace_back(m);
                child_mutations.back().set_par_nuc(new_par_allele);
            child_mutation_iter++;
        } else {
            child_mutations.emplace_back(m);
        }
    }
    while (child_mutation_iter != child_mutation_end) {
        child_mutations.push_back(*child_mutation_iter);
        child_mutations.back().set_change(
            0, child_mutation_iter->get_par_one_hot(),child_mutation_iter->get_par_one_hot());
        child_mutations.back().set_par_nuc(
            child_mutation_iter->get_mut_one_hot());
        child_mutation_iter++;
    }
}
static void search_subtree_not_LCA(
    MAT::Node *root, MAT::Node *child_exclude, MAT::Node *LCA, MAT::Node *src,
    Mutation_Count_Change_Collection mutations,
    output_t &profitable_moves,
    const Mutation_Count_Change_Collection &root_mutations_altered,int radius,
    int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<state_change_hist_dbg> &debug_from_src
#endif
) {
    if(!(root->children.size()==1&&root->children[0]->is_leaf())&&radius>0){
    for (MAT::Node *child : root->children) {
        if (child != child_exclude) {
            Mutation_Count_Change_Collection child_mutations;
            merge_mutation_LCA_to_dst(child, mutations, child_mutations);
            search_subtree_not_LCA(child, nullptr, LCA, src, child_mutations,
                           profitable_moves, root_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src
#endif
            );
        }
    }}
        auto score_change = check_move_profitable(
            src, root, LCA, mutations, root_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src,false,Mutation_Count_Change_Collection()
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src
#endif
        );
}

void search_subtree(
    MAT::Node *LCA, MAT::Node *src,Mutation_Count_Change_Collection mutations,output_t &profitable_moves,const Mutation_Count_Change_Collection &src_branch_mutations_altered,int radius,int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src, Mutation_Count_Change_Collection& node_unique_mutation, bool try_place_at_LCA,const Mutation_Count_Change_Collection &root_mutation_altered
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<state_change_hist_dbg> &debug_from_src
#endif
) {
    if (try_place_at_LCA) {
        auto score_change = check_move_profitable(
            src, LCA, LCA, mutations, src_branch_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src,true,node_unique_mutation
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src
#endif
        );
    }
    if(!(LCA->children.size()==1&&LCA->children[0]->is_leaf())&&radius>0){
    for (MAT::Node *child : LCA->children) {
        if (child != src) {
            Mutation_Count_Change_Collection child_mutations;
            merge_mutation_LCA_to_dst(child, mutations, child_mutations);
            search_subtree_not_LCA(child, nullptr, LCA, src, child_mutations,
                           profitable_moves, src_branch_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src
#endif
            );
        }
    }}

}


static bool
merge_mutation_src_to_LCA(const MAT::Node *root,
                          Mutation_Count_Change_Collection &mutations,Mutation_Count_Change_Collection &dst_unique_mutation) {
    Mutation_Count_Change_Collection merged_mutations;
    bool have_shared=false;
    merged_mutations.reserve(mutations.size() + root->mutations.size());
    dst_unique_mutation.reserve(mutations.size());
    auto iter = mutations.begin();
    auto end = mutations.end();
    for (const auto &m : root->mutations) {
        if (!m.is_valid()) {
            continue;
        }
        while (iter != end && iter->get_position() < m.get_position()) {
            merged_mutations.push_back(*iter);
            dst_unique_mutation.push_back(*iter);
            iter++;
        }
        if (iter != end && iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_nuc = m.get_mut_one_hot();
            if (new_par_nuc != iter->get_incremented()) {
                merged_mutations.push_back(*iter);
                merged_mutations.back().set_par_nuc(m.get_par_one_hot());
                dst_unique_mutation.push_back(*iter);
            }else {
                have_shared=true;
            }
            iter++;
        } else {
            merged_mutations.emplace_back(m);
            merged_mutations.back().set_change(0, m.get_mut_one_hot(), m.get_mut_one_hot());
        }
    }
    while (iter!=end) {
        merged_mutations.push_back(*iter);
        dst_unique_mutation.push_back(*iter);
        iter++;
    }
    mutations = std::move(merged_mutations);
    return have_shared;
}
static void init_mutation_change(MAT::Node* src, Mutation_Count_Change_Collection& mutations){
    mutations.reserve(src->mutations.size());
    for (const auto &m : src->mutations) {
        if (m.is_valid()||m.get_all_major_allele()!=m.get_mut_one_hot()) {
            mutations.emplace_back(m);
            mutations.back().set_change(0, m.get_all_major_allele(),m.get_all_major_allele());
        }
    }
}
void find_profitable_moves(MAT::Node *src, output_t &out,int radius) {
    if (src->is_root()) {
        return;
    }
    if(src->parent->children.size()<=1){
        return;
    }
    MAT::Node *root = src->parent;
    MAT::Node *exclude = src;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection alter_mutations;
    Mutation_Count_Change_Collection new_alter_mutations;
    // New_Tie_Collection_t unresolved_ties_from_LCA;
    // New_Tie_Collection_t new_tie;
    std::vector<MAT::Node *> node_stack;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<state_change_hist_dbg> debug;
#endif
    Mutation_Count_Change_Collection enroute_unique;
    // Do searching on the current parent of src directly
    {
        search_subtree(root, src, mutations, out,
                       Mutation_Count_Change_Collection(),radius, 0,
                       std::vector<MAT::Node *>({}),enroute_unique,false,Mutation_Count_Change_Collection()
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,std::vector<state_change_hist_dbg>()
                       #endif
    );
    }

    get_parent_altered_remove(alter_mutations, src, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                              ,
                              debug, node_stack
#endif
    );

    node_stack.push_back(src->parent);
    bool have_shared=merge_mutation_src_to_LCA(root, mutations,enroute_unique);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(std::max(src->mutations.size(), src->parent->mutations.size()) <=
           debug.size());
    size_t old_debug_size = debug.size();
#endif
    exclude=root;
    root=root->parent;
    while (root&&radius) {
        search_subtree(root, src, mutations, out,
                       alter_mutations,radius, parsimony_score_change,
                       std::vector<MAT::Node *>(node_stack),enroute_unique,have_shared
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,
                       std::vector<state_change_hist_dbg>(debug)
#endif
        );

        get_intermediate_nodes_mutations(
            root, alter_mutations, new_alter_mutations, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,
            debug, node_stack
#endif
        );
        node_stack.push_back(root);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        assert(debug.size() == old_debug_size);
#endif
        // mutations.merge(root->mutations,
        // MAT::Mutations_Collection::KEEP_SELF);
        have_shared=merge_mutation_src_to_LCA(root, mutations,enroute_unique);
        alter_mutations = std::move(new_alter_mutations);
        // unresolved_ties_from_LCA = std::move(new_tie);
        exclude = root;
        root = root->parent;
        radius--;
    }
}
int individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA){
MAT::Node *root = src->parent;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection root_mutations_altered;
    Mutation_Count_Change_Collection new_alter_mutations;
    // New_Tie_Collection_t unresolved_ties_from_LCA;
    // New_Tie_Collection_t new_tie;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<state_change_hist_dbg> debug_from_src;
#endif
    std::vector<MAT::Node *> node_stack_from_src;
    size_t old_debug_size;
    if(src->parent->children.size()<=1){
        return 0;
    }
    if(root!=LCA){
    get_parent_altered_remove(root_mutations_altered, src, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                              ,debug_from_src, node_stack_from_src
#endif
    );
    node_stack_from_src.push_back(src->parent);
    #ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(std::max(src->mutations.size(), src->parent->mutations.size()) <=
           debug_from_src.size());
    old_debug_size = debug_from_src.size();
    merge_mutation_src_to_LCA(root, mutations);
    root=root->parent;
#endif
    }


    while (root&&root!=LCA) {
        get_intermediate_nodes_mutations(
            root, root_mutations_altered, new_alter_mutations, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,
            debug_from_src, node_stack_from_src
#endif
        );
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        assert(debug_from_src.size() == old_debug_size);
#endif
        node_stack_from_src.push_back(root);
        // mutations.merge(root->mutations,
        // MAT::Mutations_Collection::KEEP_SELF);
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
    for (int i=dst_path.size()-1; i>=0; i--) {
        Mutation_Count_Change_Collection child_mutations;
        merge_mutation_LCA_to_dst(dst_path[i], mutations, child_mutations);
        mutations=std::move(child_mutations);
    }
    output_t out;
    return check_move_profitable(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change,out, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    , debug_from_src
#endif
    );
}