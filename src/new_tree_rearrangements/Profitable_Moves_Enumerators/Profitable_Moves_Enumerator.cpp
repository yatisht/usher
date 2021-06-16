#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>

int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change, output_t &output,
    const std::vector<MAT::Node *> node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> debug_from_src,
    const MAT::Tree* tree
#endif
);
static void
merge_mutation_LCA_to_dst(MAT::Node *child,
                          const Mutation_Count_Change_Collection &mutations,
                          Mutation_Count_Change_Collection &child_mutations) {
    child_mutations.reserve(mutations.size()+child->mutations.size());
    auto child_mutation_iter = child->mutations.begin();
    auto child_mutation_end = child->mutations.end();
    for (const Mutation_Count_Change &m : mutations) {
        while (child_mutation_iter != child_mutation_end &&
               !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }
        while (child_mutation_iter != child_mutation_end &&
               *child_mutation_iter < m) {
            child_mutations.emplace_back(*child_mutation_iter,0,0,0,true);
            child_mutations.back().set_change(
                0, child_mutation_iter->get_par_one_hot(),child_mutation_iter->get_par_one_hot(),true);
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
        child_mutations.emplace_back(*child_mutation_iter,0,0,0,true);
        child_mutations.back().set_change(
            0, child_mutation_iter->get_par_one_hot(),child_mutation_iter->get_par_one_hot(),true);
        child_mutations.back().set_par_nuc(
            child_mutation_iter->get_mut_one_hot());
        child_mutation_iter++;
    }
}
static void search_subtree_not_LCA(
    MAT::Node *root, MAT::Node *LCA, MAT::Node *src,
    const Mutation_Count_Change_Collection& mutations,
    output_t &profitable_moves,
    const Mutation_Count_Change_Collection &root_mutations_altered,int radius,
    int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<Mutation_Count_Change_Collection> &debug_from_src,
    MAT::Tree* tree
#endif
) {
    if(!(root->children.size()==1&&root->children[0]->is_leaf())&&radius>0){
        Mutation_Count_Change_Collection child_mutations;
        merge_mutation_LCA_to_dst(root, mutations, child_mutations);
    for (MAT::Node *child : root->children) {
            search_subtree_not_LCA(child, LCA, src, child_mutations,
                           profitable_moves, root_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src,tree
#endif
            );
    }}
            check_move_profitable(
            src, root, LCA, mutations, root_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src,tree
#endif
        );
}

void search_subtree(
    MAT::Node *LCA, MAT::Node *src,Mutation_Count_Change_Collection mutations,output_t &profitable_moves,const Mutation_Count_Change_Collection &src_branch_mutations_altered,int radius,int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<Mutation_Count_Change_Collection> &debug_from_src,MAT::Tree* tree
#endif
) {
    if (src->parent!=LCA) {
            check_move_profitable(
            src, LCA, LCA, mutations, src_branch_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src,tree
#endif
        );
    }
    MAT::Node* exclude=node_stack_from_src.empty()?src:node_stack_from_src.back();
    if(!(LCA->children.size()==1&&LCA->children[0]->is_leaf())&&radius>0){
    for (MAT::Node *child : LCA->children) {
        if (child !=exclude) {
            search_subtree_not_LCA(child, LCA, src,mutations,
                           profitable_moves, src_branch_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src,tree
#endif
            );
        }
    }}

}


static void
merge_mutation_src_to_LCA(const MAT::Node *root,
                          Mutation_Count_Change_Collection &mutations) {
    Mutation_Count_Change_Collection merged_mutations;
    merged_mutations.reserve(mutations.size() + root->mutations.size());
    auto iter = mutations.begin();
    auto end = mutations.end();
    for (const auto &m : root->mutations) {
        if (!m.is_valid()) {
            continue;
        }
        while (iter != end && iter->get_position() < m.get_position()) {
            merged_mutations.push_back(*iter);
            iter++;
        }
        if (iter != end && iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_nuc = m.get_par_one_hot();
            if (new_par_nuc != iter->get_incremented()) {
                merged_mutations.push_back(*iter);
                merged_mutations.back().set_par_nuc(m.get_par_one_hot());
            }
            iter++;
        } else {
            merged_mutations.emplace_back(m,0, m.get_mut_one_hot(), m.get_mut_one_hot(),true);
        }
    }
    while (iter!=end) {
        merged_mutations.push_back(*iter);
        iter++;
    }
    mutations = std::move(merged_mutations);
}
static void init_mutation_change(MAT::Node* src, Mutation_Count_Change_Collection& mutations){
    mutations.reserve(src->mutations.size());
    for (const auto &m : src->mutations) {
        if (m.is_valid()||m.get_all_major_allele()!=m.get_mut_one_hot()) {
            mutations.emplace_back(m,0, m.get_all_major_allele(),m.get_all_major_allele());
        }
    }
}
void find_profitable_moves(MAT::Node *src, output_t &out,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,MAT::Tree* tree
#endif
) {
    if (src->is_root()) {
        return;
    }
    if(src->parent->children.size()<=1){
        return;
    }
    MAT::Node *root = src->parent;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection alter_mutations;
    Mutation_Count_Change_Collection new_alter_mutations;
    std::vector<MAT::Node *> node_stack;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<Mutation_Count_Change_Collection> debug;
#endif
    // Do searching on the current parent of src directly
    {
        search_subtree(root, src, mutations, out,
                       Mutation_Count_Change_Collection(),radius, 0,
                       std::vector<MAT::Node *>({})
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,std::vector<Mutation_Count_Change_Collection>(),tree
                       #endif
    );
    }

    get_parent_altered_remove(alter_mutations, src, parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    debug.push_back(alter_mutations);
#endif
    node_stack.push_back(src->parent);
    merge_mutation_src_to_LCA(root, mutations);
    root=root->parent;
    while (root&&radius) {
        search_subtree(root, src, mutations, out,
                       alter_mutations,radius, parsimony_score_change,
                       std::vector<MAT::Node *>(node_stack)
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,
                       std::vector<Mutation_Count_Change_Collection>(debug),tree
#endif
        );

        get_intermediate_nodes_mutations(root,node_stack.back(),alter_mutations, new_alter_mutations, parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        debug.push_back(new_alter_mutations);
#endif
        node_stack.push_back(root);
        merge_mutation_src_to_LCA(root, mutations);
        alter_mutations = std::move(new_alter_mutations);
        root = root->parent;
        radius--;
    }
}
int individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA,output_t& out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,MAT::Tree* tree
#endif
){
MAT::Node *root = src->parent;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection root_mutations_altered;
    Mutation_Count_Change_Collection new_alter_mutations;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<Mutation_Count_Change_Collection> debug;
#endif
    std::vector<MAT::Node *> node_stack_from_src;
    if(src->parent->children.size()<=1){
        return 0;
    }
    if(root!=LCA){
    get_parent_altered_remove(root_mutations_altered, src, parsimony_score_change);
    node_stack_from_src.push_back(src->parent);
    merge_mutation_src_to_LCA(root, mutations);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    debug.push_back(root_mutations_altered);
#endif
    root=root->parent;
    }


    while (root&&root!=LCA) {
        get_intermediate_nodes_mutations(root,node_stack_from_src.back(),root_mutations_altered, new_alter_mutations, parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        debug.push_back(new_alter_mutations);
#endif
        assert(node_stack_from_src.back()!=root);
        node_stack_from_src.push_back(root);
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
    return check_move_profitable(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change,out, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    , debug,tree
#endif
    );
}