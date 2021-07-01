#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>
//defined and documented at check_move_profitable.cpp
int check_move_profitable_dst_not_LCA(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const range<Mutation_Count_Change>  &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change, output_t &output,
    const std::vector<MAT::Node *>& node_stack_from_src,int radius_left
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> debug_from_src,
    const MAT::Tree* tree
#endif
);
int check_move_profitable_LCA(
    MAT::Node *src, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change,
    const std::vector<MAT::Node *> &node_stack_from_src, output_t &output,int radius_left
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> &debug_above_LCA,
    const MAT::Tree* tree
#endif
);
template<typename output_type>
static void
merge_mutation_LCA_to_dst(MAT::Node *child,
                          range<Mutation_Count_Change> mutations,
                          output_type &child_mutations) {
    //probably a gross overestimate, as we don't merge invalid mutations
    child_mutations.reserve(mutations.size()+child->mutations.size());
    auto child_mutation_iter = child->mutations.begin();
    auto child_mutation_end = child->mutations.end();
    for ( ; mutations;mutations++) {
        const Mutation_Count_Change &m=*mutations;
        //skip invalid mutations, src don't have these allele distribution among its children
        while (child_mutation_iter != child_mutation_end &&
               !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }
        while (child_mutation_iter != child_mutation_end &&
               *child_mutation_iter < m) {
            //newly encountered mutation to merge, add as back mutation
            child_mutations.emplace_back(*child_mutation_iter,0,child_mutation_iter->get_par_one_hot(),child_mutation_iter->get_par_one_hot(),true);
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
                child_mutations.emplace_back(m);
                child_mutations.back().set_par_nuc(new_par_allele);
            child_mutation_iter++;
        } else {
            //copy over other undisturbed mutations
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
//basically the same as such subtree, and is called by it, but it explore subtree rooted at 'root' for dst
static void search_subtree_not_LCA(
    MAT::Node *root, MAT::Node *LCA, MAT::Node *src,
    const range<Mutation_Count_Change>& mutations,
    output_t &profitable_moves,
    const Mutation_Count_Change_Collection &root_mutations_altered,int radius,
    int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src
    ,stack_allocator<Mutation_Count_Change>& allocator
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<Mutation_Count_Change_Collection> &debug_from_src,
    MAT::Tree* tree
#endif
) {
    if(!(root->children.size()==1&&root->children[0]->is_leaf())&&radius>0){
        Mutation_Count_Change_Collection_FIFO child_mutations(allocator);
        merge_mutation_LCA_to_dst(root, mutations, child_mutations);
        //as root is not LCA, don't need to worry about a node is moved under its original parent again
        //Pre-order, search when mutation vector of dst is slightly more likely to be hot
        check_move_profitable_dst_not_LCA(
            src, root, LCA, mutations, root_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src,radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src,tree
#endif
        );
        //Haven't exhausted hops, recurse down to explore more dsts
    for (MAT::Node *child : root->children) {
            search_subtree_not_LCA(child, LCA, src, child_mutations,
                           profitable_moves, root_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src,allocator
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src,tree
#endif
            );
    }}

}
/**
 * @brief Find most profitable dst to move src to, under the subtree rooted at LCA
 * @param LCA
 * @param src
 * @param mutations  Mutations moved src needs to maintain its state if moved as child of LCA
 * @param profitable_moves output
 * @param src_branch_mutations_altered change in Fitch set of src_branch_node from removal of src 
 * @param radius number of hops left to go down the tree
 * @param parsimony_score_change_from_removal
 * @return (void)
 */
void search_subtree(
    MAT::Node *LCA, MAT::Node *src,
    const Mutation_Count_Change_Collection &mutations,
    output_t &profitable_moves,
    const Mutation_Count_Change_Collection &src_branch_mutations_altered,
    int radius, int parsimony_score_change_from_removal,
    const std::vector<MAT::Node *> &node_stack_from_src,
    stack_allocator<Mutation_Count_Change>& allocator
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> &debug_from_src,
    MAT::Tree *tree
#endif
) {
    if (src->parent!=LCA) {
            //moving src to LCA, if src's parent was not originally LCA 
            check_move_profitable_LCA(
            src, LCA, mutations, src_branch_mutations_altered,
            parsimony_score_change_from_removal, node_stack_from_src,profitable_moves,radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src,tree
#endif
        );
    }
    //Avoid moving a node to its descendents.
    MAT::Node* exclude=node_stack_from_src.empty()?src:node_stack_from_src.back();
    if(!(LCA->children.size()==1&&LCA->children[0]->is_leaf())&&radius>0){
    for (MAT::Node *child : LCA->children) {
        if (child !=exclude) {
            search_subtree_not_LCA(child, LCA, src,mutations,
                           profitable_moves, src_branch_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src,
                           allocator
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
            merged_mutations.emplace_back(m,0, m.get_mut_one_hot(), m.get_mut_one_hot(),true);
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
static void init_mutation_change(MAT::Node* src, Mutation_Count_Change_Collection& mutations){
    mutations.reserve(src->mutations.size());
    for (const auto &m : src->mutations) {
        if (m.is_valid()||m.get_all_major_allele()!=m.get_mut_one_hot()) {
            mutations.emplace_back(m,0, m.get_all_major_allele(),m.get_all_major_allele());
        }
    }
}
void find_profitable_moves(MAT::Node *src, output_t &out,int radius,
    stack_allocator<Mutation_Count_Change>& allocator
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
                       node_stack,allocator
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
                       node_stack,allocator
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
    for (int i=dst_path.size()-1; i>0; i--) {
        Mutation_Count_Change_Collection child_mutations;
        merge_mutation_LCA_to_dst(dst_path[i], mutations, child_mutations);
        mutations=std::move(child_mutations);
    }
    if (dst==LCA) {
        return check_move_profitable_LCA(src, LCA, mutations, root_mutations_altered,
            parsimony_score_change, node_stack_from_src,out,1
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug,tree
#endif
        );;
    }else{
    return check_move_profitable_dst_not_LCA(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change,out, node_stack_from_src,1
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    , debug,tree
#endif
    );}
}