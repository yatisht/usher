#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "process_each_node.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
static void check_LCA(MAT::Node* src, MAT::Node* dst,MAT::Node* LCA_to_check) {
    std::unordered_set<int> path;
    while (src!=LCA_to_check) {
        path.insert(src->bfs_index);
        src=src->parent;
    }
    while (dst!=LCA_to_check) {
        assert(!path.count(dst->bfs_index));
        dst=dst->parent;
    }
    assert(dst==LCA_to_check);
}
#endif

bool output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,int radius_left) {
    if (parsimony_score_change <= output.score_change) {
        if (parsimony_score_change < output.score_change) {
            output.score_change = parsimony_score_change;
            output.moves->clear();
        }
        Profitable_Moves_ptr_t new_move = Profitable_Moves_ptr_t(new Profitable_Moves);
        new_move->score_change = parsimony_score_change;
        new_move->src = src;
        new_move->dst=dst;
        new_move->LCA = LCA;
        new_move->radius_left=radius_left;
        //assert(dst == LCA || new_move->dst_to_LCA.back()->parent == LCA);
        //assert(src->parent == LCA ||new_move->src_to_LCA.back()->parent == LCA);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        std::unordered_set<int> node_idx_set;
        node_idx_set.insert(src->bfs_index);
        for (auto node : new_move->src_to_LCA) {
            assert(node_idx_set.insert(node->bfs_index).second);
        }
        for (auto node : new_move->dst_to_LCA) {
            assert(node_idx_set.insert(node->bfs_index).second);
        }
#endif
        output.moves->push_back(new_move);
        return true;
    }
    return false;
}
int get_parsimmony_score_only(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,const MAT::Tree* ori_tree);
bool LCA_place_mezzanine(
    const MAT::Node *src_branch_node,
    const Mutation_Count_Change_Collection &dst_mutations,
    const Mutation_Count_Change_Collection &src_branch_node_mutations_altered,
    Mutation_Count_Change_Collection &out, int &parsimony_score_change
) ;
extern unsigned int early_stop_saving;
/**
 * @brief Calculate Parsimony score change of this individual move without applying it
 * @param src
 * @param dst
 * @param LCA
 * @param mutations mutations needed above src node for moving src under parent of dst
 * @param root_mutations_altered fitch set change of src_branch_node from removing src
 * @param parsimony_score_change
 * @param debug_above_LCA
 * @param output Most profitable moves from this src node
 * @return The highest node had their Fitch set changed
 */
int check_move_profitable_dst_not_LCA(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const range<Mutation_Count_Change>  &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change, output_t &output,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> debug_from_src,
    const MAT::Tree* tree
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    check_LCA(src, dst, LCA);
    assert(src->parent!=dst);
    assert(src!=dst);
    std::vector<Mutation_Count_Change_Collection> debug_from_dst;
#endif

    //assert(dst);
    Mutation_Count_Change_Collection dst_added;
    //Going up from dst node to LCA node to adjust state assignment
    dst_branch(LCA, mutations, parsimony_score_change,
                    dst, dst_added,src->parent==LCA?src->mutations.size():0
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                   , debug_from_dst
#endif
                  );
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<Mutation_Count_Change_Collection> debug_above_LCA;
#endif
    Mutation_Count_Change_Collection parent_of_parent_added;
    parent_of_parent_added.reserve(mutations.size()+root_mutations_altered.size());
    //Adjust LCA node and above

    bool is_src_terminal = src->parent == LCA;
    if ((!(root_mutations_altered.empty() && dst_added.empty())) ||
            is_src_terminal ) {
        get_LCA_mutation(LCA,src, root_mutations_altered, dst_added, parent_of_parent_added, parsimony_score_change);
    }
    dst_added.swap(parent_of_parent_added);
    check_parsimony_score_change_above_LCA(LCA->parent, parsimony_score_change, dst_added, parent_of_parent_added);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    //fprintf(stderr, "LCA idx: %zu",LCA_a->bfs_index);
    auto ref_score=get_parsimmony_score_only(src,dst,LCA,tree);
    assert(parsimony_score_change == ref_score);
#endif
    output_result(src, dst, LCA, parsimony_score_change, output,radius);
    return parsimony_score_change;
}

int check_move_profitable_LCA(
    MAT::Node *src, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change,
    const MAT::Node* last_src_branch_node_below_LCA,
    output_t &output,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> &debug_above_LCA,
    const MAT::Tree* tree
#endif
) {
    Mutation_Count_Change_Collection parent_of_parent_added;
    parent_of_parent_added.reserve(mutations.size()+root_mutations_altered.size());
    LCA_place_mezzanine(last_src_branch_node_below_LCA, mutations, root_mutations_altered, parent_of_parent_added, parsimony_score_change);
    //No need to go to parent, the node from branch splitting is the actual LCA
    Mutation_Count_Change_Collection parent_added;
    parent_added.reserve(parent_of_parent_added.size());
    parent_added .swap(parent_of_parent_added);
    //Go up the tree until there is no change to fitch set
    check_parsimony_score_change_above_LCA(LCA, parsimony_score_change, parent_added, parent_of_parent_added);

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    //fprintf(stderr, "LCA idx: %zu",LCA_a->bfs_index);
    auto ref_score=get_parsimmony_score_only(src,LCA,LCA,tree);
    assert(parsimony_score_change == ref_score);
#endif
    std::vector<MAT::Node*> ignored;
    output_result(src, LCA, LCA, parsimony_score_change, output,
                  radius);
    return parsimony_score_change;
}