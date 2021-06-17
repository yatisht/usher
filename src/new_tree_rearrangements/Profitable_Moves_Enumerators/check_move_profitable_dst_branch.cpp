#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"

//Push all nodes from cur to LCA onto dst_stack
static void add_remaining_dst_to_LCA_nodes(MAT::Node *cur, const MAT::Node *LCA,
                         std::vector<MAT::Node *> &dst_stack
) {
    while (cur != LCA) {
        assert(dst_stack.empty()||dst_stack.back()!=cur);
        dst_stack.push_back(cur);
        cur = cur->parent;
    }
}

/**
 * @brief Calculate Parsimony score change contributed by inserting src between dst and parent of dst, along the path to LCA
 * @param LCA 
 * @param mutations Assuming state of all nodes stay the same, the mutations on the edge from dst to src to perserve the state of moved src.
 * @param[inout] parsimony_score_change 
 * @param node_stack_from_dst
 * @param this_node dst node 
 * @param parent_added the change in major alleles of the immediate children of LCA on path from dst to LCA
 * @param debug_from_dst First element is the Fitch set change from the original dst to the intermediate node formed with children dst and src, and then its parent, until there is no change to fitch set.
 * @return Whether it is more profitable to insert on the edge between this_node->parent and this_node or this_node->parent->parent and this_node->parent 
 */
bool dst_branch(const MAT::Node *LCA,
           const Mutation_Count_Change_Collection &mutations,
           int &parsimony_score_change,
           std::vector<MAT::Node *> &node_stack_from_dst, MAT::Node *this_node,
           Mutation_Count_Change_Collection &parent_added
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
           ,std::vector<Mutation_Count_Change_Collection> &debug_from_dst
#endif
           ) {
    //Change in major alleles set of this_node, needed to update major alleles set of parent of this_node
    //Push dst to dst->LCA node stack
    node_stack_from_dst.push_back(this_node);
    //Check whether the dst have any mutation not shared by src, if not, stop here, as moving to a child of dst will be at least as profitable, if there is any. 
    //Also register parsimony score contributed by dst->src edge, and major allele state change of dst.
    if(!get_parsimony_score_change_from_add(this_node, mutations, parent_added,
                                         parsimony_score_change
    )&&(!this_node->is_leaf())){return false;}
    this_node = this_node->parent;

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    debug_from_dst.push_back(parent_added);
#endif

    //Going up until LCA to calculate change in major allele state of nodes from dst to LCA, and number of mutations on corresponding edges
    while (this_node != LCA) {
        Mutation_Count_Change_Collection parent_of_parent_added;
        get_intermediate_nodes_mutations(
                this_node, node_stack_from_dst.back(), parent_added,
                parent_of_parent_added, parsimony_score_change
            );

        //add current node dst->LCA node stack
        assert(node_stack_from_dst.empty()||node_stack_from_dst.back()!=this_node);
        node_stack_from_dst.push_back(this_node);
        //In next iteration, major allele set change of this_node become the change in count of major allele state among children of its parent 
        parent_added = std::move(parent_of_parent_added);
        //If the major allele state of this_node didn't change, nor should its ancestors, so no need to go up anymore, just complete the node stack
        if (parent_added.empty()) {
            add_remaining_dst_to_LCA_nodes(this_node->parent, LCA, node_stack_from_dst);
            break;
        }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        debug_from_dst.push_back(parent_added);
#endif
        this_node = this_node->parent;
    }
    return true;
}