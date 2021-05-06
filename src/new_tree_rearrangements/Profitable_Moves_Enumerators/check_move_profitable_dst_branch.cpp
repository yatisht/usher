#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"

static void update_debug(MAT::Node *cur, const MAT::Node *LCA,
                         std::vector<MAT::Node *> &dst_stack
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                         ,
                         std::vector<state_change_hist_dbg> &debug
#endif
) {
    while (cur != LCA) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        auto cur_mut_iter = cur->mutations.begin();
        auto cur_mut_end = cur->mutations.end();
        for (state_change_hist_dbg &ele : debug) {
            while (cur_mut_iter != cur_mut_end &&
                   cur_mut_iter->get_position() < ele.position) {
                cur_mut_iter++;
            }
            if (cur_mut_iter != cur_mut_end &&
                cur_mut_iter->get_position() == ele.position) {
                ele.count_change.push_back(Mutation_Count_Change());
                ele.mutation_score_change.push_back(0);
                ele.par_nuc.push_back(cur_mut_iter->get_par_one_hot());
                ele.major_allele.push_back(
                    cur_mut_iter->get_all_major_allele());
            } else {
                ele.count_change.push_back(Mutation_Count_Change());
                ele.mutation_score_change.push_back(0);
                ele.par_nuc.push_back(ele.par_nuc.back());
                ele.major_allele.push_back(ele.par_nuc.back());
            }
        }
#endif
        assert(dst_stack.empty()||dst_stack.back()!=cur);
        dst_stack.push_back(cur);
        cur = cur->parent;
    }
}


bool dst_branch(const MAT::Node *LCA,
           const Mutation_Count_Change_Collection &mutations,
           int &parsimony_score_change,
           std::vector<MAT::Node *> &node_stack_from_dst, MAT::Node *this_node,
           Mutation_Count_Change_Collection &parent_added
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
           ,std::vector<state_change_hist_dbg> &debug_from_dst
#endif
           ) {

    Mutation_Count_Change_Collection parent_of_parent_added;
    node_stack_from_dst.push_back(this_node);
    if(!get_parsimony_score_change_from_add(this_node, mutations, parent_added,
                                         parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                        ,
                                        debug_from_dst, node_stack_from_dst
#endif
    )){return false;}
    this_node = this_node->parent;
    while (this_node != LCA) {
        if (this_node->children.size() == 2) {
            get_two_child_intermediate_node_mutations(
                this_node, node_stack_from_dst.back(), parent_added,
                parent_of_parent_added, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,
                debug_from_dst, node_stack_from_dst
#endif
            );
        } else {
            get_intermediate_nodes_mutations(this_node, parent_added,
                                             parent_of_parent_added,
                                             parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                             ,
                                             debug_from_dst, node_stack_from_dst
#endif
            );
        }
        assert(node_stack_from_dst.empty()||node_stack_from_dst.back()!=this_node);
        node_stack_from_dst.push_back(this_node);
        parent_added = std::move(parent_of_parent_added);
        if (parent_added.empty()) {
            update_debug(this_node->parent, LCA, node_stack_from_dst
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                         ,
                         debug_from_dst
#endif
            );
            break;
        }

        this_node = this_node->parent;
    }
    return true;
}