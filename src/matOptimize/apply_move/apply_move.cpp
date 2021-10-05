#include "apply_move.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_set>
#include <vector>
std::vector<std::string> changed_nodes;
//Fix state of src node when its parent state changes to make it consistent with parent state
//, and if a state change is needed, add to altered nodes list
void clean_up_src_states(MAT::Node *src,
                         std::vector<Altered_Node_t> &out) {
    MAT::Mutations_Collection &in = src->mutations;
    //bool have_change = false;
    State_Change_Collection changed_state({state_change(NEW_SRC_MARK)});
    for (auto &mut : in) {
        nuc_one_hot new_state =
            mut.get_par_one_hot() & mut.get_all_major_allele();
        nuc_one_hot old_state = mut.get_mut_one_hot();
        //Change state if it can now follow parent state, and the old state is inconsistent
        if (new_state && old_state != new_state) {
            changed_nodes.push_back(src->identifier);
            //have_change = true;
            mut.set_mut_one_hot(new_state);
            changed_state.emplace_back(mut, old_state);
        }
    }
    //Clean up irrelevent mutations whose major allele agree with parent allele and have no sensitive alleles (allele count one less than major allele count)
    in.mutations.erase(std::remove_if(in.mutations.begin(), in.mutations.end(),
    [](const MAT::Mutation &mut) {
        return (
                   mut.get_par_one_hot() ==
                   mut.get_all_major_allele() &&
                   (!mut.get_boundary1_one_hot()));
    }),
    in.mutations.end());
    if (changed_state.size()>1) {
        out.emplace_back(src);
        out.back().changed_states = std::move(changed_state);
    }
}
#ifdef CHECK_STATE_REASSIGN
void compare_mutations(MAT::Node *old_nodes, MAT::Node *new_nodes) {
    for (size_t mut_idx = 0; mut_idx < old_nodes->mutations.size(); mut_idx++) {
        if (!(old_nodes->mutations[mut_idx] == new_nodes->mutations[mut_idx])) {
            fprintf(stderr, "%d\n",mut_idx);
            assert(false);
        }
        assert(old_nodes->mutations[mut_idx] == new_nodes->mutations[mut_idx]);
        if (old_nodes->children.size() == 2) {
            assert(old_nodes->mutations[mut_idx].get_left_child_state() ==
                   new_nodes->mutations[mut_idx].get_left_child_state());
            assert(old_nodes->mutations[mut_idx].get_right_child_state() ==
                   new_nodes->mutations[mut_idx].get_right_child_state());
        }
    }
}
#endif
void check_major_state(MAT::Node *node, const MAT::Tree &new_tree);
//Called when debugging to locate a node
const Altered_Node_t* find_altered_node(char* node_name, const std::vector<Altered_Node_t>& to_find) {
    for(const auto& temp:to_find) {
        if (temp.altered_node->identifier==node_name) {
            return &temp;
        }
    }
    return nullptr;
}

#ifdef CHECK_STATE_REASSIGN
void compare_mutation_tree(MAT::Tree &t,MAT::Tree &new_tree) {
    std::vector<MAT::Node *> old_nodes = t.breadth_first_expansion();
    std::vector<MAT::Node *> new_nodes = new_tree.breadth_first_expansion();
    for (int i = old_nodes.size() - 1; i >= 0; i--) {
        compare_mutations(old_nodes[i], new_nodes[i]);
    }
    new_tree.delete_nodes();
}
#endif
/**
 * @brief Apply all the non-conflicting moves
 * @param all_moves all non-conflicting moves to apply
 * @param t Tree where the nodes are
 * @param bfs_ordered_nodes all nodes in bfs order
 * @param to_filter nodes to search in the next round, with removed nodes filtered out
 */
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t
#ifdef CHECK_STATE_REASSIGN
                 ,
                 const Original_State_t &original_state
#endif
                ) {
    std::vector<MAT::Node *> altered_node; //nodes whose children are altered, need to do FS backward pass from them
    std::vector<Altered_Node_t> forward_pass_altered_nodes; //Their state have changed, need to make their children follow parent
    std::unordered_set<size_t> deleted_node_ptrs;//pointers of deleted nodes, delay deletion until the end
    std::vector<MAT::Node *> nodes_to_clean;//The src nodes and sibling whose parent have been condensed, whose parent state may change due to branch splitting
    //std::unordered_set<std::string> samples;
#ifdef CHECK_STATE_REASSIGN
    /*for(const auto& s:original_state){
        samples.insert(s.first);
    }*/
#ifndef SINGLE_THREAD_TEST
    FILE* log=fopen("moves", "w");
#endif
#endif
    for (const auto &move : all_moves) {
        //skip moves whose end points have been deleted
        if (deleted_node_ptrs.count((size_t)move->src) ||
                deleted_node_ptrs.count((size_t)move->get_dst())||move->src->parent==move->get_dst()) {
            continue;
        }
#ifdef CHECK_STATE_REASSIGN
#ifdef SINGLE_THREAD_TEST
        fprintf(stderr, "%s\tto\t%s\n",move->src->identifier.c_str(),move->get_dst()->identifier.c_str());
#else
        fprintf(log, "%s\tto\t%s\n",move->src->identifier.c_str(),move->get_dst()->identifier.c_str());
        fflush(log);
#endif
#endif
        //fprintf(stderr, "applying move %zu to %zu \n",move->src->dfs_index,move->dst->dfs_index);
        //Prelimary move perserving state of all nodes
        move_node(move->src, move->get_dst(), altered_node, t,
                  deleted_node_ptrs, nodes_to_clean);
        /*std::unordered_set<std::string> to_check(samples);
        for(auto node:dfs){
            to_check.erase(node->identifier);
        }
        assert(to_check.empty());*/
    }
    //need to redo dfs, as the Fitch sankoff patching need to process nodes in dfs order to have the
    // major allele set of children of a node fully updated before assigning major allele to this node
    std::vector<MAT::Node*> dfs=t.depth_first_expansion();
    //filter deleted nodes among nodes to do fitch sankoff
    auto end = std::remove_if(altered_node.begin(), altered_node.end(),
    [&deleted_node_ptrs](MAT::Node *node) {
        return deleted_node_ptrs.count((size_t)node);
    });
    altered_node.erase(end, altered_node.end());
    nodes_to_clean.erase(std::remove_if(nodes_to_clean.begin(), nodes_to_clean.end(),
    [&deleted_node_ptrs](MAT::Node *node) {
        return deleted_node_ptrs.count((size_t)node);
    }),nodes_to_clean.end());
#ifdef CHECK_STATE_REASSIGN

    //MAT::Tree new_tree;
    //new_tree.load_detatiled_mutations("After_reassign.pb");
    MAT::Tree new_tree = reassign_state_full(t);
    new_tree.save_detailed_mutations("After_reassign.pb");
#endif
    //backward pass, patch major allele set
    if (!altered_node.empty()) {
        reassign_backward_pass(altered_node, forward_pass_altered_nodes
#ifdef CHECK_STATE_REASSIGN
                               ,
                               new_tree
#endif
                              );
    }
    for (const auto node : nodes_to_clean) {
        clean_up_src_states(node, forward_pass_altered_nodes);
    }
    //forward pass patch state from parent state and calibrated major allele set
    if (!forward_pass_altered_nodes.empty()) {
        forward_pass(forward_pass_altered_nodes
#ifdef CHECK_STATE_REASSIGN
                     ,
                     new_tree
#endif
                    );
    }
#ifdef CHECK_STATE_REASSIGN
    compare_mutation_tree(t, new_tree);
#endif
    //filter out deleted nodes from nodes to search in the next round
    //delayed deletion of removed nodes (these are identified by there memory location, if freed to early, it can be reused)
    for(auto node:deleted_node_ptrs) {
        t.all_nodes.erase(((MAT::Node*) node)->identifier);
        delete ((MAT::Node*) node);
    }
#ifdef CHECK_STATE_REASSIGN
#ifndef SINGLE_THREAD_TEST
    fclose(log);
#endif
#endif
}