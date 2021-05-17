#include "apply_move.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstdio>
#include <unordered_set>
#include <vector>
void move_node(MAT::Node *src, MAT::Node *dst,
               std::vector<MAT::Node *> &altered_node, MAT::Tree &tree,
               std::unordered_set<size_t> &deleted,
               std::vector<MAT::Node *> &nodes_to_clean
#ifdef CHECK_PRIMARY_MOVE
               ,
               Original_State_t original_state
#endif
);

#ifdef CHECK_PRIMARY_MOVE
MAT::Tree reassign_state_full(MAT::Tree &tree_in);
#endif

void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out);

void forward_pass(std::vector<Altered_Node_t> &in
#ifdef CHECK_STATE_REASSIGN
                  ,
                  MAT::Tree &new_tree
#endif
);
void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
    ,
    MAT::Tree &new_tree
#endif
);
static void clean_up_src_states(MAT::Node *src,
                                std::vector<Altered_Node_t> &out) {
    MAT::Mutations_Collection &in = src->mutations;
    bool have_change = false;
    State_Change_Collection changed_state;
    for (auto &mut : in) {
        nuc_one_hot new_state =
            mut.get_par_one_hot() & mut.get_all_major_allele();
        nuc_one_hot old_state = mut.get_mut_one_hot();
        if (new_state && old_state != new_state) {
            have_change = true;
            mut.set_mut_one_hot(new_state);
            changed_state.emplace_back(mut, old_state);
        }
    }
    in.mutations.erase(std::remove_if(in.mutations.begin(), in.mutations.end(),
                                      [](const MAT::Mutation &mut) {
                                          return (
                                              mut.get_par_one_hot() ==
                                                  mut.get_all_major_allele() &&
                                              (!mut.get_boundary1_one_hot()));
                                      }),
                       in.mutations.end());
    if (!changed_state.empty()) {
        out.emplace_back(src);
        out.back().changed_states = std::move(changed_state);
    }
}
void compare_mutations(MAT::Node *old_nodes, MAT::Node *new_nodes) {
    for (int mut_idx = 0; mut_idx < old_nodes->mutations.size(); mut_idx++) {
        assert(old_nodes->mutations[mut_idx] == new_nodes->mutations[mut_idx]);
        if (old_nodes->children.size() == 2) {
            assert(old_nodes->mutations[mut_idx].get_left_child_state() ==
                   new_nodes->mutations[mut_idx].get_left_child_state());
            assert(old_nodes->mutations[mut_idx].get_right_child_state() ==
                   new_nodes->mutations[mut_idx].get_right_child_state());
        }
    }
}
void check_major_state(MAT::Node *node, const MAT::Tree &new_tree);
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_PRIMARY_MOVE
                 ,
                 Original_State_t original_state
#endif
) {
    std::vector<MAT::Node *> altered_node;
    std::vector<Altered_Node_t> forward_pass_altered_nodes;
    std::unordered_set<size_t> deleted_node_ptrs;
    std::vector<MAT::Node *> nodes_to_clean;
    for (const auto &move : all_moves) {
        if (move->src->identifier == "node_32_condensed_2_leaves" &&
            move->get_dst()->identifier == "338") {
            fputc('a', stderr);
        }
        if (deleted_node_ptrs.count((size_t)move->src) ||
            deleted_node_ptrs.count((size_t)move->get_dst())) {
            continue;
        }
        move_node(move->src, move->get_dst(), altered_node, t,
                  deleted_node_ptrs, nodes_to_clean, original_state);
        t.depth_first_expansion();
    }
    auto end = std::remove_if(altered_node.begin(), altered_node.end(),
                              [&deleted_node_ptrs](MAT::Node *node) {
                                  return deleted_node_ptrs.count((size_t)node);
                              });
    altered_node.erase(end, altered_node.end());
#ifdef CHECK_PRIMARY_MOVE
    std::vector<MAT::Node *> old_nodes = t.breadth_first_expansion();
    for (auto node : old_nodes) {
        assert(node->is_root() || node->is_leaf() || node->children.size() > 1);
        // assert(node->is_root()||node->mutations.size()>0);
    }
    MAT::Tree new_tree = reassign_state_full(t);
#endif
    reassign_backward_pass(altered_node, forward_pass_altered_nodes
#ifdef CHECK_STATE_REASSIGN
                           ,
                           new_tree
#endif
    );
#ifdef CHECK_STATE_REASSIGN
    for (const auto move : all_moves) {
        check_major_state(move->src, new_tree);
        auto p = move->src->parent;
        if (p != move->get_dst()) {
            check_major_state(p, new_tree);
        }
    }
#endif
    for (const auto node : nodes_to_clean) {
        clean_up_src_states(node, forward_pass_altered_nodes);
    }
    if (!forward_pass_altered_nodes.empty()) {
        forward_pass(forward_pass_altered_nodes
#ifdef CHECK_STATE_REASSIGN
                     ,
                     new_tree
#endif
        );
    }
#ifdef CHECK_STATE_REASSIGN
    std::vector<MAT::Node *> new_nodes = new_tree.breadth_first_expansion();
    for (int i = old_nodes.size() - 1; i >= 0; i--) {
        compare_mutations(old_nodes[i], new_nodes[i]);
    }
    new_tree.delete_nodes();
#endif
}
