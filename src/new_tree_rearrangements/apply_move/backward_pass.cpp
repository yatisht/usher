#include "apply_move.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <vector>

bool get_major_allele_polytomy(MAT::Node *node,
                               MAT::Mutations_Collection &new_major_alleles_out,
                               State_Change_Collection &changed_positions);

bool get_major_allele_binary(const MAT::Node *node,
                             MAT::Mutations_Collection &out,
                             State_Change_Collection &changed_states);

struct State_Assign_Comparator {
    bool operator()(const Altered_Node_t &first, const Altered_Node_t &second) {
        return first.altered_node->dfs_index < second.altered_node->dfs_index;
    }
};

class Backward_Pass_Heap {
    std::vector<Altered_Node_t> altered_nodes;
    std::vector<Altered_Node_t> &processed_altered_nodes;
    void increment() {
        if (altered_nodes.size() <= 1) {
            return;
        }
        std::pop_heap(altered_nodes.begin(), altered_nodes.end(),
                      State_Assign_Comparator());
        while (altered_nodes.front().altered_node ==
               altered_nodes.back().altered_node) {
            std::pop_heap(altered_nodes.begin(), altered_nodes.end() - 1,
                          State_Assign_Comparator());
#ifdef CHECK_STATE_REASSIGN
            size_t last_idx = altered_nodes.size() - 1;
            assert(altered_nodes[last_idx].altered_node ==
                   altered_nodes[last_idx - 1].altered_node);
#endif
            altered_nodes.pop_back();
            if (altered_nodes.size() <= 1) {
                return;
            }
        }
    }

  public:
    Backward_Pass_Heap(const std::vector<MAT::Node *> &in,
                       std::vector<Altered_Node_t> &processed_altered_nodes)
        : processed_altered_nodes(processed_altered_nodes) {
        altered_nodes.reserve(in.size());
        for (auto node : in) {
            altered_nodes.emplace_back(node);
        }
        std::make_heap(altered_nodes.begin(), altered_nodes.end(),
                       State_Assign_Comparator());
        increment();
    }

    Altered_Node_t &operator*() { return altered_nodes.back(); }
    Altered_Node_t *operator->() { return &altered_nodes.back(); }

    bool next(bool major_allele_changed) {
        auto &changed = altered_nodes.back();
        if (!changed.changed_states.empty()) {
            processed_altered_nodes.push_back(changed);
        }
        if (major_allele_changed) {
            MAT::Node *parent_node = changed.altered_node->parent;
            if (parent_node) {
                changed.altered_node = parent_node;
                changed.changed_states.clear();
                std::push_heap(altered_nodes.begin(), altered_nodes.end(),
                               State_Assign_Comparator());
            } else {
                altered_nodes.pop_back();
            }
        } else {
            altered_nodes.pop_back();
        }
        if (altered_nodes.empty()) {
            return false;
        }
        increment();
        return true;
    }
};
#ifdef CHECK_STATE_REASSIGN
static void check_changed(const MAT::Mutations_Collection &ori,
                          const MAT::Mutations_Collection &changed_mut,
                          bool changed,
                          const State_Change_Collection &state_changed,bool no_check) {
    auto state_change_iter = state_changed.begin();
    auto iter = ori.begin();
    auto end = ori.end();
    bool checked = false;
    for (const auto &mut : changed_mut) {
        while (iter != end && iter->get_position() < mut.get_position()) {
            if (iter->get_all_major_allele() != iter->get_par_one_hot()) {
                assert(changed||no_check);
                checked = true;
                return;
            }
            if (iter->get_par_one_hot() != iter->get_mut_one_hot()) {
                assert(state_change_iter->position == iter->get_position());
                assert(state_change_iter->new_state == iter->get_par_one_hot());
                assert(state_change_iter->old_state == iter->get_mut_one_hot());
                state_change_iter++;
            }
            iter++;
        }
        if (iter != end && iter->get_position() == mut.get_position()) {
            if (iter->get_all_major_allele() != mut.get_all_major_allele()) {
                assert(changed||no_check);
                checked = true;
                return;
            }
            if (iter->get_mut_one_hot() != mut.get_mut_one_hot()) {
                assert(state_change_iter->position == iter->get_position());
                assert(state_change_iter->new_state == mut.get_mut_one_hot());
                assert(state_change_iter->old_state == iter->get_mut_one_hot());
                state_change_iter++;
            }
            iter++;
        } else {
            if (mut.get_all_major_allele() != mut.get_par_one_hot()) {
                assert(changed||no_check);
                checked = true;
                return;
            }
            if (mut.get_mut_one_hot() != mut.get_par_one_hot()) {
                assert(state_change_iter->position == mut.get_position());
                assert(state_change_iter->new_state == mut.get_mut_one_hot());
                assert(state_change_iter->old_state == mut.get_par_one_hot());
                state_change_iter++;
            }
        }
    }
    while (iter != end) {
        if (iter->get_all_major_allele() != iter->get_par_one_hot()) {
            assert(changed||no_check);
            checked = true;
            return;
        }
        if (iter->get_par_one_hot() != iter->get_mut_one_hot()) {
            assert(state_change_iter->position == iter->get_position());
            assert(state_change_iter->new_state == iter->get_par_one_hot());
            assert(state_change_iter->old_state == iter->get_mut_one_hot());
            state_change_iter++;
        }
        iter++;
    }
    assert(!changed || checked||no_check);
    assert(state_change_iter == state_changed.end());
}
#endif
void check_major_state(MAT::Node *node, const MAT::Tree &new_tree) {
    auto ref_node = new_tree.get_node(node->identifier);
    auto ref_mut_iter = ref_node->mutations.begin();
    auto ref_mut_end = ref_node->mutations.end();
    for (const auto &mut : node->mutations) {
        while (ref_mut_iter != ref_mut_end &&
               ref_mut_iter->get_position() < mut.get_position()) {
            assert(ref_mut_iter->get_mut_one_hot() ==
                       ref_mut_iter->get_all_major_allele() &&
                   (!ref_mut_iter->get_boundary1_one_hot()||node->children.size()<=1));
            assert(
                ref_mut_iter->get_mut_one_hot() ==
                get_parent_state(node, ref_mut_iter->get_position()));
            ref_mut_iter++;
        }
        if (ref_mut_iter != ref_mut_end &&
            ref_mut_iter->get_position() == mut.get_position()) {
            assert(ref_mut_iter->get_all_major_allele() ==
                   mut.get_all_major_allele());
            assert(ref_mut_iter->get_boundary1_one_hot() ==
                   mut.get_boundary1_one_hot());
            if (node->children.size() == 2) {
                assert(ref_mut_iter->get_left_child_state() ==
                       mut.get_left_child_state());
                assert(ref_mut_iter->get_right_child_state() ==
                       mut.get_right_child_state());
            }
            ref_mut_iter++;
        } else {
            assert(
                mut.get_mut_one_hot() == mut.get_all_major_allele() &&
                (!mut.get_boundary1_one_hot() || node->children.size() <= 1));
            assert(mut.get_mut_one_hot() ==
                   get_parent_state(ref_node, mut.get_position()));
        }
    }
    while (ref_mut_iter != ref_mut_end) {
        assert(ref_mut_iter->get_mut_one_hot() ==
                   ref_mut_iter->get_all_major_allele() &&
               (!ref_mut_iter->get_boundary1_one_hot()));
        assert(ref_mut_iter->get_mut_one_hot() ==
               get_parent_state(node, ref_mut_iter->get_position()));
        ref_mut_iter++;
    }
}
static bool adjust_single_child_node(MAT::Node *node,MAT::Mutations_Collection& new_mut,State_Change_Collection& state_changes) {
    bool changed = false;
    const MAT::Mutations_Collection &child_mut = node->children[0]->mutations;
    auto child_mut_iter = child_mut.begin();
    auto child_mut_end = child_mut.end();
    for (const auto &ori_mut : node->mutations) {
        while (child_mut_iter != child_mut_end &&
               (child_mut_iter->get_position() < ori_mut.get_position())) {
            if (child_mut_iter->get_par_one_hot() !=
                child_mut_iter->get_all_major_allele()) {
                changed = true;
                new_mut.push_back(*child_mut_iter);
                new_mut.back().set_boundary_one_hot(
                    0xf & (~child_mut_iter->get_all_major_allele()));
                if (child_mut_iter->is_valid()) {
                    state_changes.emplace_back(new_mut.back(),child_mut_iter->get_par_one_hot());
                }
            }
            child_mut_iter++;
        }
        if ((child_mut_iter != child_mut_end) &&
            (child_mut_iter->get_par_one_hot() !=
             child_mut_iter->get_all_major_allele()) &&
            (child_mut_iter->get_position() == ori_mut.get_position())) {
            new_mut.push_back(ori_mut);
            changed = true;
            nuc_one_hot new_state=child_mut_iter->get_mut_one_hot();
            if (ori_mut.get_par_one_hot()&child_mut_iter->get_all_major_allele()) {
                new_state=ori_mut.get_par_one_hot()&child_mut_iter->get_all_major_allele();
            }
            new_mut.back().set_mut_one_hot(new_state);
            new_mut.back().set_auxillary(
                child_mut_iter->get_all_major_allele(),
                0xf & (~child_mut_iter->get_all_major_allele()));
                if (ori_mut.get_mut_one_hot()!=new_mut.back().get_mut_one_hot()) {
                    state_changes.emplace_back(new_mut.back(),ori_mut.get_mut_one_hot());
                }
            } else {
            if (ori_mut.is_valid()) {
                new_mut.push_back(ori_mut);
            }else {
                changed=true;
            }
        }
        if ((child_mut_iter != child_mut_end) &&
            (child_mut_iter->get_position() == ori_mut.get_position())) {
            child_mut_iter++;
        }
    }
    while (child_mut_iter != child_mut_end) {
        if (child_mut_iter->get_par_one_hot() !=
            child_mut_iter->get_all_major_allele()) {
            changed = true;
            new_mut.push_back(*child_mut_iter);
            new_mut.back().set_boundary_one_hot(
                0xf & (~child_mut_iter->get_all_major_allele()));
            if (child_mut_iter->is_valid()) {
                state_changes.emplace_back(new_mut.back(),child_mut_iter->get_par_one_hot());
            }
        }
        child_mut_iter++;
    }
    return changed;
}
static bool adjust_node(MAT::Node *node,
                        State_Change_Collection &changed_mutation
#ifdef CHECK_STATE_REASSIGN
                        ,
                        MAT::Tree &new_tree
#endif
) {
    MAT::Mutations_Collection new_major_alleles;
    bool changed;
    if (node->children.size() == 2) {
        changed =
            get_major_allele_binary(node, new_major_alleles, changed_mutation);
    } else if (node->children.size() == 1) {
        changed = adjust_single_child_node(node,new_major_alleles, changed_mutation);
    } else {
        changed = get_major_allele_polytomy(node, new_major_alleles,
                                            changed_mutation);
    }
#ifdef CHECK_STATE_REASSIGN
    check_changed(node->mutations, new_major_alleles, changed,
                  changed_mutation,node->is_root());
#endif
    node->mutations.swap(new_major_alleles);
#ifdef CHECK_STATE_REASSIGN
    check_major_state(node, new_tree);
#endif
    return changed;
}

void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
    ,
    MAT::Tree &new_tree
#endif
) {
    Backward_Pass_Heap heap(altered_nodes_in, nodes_with_changed_states_out);
#ifdef CHECK_STATE_REASSIGN
    size_t last_idx = 0xffffffff;
#endif
    bool changed;
    do {
#ifdef CHECK_STATE_REASSIGN
        assert(heap->altered_node->dfs_index < last_idx);
        last_idx = heap->altered_node->dfs_index;
#endif
        changed = adjust_node(heap->altered_node, heap->changed_states
#ifdef CHECK_STATE_REASSIGN
                              ,
                              new_tree
#endif
        );
        if (changed) {
            heap->altered_node->changed=true;
        }
    } while (heap.next(changed));
}