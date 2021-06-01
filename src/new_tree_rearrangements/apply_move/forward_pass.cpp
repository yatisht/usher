#include "apply_move.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstdio>
#include <vector>
State_Change_Collection merge(const State_Change_Collection &in1,
                              const State_Change_Collection &in2) {
    auto in1_iter = in1.begin();
    auto in1_end = in1.end();
    bool in1_new = false;
    State_Change_Collection out;
    if (in1_iter->position == NEW_MARK) {
        in1_new = true;
        in1_iter++;
        assert(in2.begin()->position != NEW_MARK);
    } else {
        assert(in2.begin()->position == NEW_MARK);
    }
    for (const auto &in2_change : in2) {
        if (in2_change.position==NEW_MARK) {
            continue;
        }
        while (in1_iter != in1_end &&
               in1_iter->position < in2_change.position) {
            out.push_back(*in1_iter);
            in1_iter++;
        }
        out.push_back(in2_change);
        if (in1_iter != in1_end && in1_iter->position == in2_change.position) {
            if (in1_new) {
                out.back().new_state = in1_iter->new_state;
            } else {
                out.back().old_state = in1_iter->old_state;
            }
            if (out.back().new_state == out.back().old_state) {
                out.pop_back();
            }
            in1_iter++;
        }
    }
    while (in1_iter != in1_end) {
        out.push_back(*in1_iter);
        in1_iter++;
    }
    return out;
}

struct Forward_Pass_Comparator {
    bool operator()(const Altered_Node_t &first, const Altered_Node_t &second) {
        return first.altered_node->dfs_index > second.altered_node->dfs_index;
    }
};

class Forward_Pass_Heap {
    std::vector<Altered_Node_t> altered_nodes;
    void increment() {
        if (altered_nodes.size() == 1) {
            return;
        }
        std::pop_heap(altered_nodes.begin(), altered_nodes.end(),
                      Forward_Pass_Comparator());
        while (altered_nodes.front().altered_node ==
               altered_nodes.back().altered_node) {
            std::pop_heap(altered_nodes.begin(), altered_nodes.end() - 1,
                          Forward_Pass_Comparator());
            size_t last_idx = altered_nodes.size() - 1;
            assert(altered_nodes[last_idx].altered_node ==
                   altered_nodes[last_idx - 1].altered_node);
            altered_nodes[last_idx - 1].changed_states =
                merge(altered_nodes[last_idx - 1].changed_states,
                      altered_nodes[last_idx].changed_states);
            altered_nodes.pop_back();
            if (altered_nodes.back().changed_states.empty()) {
                altered_nodes.pop_back();
                std::pop_heap(altered_nodes.begin(), altered_nodes.end(),
                    Forward_Pass_Comparator());
            }
            if (altered_nodes.size() <= 1) {
                return;
            }
        }
    }

  public:
    Forward_Pass_Heap(std::vector<Altered_Node_t> &in) {
        altered_nodes = std::move(in);
        std::make_heap(altered_nodes.begin(), altered_nodes.end(),
                       Forward_Pass_Comparator());
        increment();
    }

    Altered_Node_t &operator*() { return altered_nodes.back(); }
    void pop_back() {
         altered_nodes.pop_back(); }

    void push_back(Altered_Node_t &altered) {
        altered_nodes.push_back(altered);
        std::push_heap(altered_nodes.begin(), altered_nodes.end(),
                       Forward_Pass_Comparator());
    }

    bool next() {
        if (altered_nodes.empty()) {
            return false;
        }
        increment();
        return !altered_nodes.empty();
    }
};

void unmatched_parent_state_change(
    MAT::Node *&node,
    MAT::Mutations_Collection &new_mut,
    const state_change &change
#ifdef CHECK_STATE_REASSIGN
    ,MAT::Mutations_Collection::const_iterator
        &ref_iter,
    MAT::Mutations_Collection::const_iterator
        &ref_end
#endif
        ) {
    MAT::Mutation mut(change.chr_idx, change.position, change.new_state,
                      change.old_state);
    if (node->children.size() <= 1) {
        mut.set_boundary_one_hot(0xf & (~mut.get_mut_one_hot()));
        if (mut.get_all_major_allele() !=
            mut.get_par_one_hot()) {
            new_mut.push_back(mut);
#ifdef CHECK_STATE_REASSIGN
            assert(ref_iter != ref_end);
            assert(*ref_iter == new_mut.back());
            ref_iter++;
#endif
        }
    } else if (mut.get_all_major_allele() != mut.get_par_one_hot() ||
               (mut.get_boundary1_one_hot())) {
        new_mut.push_back(mut);
#ifdef CHECK_STATE_REASSIGN
        assert(ref_iter != ref_end);
        assert(*ref_iter == new_mut.back());
        ref_iter++;
#endif
    }
}
void set_state_from_parent(MAT::Node *node,
                           const State_Change_Collection &parent_altered,
                           State_Change_Collection &this_state_altered_out
#ifdef CHECK_STATE_REASSIGN
                           ,
                           MAT::Tree &new_tree
#endif
) {
    node->changed=true;
    MAT::Mutations_Collection new_mut;
    auto iter = parent_altered.begin();
    auto end = parent_altered.end();
    if (iter->position == NEW_MARK) {
        iter++;
    }
#ifdef CHECK_STATE_REASSIGN
    auto ref_node = new_tree.get_node(node->identifier);
    MAT::Mutations_Collection::const_iterator ref_iter = ref_node->mutations.begin();
    MAT::Mutations_Collection::const_iterator ref_end = ref_node->mutations.end();
#endif
    for (auto &node_mut : node->mutations) {
        while (iter != end && iter->position < node_mut.get_position()) {
            unmatched_parent_state_change(node, new_mut, *iter
#ifdef CHECK_STATE_REASSIGN
            , ref_iter, ref_end
#endif
            );
            iter++;
        }
        if (iter != end && iter->position == node_mut.get_position()) {
            nuc_one_hot par_state = iter->new_state;
            nuc_one_hot new_state = node_mut.get_all_major_allele() & par_state;
            nuc_one_hot old_state = node_mut.get_mut_one_hot();
            if (!new_state) {
                new_state = old_state;
            }
            node_mut.set_par_mut(par_state, new_state);
            if (par_state != node_mut.get_all_major_allele() ||
                (node_mut.get_boundary1_one_hot() &&
                 node->children.size() > 1)) {
                new_mut.push_back(node_mut);
#ifdef CHECK_STATE_REASSIGN
                assert(ref_iter != ref_end);
                assert(*ref_iter == new_mut.back());
                ref_iter++;
#endif
            }
            if (new_state != old_state) {
                this_state_altered_out.emplace_back(new_mut.back(), old_state);
            }
            iter++;
        } else {
            new_mut.push_back(node_mut);
#ifdef CHECK_STATE_REASSIGN
            assert(ref_iter != ref_end);
            assert(*ref_iter == new_mut.back());
            ref_iter++;
#endif
        }
    }
    while (iter != end) {
        unmatched_parent_state_change(node, new_mut, *iter
#ifdef CHECK_STATE_REASSIGN
        , ref_iter, ref_end
#endif
);
        iter++;
    }
#ifdef CHECK_STATE_REASSIGN
    assert(ref_iter == ref_end);
#endif
    node->mutations.swap(new_mut);
}

void forward_pass(std::vector<Altered_Node_t> &in
#ifdef CHECK_STATE_REASSIGN
                  ,
                  MAT::Tree &new_tree
#endif
) {
#ifdef CHECK_STATE_REASSIGN
    int last_dfs_idx=0;
#endif
    Forward_Pass_Heap heap(in);
    do {
        Altered_Node_t altered = *heap;
#ifdef CHECK_STATE_REASSIGN
        assert(last_dfs_idx<altered.altered_node->dfs_index);
        last_dfs_idx=altered.altered_node->dfs_index;
#endif
        heap.pop_back();
        for (auto child : altered.altered_node->children) {
            Altered_Node_t out(child);
            out.changed_states.emplace_back();
            set_state_from_parent(child, altered.changed_states,
                                  out.changed_states
#ifdef CHECK_STATE_REASSIGN
                                  ,
                                  new_tree
#endif
            );
            if (out.changed_states.size() > 1) {
                heap.push_back(out);
            }
        }
    } while (heap.next());
}