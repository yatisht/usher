#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement.hpp"
#include "usher_graph.hpp"
#include <algorithm>
#include <cctype>
#include <string>
#include <utility>
#include <vector>
/* This is for non-adjacent moves
class Merged_Mutation_Iterator {
    class range {
        std::vector<MAT::Mutation>::const_iterator start;
        std::vector<MAT::Mutation>::const_iterator end;

      public:
        range(std::vector<MAT::Mutation> in)
            : start(in.begin()), end(in.end()) {}
        const std::vector<MAT::Mutation>::const_iterator &operator->() const {
            return start;
        }
        const std::vector<MAT::Mutation>::const_iterator &operator*() const {
            return start;
        }
        void operator++(int) { start++; }
        operator bool() { return start != end; }
    };
    struct comparator {
        bool operator()(const range &lhs, const range rhs) {
            return lhs->position > rhs->position;
        }
    };
    std::vector<range> heap;
    static std::vector<range>
    create_from_node_helper(const MAT::Node *parent, const MAT::Node *exclude) {
        std::vector<range> iterators;

        return iterators;
    }

  public:
    Merged_Mutation_Iterator(const MAT::Node *parent,
                             const MAT::Node *exclude) {
        heap.reserve(parent->children.size());
        for (MAT::Node *child : parent->children) {
            if (child == exclude)
                continue;
            if (child->mutations.size()) {
                heap.emplace_back(child->mutations);
            }
        }
        std::make_heap(heap.begin(), heap.end(), comparator());
    }
    void operator++(int) {
        std::pop_heap(heap.begin(), heap.end(), comparator());
        range &incremented_iter = heap.back();
        if (incremented_iter) {
            std::push_heap(heap.begin(), heap.end(), comparator());
        } else {
            heap.pop_back();
        }
    }
    const std::vector<MAT::Mutation>::const_iterator &operator->() const {
        return *(heap.front());
    }
    operator bool() { return !heap.empty(); }
};
*/
typedef std::pair<int, std::vector<Fitch_Sankoff::States_Type>>
    pending_change_type;
static MAT::Mutation reverse_mutation(MAT::Mutation to_reverse) {
    auto old_par_nuc = to_reverse.par_nuc;
    to_reverse.par_nuc = to_reverse.mut_nuc;
    to_reverse.mut_nuc = old_par_nuc;
    return to_reverse;
}
template <bool move_to_parent>
static void
set_mutation_iterator(MAT::Node *this_node, MAT::Node *new_parent,
                      MAT::Mutations_Collection::iterator &mutations_begin,
                      MAT::Mutations_Collection::iterator &mutations_end) {
    if (move_to_parent) {
        assert(new_parent == this_node->parent->parent);
        auto parent_of_parent = this_node->parent->parent;
        mutations_begin = parent_of_parent->mutations.begin();
        mutations_end = parent_of_parent->mutations.end();
    } else {
        assert(new_parent->parent == this_node->parent);
        mutations_begin = new_parent->mutations.begin();
        mutations_end = new_parent->mutations.end();
    }
}

// HACK: try to identify whether a node is a sample from name
static bool not_sample(MAT::Node *node) {
    for (auto c : node->identifier) {
        if (!std::isdigit(c))
            return false;
    }
    return true;
}

static void insert_node(MAT::Node *parent, MAT::Node *to_insert,
                        MAT::Tree &tree) {
    std::vector<MAT::Mutation> best_common;
    std::vector<MAT::Mutation> best_sibling_unique;
    std::vector<MAT::Mutation> best_insert_unique;
    MAT::Node *best_sibling;
    for (auto child : parent->children) {
        std::vector<MAT::Mutation> this_common;
        std::vector<MAT::Mutation> this_sibling_unique;
        std::vector<MAT::Mutation> this_insert_unique;
        child->mutations.set_difference(to_insert->mutations,
                                        this_sibling_unique, this_insert_unique,
                                        this_common);
        if (this_common.size() > best_common.size()) {
            best_sibling = child;
            best_common.swap(this_common);
            best_sibling_unique.swap(this_sibling_unique);
            best_insert_unique.swap(this_insert_unique);
        }
    }

    auto &parent_children = parent->children;
    if (best_common.size() == 0) {
        parent_children.push_back(to_insert);
    } else {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              best_sibling);
        parent_children.erase(iter);

        auto new_node = tree.create_node(
            std::to_string(tree.curr_internal_node++), parent->identifier);
        new_node->mutations.swap(best_common);

        to_insert->mutations.swap(best_insert_unique);
        new_node->children.push_back(to_insert);

        best_sibling->mutations.swap(best_sibling_unique);
        new_node->children.push_back(best_sibling);

        parent_children.push_back(new_node);
    }
}

template <bool move_to_parent>
static void
apply_move(MAT::Node *this_node, const std::pair<size_t, size_t> &range,
           MAT::Node *new_parent, std::vector<MAT::Node *> &dfs_ordered_nodes,
           std::vector<Fitch_Sankoff::States_Type> states_all_pos,
           MAT::Tree &tree) {
    MAT::Mutations_Collection::iterator mutations_begin;
    MAT::Mutations_Collection::iterator mutations_end;
    set_mutation_iterator<move_to_parent>(this_node, new_parent,
                                          mutations_begin, mutations_end);
    // apply adjustment
    for (auto states : states_all_pos) {
        Fitch_Sankoff::sankoff_forward_pass(range, states, dfs_ordered_nodes,
                                            *mutations_begin);
        mutations_begin++;
    }
    assert(mutations_begin == mutations_end);
    // start moving
    auto parent = this_node->parent;
    auto &parent_children = parent->children;
    auto iter =
        std::find(parent_children.begin(), parent_children.end(), this_node);
    parent_children.erase(iter);
    if (parent->children.size() <= 1 && not_sample(parent)) {
        tree.remove_node(parent->identifier, false);
    }
    // Try merging with sibling in the new location
}

template <bool move_to_parent>
static pending_change_type check_move_profitable(
    MAT::Node *this_node, const std::pair<size_t, size_t> &range,
    MAT::Node *new_parent, const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<Fitch_Sankoff::States_Type> states_all_pos;

    MAT::Mutations_Collection::iterator mutations_begin;
    MAT::Mutations_Collection::iterator mutations_end;
    set_mutation_iterator<move_to_parent>(this_node, new_parent,
                                          mutations_begin, mutations_end);

    states_all_pos.reserve(mutations_end - mutations_begin);
    int score_change = 0; // change in parsimony score if moved to new location
    for (; mutations_begin < mutations_end; mutations_begin++) {
        Fitch_Sankoff::Scores_Type scores;
        scores.reserve(range.second - range.first);
        Fitch_Sankoff::States_Type states;
        states.reserve(range.second - range.first);
        Fitch_Sankoff::sankoff_backward_pass(range, *mutations_begin,
                                             dfs_ordered_nodes, scores, states);
        char cur_parent_nuc = move_to_parent ? mutations_begin->mut_nuc
                                             : mutations_begin->par_nuc;
        char new_parent_nuc = move_to_parent ? mutations_begin->par_nuc
                                             : mutations_begin->mut_nuc;
        // score for current parent
        auto curr_result = Fitch_Sankoff::get_child_score_on_par_nuc(
            cur_parent_nuc, scores.back());
        auto new_result = Fitch_Sankoff::get_child_score_on_par_nuc(
            new_parent_nuc, scores.back());
        // always prepare states for new parent, as nothing will be done if
        // there is no improvement and this_node stay with its original parent.
        states.back() = new_result.second;
        score_change += (new_result.first - curr_result.first);
        states_all_pos.push_back(states);
    }
    if (score_change >= 0) {
        states_all_pos.clear();
    }
    return std::make_pair(score_change, states_all_pos);
}

/**
 * @brief See whether there will be any improvement by moving current node as a
 * child of parent of parent of a child of its sibling. NNI can be achieved with
 * such moves:
 *       parent of parent
 *          /       \
 *         A        parent
 *                   /  \
 *                  B   this_node
 * Move this_node to parent of parent (parent and B fused):
 *      parent of parent
 *          /    |      \
 *         A     B      this_node
 * Make a new parent for A and B if they share mutations:
 *      parent of parent
 *         /      \
 *        *      this_node
 *       / \
 *      B   A
 * Finish a NNI of A and this_node
 * Moving under siblings unlikely be NNI
 *          parent
 *      /      |    \
 * this_node   B     sibling
 *                      \
 *                        A
 * Move this under sibling
 *    parent
 *  /       \
 *  B     sibling
 *        /    \
 *       A   this_node
 * @param this_node
 * @return bool Whether there is an improvement from moving this_node
 */
bool Tree_Rearrangement::move_nearest(
    MAT::Node *this_node, std::vector<MAT::Node *> &dfs_ordered_nodes,
    MAT::Tree &tree) {
    if (!this_node->parent || !this_node->parent->parent) {
        return false;
    }
    auto range = Fitch_Sankoff::dfs_range(this_node);
    // moving to parent
    auto best_new_parent = this_node->parent->parent;
    pending_change_type pending_change = check_move_profitable<true>(
        this_node, range, best_new_parent, dfs_ordered_nodes);
    auto best_score = pending_change.first;
    std::vector<Fitch_Sankoff::States_Type> best_states = pending_change.second;
    for (auto child : this_node->parent->children) {
        if (child == this_node)
            continue;
        pending_change = check_move_profitable<false>(this_node, range, child,
                                                      dfs_ordered_nodes);
        if (pending_change.first <= best_score) {
            best_score = pending_change.first;
            best_new_parent = this_node;
            best_states.swap(pending_change.second);
        }
    }
    if (best_score >= 0) {
        return false;
    }
    if (best_new_parent == this_node->parent->parent) {
        apply_move<true>(this_node, range, best_new_parent, dfs_ordered_nodes,
                         best_states, tree);
    } else {
        apply_move<false>(this_node, range, best_new_parent, dfs_ordered_nodes,
                          best_states, tree);
    }
    return true;
}