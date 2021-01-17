#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement.hpp"
#include "usher_graph.hpp"
#include <algorithm>
#include <cctype>
#include <string>
#include <utility>
#include <vector>
#include "check_samples.hpp"
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
typedef std::tuple<int, std::vector<Fitch_Sankoff::States_Type>,std::vector<Fitch_Sankoff::Scores_Type>>
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



static void insert_node(MAT::Node *parent, MAT::Node *to_insert,
                        MAT::Tree &tree) {
    Mutation_Annotated_Tree::Mutations_Collection best_common;
    Mutation_Annotated_Tree::Mutations_Collection best_sibling_unique;
    Mutation_Annotated_Tree::Mutations_Collection best_insert_unique;
    MAT::Node *best_sibling;
    for (auto child : parent->children) {
        Mutation_Annotated_Tree::Mutations_Collection this_common;
        Mutation_Annotated_Tree::Mutations_Collection this_sibling_unique;
        Mutation_Annotated_Tree::Mutations_Collection this_insert_unique;
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
            std::to_string(++tree.curr_internal_node), parent->identifier);
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
           MAT::Tree &tree,std::vector<Fitch_Sankoff::Scores_Type> scores_all_pos) {
    MAT::Mutations_Collection::iterator mutations_begin;
    MAT::Mutations_Collection::iterator mutations_end;
    set_mutation_iterator<move_to_parent>(this_node, new_parent,
                                          mutations_begin, mutations_end);
    Sample_Mut_Type ori;
    auto node_check=this_node->parent;
    check_samples(node_check,ori);
    // start moving
    MAT::Node* parent = this_node->parent;
    std::vector<MAT::Node*>& parent_children = parent->children;
    auto iter =
        std::find(parent_children.begin(), parent_children.end(), this_node);
    parent_children.erase(iter);
    if (parent->children.size() <= 1 && parent->not_sample()) {
        tree.remove_node(parent->identifier, false);
    }
    // Try merging with sibling in the new location
    insert_node(new_parent,this_node,tree);
    this_node->parent=new_parent;
    // apply adjustment to subtree
    auto score_iter=scores_all_pos.begin();
    assert(states_all_pos.size()==mutations_end-mutations_begin);
    for (auto states : states_all_pos) {
        assert(states.back().node==this_node);
        Fitch_Sankoff::sankoff_forward_pass(range, states, dfs_ordered_nodes,
                                            *mutations_begin,move_to_parent?mutations_begin->par_nuc:mutations_begin->mut_nuc,*score_iter);
        mutations_begin++;
        score_iter++;
    }
    assert(mutations_begin == mutations_end);
    check_samples(node_check,ori);
}

template <bool move_to_parent>
static pending_change_type check_move_profitable(
    MAT::Node *this_node, const std::pair<size_t, size_t> &range,
    MAT::Node *new_parent, const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<Fitch_Sankoff::States_Type> states_all_pos;
    std::vector<Fitch_Sankoff::Scores_Type> scores_all_pos;

    MAT::Mutations_Collection::iterator mutations_begin;
    MAT::Mutations_Collection::iterator mutations_end;
    set_mutation_iterator<move_to_parent>(this_node, new_parent,
                                          mutations_begin, mutations_end);

    states_all_pos.reserve(mutations_end - mutations_begin);
    int score_change = 0; // change in parsimony score if moved to new location
    for (auto iter=mutations_begin; iter < mutations_end; iter++) {
        Fitch_Sankoff::Scores_Type scores;
        scores.reserve(range.second - range.first);
        Fitch_Sankoff::States_Type states;
        states.reserve(range.second - range.first);
        Fitch_Sankoff::sankoff_backward_pass(range, *iter,
                                             dfs_ordered_nodes, scores, states);
        char cur_parent_nuc = 31-__builtin_clz((unsigned int)(move_to_parent ? iter->mut_nuc
                                             : iter->par_nuc));
        char new_parent_nuc = 31-__builtin_clz((unsigned int)(move_to_parent ? iter->par_nuc
                                             : iter->mut_nuc));
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
        scores_all_pos.push_back(scores);
    }
    assert(states_all_pos.size()==mutations_end-mutations_begin);
    if (score_change >= 0) {
        states_all_pos.clear();
        scores_all_pos.clear();
    }
    return std::make_tuple(score_change, states_all_pos,scores_all_pos);
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
    auto best_score = std::get<0>(pending_change); 
    std::vector<Fitch_Sankoff::States_Type> best_states = std::get<1>(pending_change);
    std::vector<Fitch_Sankoff::Scores_Type> best_scores = std::get<2>(pending_change);
    for (auto sibling : this_node->parent->children) {
        if (sibling == this_node)
            continue;
        pending_change = check_move_profitable<false>(this_node, range, sibling,
                                                      dfs_ordered_nodes);
        if (std::get<0>(pending_change) < best_score) {
            best_score = std::get<0>(pending_change);
            best_new_parent = sibling;
            best_states.swap(std::get<1>(pending_change));
            best_scores.swap(std::get<2>(pending_change));
        }
    }
    if (best_score >= 0) {
        return false;
    }
    if (best_new_parent == this_node->parent->parent) {
        apply_move<true>(this_node, range, best_new_parent, dfs_ordered_nodes,
                         best_states, tree,best_scores);
    } else {
        apply_move<false>(this_node, range, best_new_parent, dfs_ordered_nodes,
                          best_states, tree,best_scores);
    }
    return true;
}