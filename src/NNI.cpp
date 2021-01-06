#include "src/mutation_annotated_tree.hpp"
#include "usher_graph.hpp"
#include <algorithm>
#include <utility>
#include <vector>
/*
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
static MAT::Mutation reverse_mutation(MAT::Mutation to_reverse) {
    auto old_par_nuc = to_reverse.par_nuc;
    to_reverse.par_nuc = to_reverse.mut_nuc;
    to_reverse.mut_nuc = old_par_nuc;
    return to_reverse;
}

/**
 * @brief Merge a set of new mutations into a set of old mutations
 *
 * @param older set of mutations happended earlier sorted in increasing order of
 * position
 * @param newer set of mutations happended later sorted in increasing order of
 * position
 * @return std::vector<MAT::Mutation> Merged set of mutations sorted in
 * increasing order of position
 */
static std::vector<MAT::Mutation>
merge_mutations(const std::vector<MAT::Mutation> &older,
                const std::vector<MAT::Mutation> &newer) {
    std::vector<MAT::Mutation> output;
    output.reserve(newer.size() + older.size());

    auto newer_iter = newer.begin();
    // basically merge sort
    for (auto old_mutations : older) {
        while (newer_iter != newer.end() &&
               old_mutations.position > newer_iter->position) {
            output.push_back(*newer_iter);
            newer_iter++;
        }
        if (newer_iter == newer.end() ||
            old_mutations.position < newer_iter->position) {
            output.push_back(old_mutations);
        } else {
            assert(newer_iter->position == old_mutations.position);
            assert(newer_iter->par_nuc == old_mutations.mut_nuc);

            if (newer_iter->mut_nuc != old_mutations.par_nuc) {
                MAT::Mutation temp(old_mutations);
                temp.mut_nuc = newer_iter->mut_nuc;
                output.push_back(temp);
            }
        }
    }
    return output;
}

static std::vector<MAT::Mutation>
get_mutations_after_moving_to_child(MAT::Node *to_move, MAT::Node *dest) {
    assert(dest->parent == to_move->parent);
    std::vector<MAT::Mutation> reversed_parent_mutations;
    reversed_parent_mutations.reserve(to_move->mutations.size());
    for (auto mutation : to_move->mutations) {
        reversed_parent_mutations.push_back(reverse_mutation(mutation));
    }
    return merge_mutations(reversed_parent_mutations, to_move->mutations);
}
/*
static unsigned int min_pos(Merged_Mutation_Iterator &in1,
                            Merged_Mutation_Iterator &in2) {
    if (!in1)
        return in2 ? in2->position : 0x7fffffff;
    if (!in2)
        return in1 ? in1->position : 0x7fffffff;
    return in1->position < in2->position ? in1->position : in2->position;
}

static std::pair<unsigned int,unsigned int> count_mutations(Merged_Mutation_Iterator &iter,
                                    int8_t nuc_to_count) {
    auto old_pos = iter->position;
    unsigned int count = 0;
    unsigned int other_count = 0;
    while (old_pos == iter->position) {
        if (iter->mut_nuc == nuc_to_count) {
            count++;
        } else {
            other_count++;
        }
        iter++;
    }
    return std::make_pair(count, other_count);
}

static int reassign_state(Merged_Mutation_Iterator &added,
                          Merged_Mutation_Iterator &removed,
                          Merged_Mutation_Iterator &rest,int num_chlidren,
                          std::vector<MAT::Mutation>& may_remove,
                          std::vector<MAT::Mutation>& may_add,
                          std::vector<MAT::Mutation>& add) {
    int parsimony_score = 0;
    while (rest) {
        if (rest->position < min_pos(added, removed))
            continue;
        while (added && rest->position > added->position) {
            added++;
            parsimony_score++;
        }

        while (removed && rest->position > removed->position) {
            removed++;
            parsimony_score--;
        }
        if (removed && removed->position == rest->position) {
            auto changes=count_mutations(rest, removed->mut_nuc);
            //parsimony score change if this mutation is removed from the node
            int parsimony_change=changes.first //parsimony score increase from pushing this mutation to children
            -(num_chlidren-changes.first-changes.second);
            //just fast forward removed to a new position
            count_mutations(removed, removed->mut_nuc);

        }
    }
    return parsimony_score;
}
*/
static int check_exchange_improvement_after_move_to_parent_of_parent(
    MAT::Node *to_check, const MAT::Node *exchanger,
    std::vector<MAT::Node *> &mutation_wrt_cur_parent_of_parent) {
    // parsimony score change accumulator
    int parsimony_change = 0;
    std::vector<MAT::Mutation> to_check_mutation_wrt_parent =
        get_mutations_after_moving_to_child(to_check, exchanger->parent);

    // Do sankoff on parent and see whether its state need to be changed
    // Only mutations in exchanger(originally a child of parent) may be removed,
    // and mutations in to_check_mutation_wrt_parent (a new child of parent) may
    // be added
    Merged_Mutation_Iterator other_children_mutations(exchanger->parent,
                                                      exchanger);

}
/*
        parent_of_parent
        /           \
    to_check        parent
                        \
                        this_node
*/
static std::pair<MAT::Node *, int>
check_children_of_parent_of_parent(MAT::Node *this_node) {
    unsigned int improvement = 0;
    MAT::Node *replacer = nullptr;
    std::vector<MAT::Mutation> mutations_wrt_parent_of_parent_before_move =
        merge_mutations(this_node->parent->mutations, this_node->mutations);
    for (MAT::Node *to_check : this_node->parent->children) {
        // Not particularly interesting to exchange this_node with itself
        if (to_check == this_node) {
            continue;
        }
        int this_improvement = check_exchange_improvement(
            to_check, this_node, mutations_wrt_parent_of_parent_before_move);
        if (this_improvement > improvement) {
            improvement = this_improvement;
            replacer = to_check;
        }
    }
    return std::make_pair(replacer, improvement);
}
bool serial_NNI(MAT::Node *this_node) {}