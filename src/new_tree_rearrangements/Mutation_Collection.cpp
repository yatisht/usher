#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>
using namespace Mutation_Annotated_Tree;
void Mutations_Collection::merge_out(const Mutations_Collection &other,
                                     Mutations_Collection &out,
                                     char keep_self) const {
// used for checking whether the two vectors are sorted while merging
#ifdef DETAIL_DEBUG_MUTATION_SORTED
#define mutation_vector_check_order(newly_inserted)                            \
    assert(last_pos_inserted < (newly_inserted));                              \
    last_pos_inserted = (newly_inserted);
    int last_pos_inserted = -1;
#else
#define mutation_vector_check_order(newly_inserted)
#endif
    out.mutations.reserve(other.mutations.size() + mutations.size());
    auto other_iter = other.mutations.begin();
    for (auto this_mutation : mutations) {
        while (other_iter != other.mutations.end() &&
               other_iter->get_position() < this_mutation.get_position()) {
            mutation_vector_check_order(other_iter->get_position());
            out.mutations.push_back(*other_iter);
            if (keep_self == INVERT_MERGE) {
                auto temp = out.mutations.back().get_mut_one_hot();
                out.mutations.back().set_mut_one_hot(out.mutations.back().get_par_one_hot());
                out.mutations.back().set_par_one_hot(temp);
            }
            other_iter++;
        }
        if (other_iter == other.mutations.end() ||
            this_mutation.get_position() < other_iter->get_position()) {
            mutation_vector_check_order(this_mutation.get_position());
            out.mutations.push_back(this_mutation);
        } else {
            mutation_vector_check_order(this_mutation.get_position());

            assert(this_mutation.get_position() == other_iter->get_position());
            switch (keep_self) {
            case NO_DUPLICATE:
                assert(false);
            case KEEP_OTHER:
                out.mutations.push_back(*other_iter);
                break;
            case KEEP_SELF:
                out.mutations.push_back(this_mutation);
                break;
            case MERGE:
                assert(other_iter->get_par_one_hot() == this_mutation.get_mut_one_hot());
                if (other_iter->get_mut_one_hot() != this_mutation.get_ref_one_hot()) {
                    out.mutations.push_back(*other_iter);
                    out.mutations.back().set_par_one_hot(this_mutation.get_par_one_hot());
                }
                break;
            default:
                assert(keep_self == INVERT_MERGE);
                assert(other_iter->get_mut_one_hot() == this_mutation.get_mut_one_hot());
                if (other_iter->get_par_one_hot() != this_mutation.get_ref_one_hot()) {
                    out.mutations.push_back(this_mutation);
                    out.mutations.back().get_mut_one_hot() = other_iter->get_par_one_hot();
                }
                break;
            }
            other_iter++;
        }
    }
    while (other_iter < other.mutations.end()) {
        mutation_vector_check_order(other_iter->get_position());
        out.mutations.push_back(*other_iter);
        other_iter++;
    }
}

void Mutations_Collection::set_difference(const Mutations_Collection &other,
                                          Mutations_Collection &this_unique,
                                          Mutations_Collection &other_unique,
                                          Mutations_Collection &common) const {
    this_unique.mutations.reserve(mutations.size());
    other_unique.mutations.reserve(other.mutations.size());
    common.mutations.reserve(
        std::min(mutations.size(), other.mutations.size()));
    auto other_iter = other.mutations.begin();
#ifdef DETAIL_DEBUG_MUTATION_SORTED
    int last_pos_inserted = -1;
#endif
    // merge sort again
    for (auto this_mutation : mutations) {
        while (other_iter != other.mutations.end() &&
               other_iter->get_position() < this_mutation.get_position()) {
            mutation_vector_check_order(other_iter->get_position());
            other_unique.mutations.push_back(*other_iter);
            other_iter++;
        }
        if (other_iter == other.mutations.end() ||
            this_mutation.get_position() < other_iter->get_position()) {
            mutation_vector_check_order(this_mutation.get_position());
            this_unique.mutations.push_back(this_mutation);
        } else {
            mutation_vector_check_order(this_mutation.get_position());
            assert(this_mutation.get_position() == other_iter->get_position());

            if (other_iter->get_mut_one_hot() == this_mutation.get_mut_one_hot()) {
                assert(other_iter->get_par_one_hot() == this_mutation.get_par_one_hot());
                common.mutations.push_back(this_mutation);
            } else {
                this_unique.mutations.push_back(this_mutation);
                other_unique.mutations.push_back(*other_iter);
            }
            other_iter++;
        }
    }
    while (other_iter < other.mutations.end()) {
        mutation_vector_check_order(other_iter->get_position());
        other_unique.mutations.push_back(*other_iter);
        other_iter++;
    }
    assert(this_unique.size() + other_unique.size() + 2 * common.size() ==
           mutations.size() + other.mutations.size());
}
void Mutations_Collection::batch_find(Mutations_Collection &target) {
    // yet another merge sort
#ifdef DETAIL_DEBUG_MUTATION_SORTED
    int last_pos_inserted = -1;
#endif
    auto target_iter = target.mutations.begin();
    for (auto this_mutation : mutations) {
        if (target_iter == target.mutations.end()) {
            return;
        }
        while (target_iter->get_position() < this_mutation.get_position()) {
            mutation_vector_check_order(target_iter->get_position());
            target_iter++;
            if (target_iter == target.mutations.end()) {
                return;
            }
        }
        if (this_mutation.get_position() < target_iter->get_position()) {
            mutation_vector_check_order(this_mutation.get_position());
        } else {
            mutation_vector_check_order(this_mutation.get_position());
            assert(this_mutation.get_position() == target_iter->get_position());
            *target_iter = this_mutation;
            target_iter++;
        }
    }
    if (target_iter != target.mutations.end()) {
        mutation_vector_check_order(target_iter->get_position());
    }
}

Mutations_Collection::iterator Mutations_Collection::find_next(int pos) {
#ifdef DETAIL_DEBUG_MUTATION_SORTED
    int lastPos = 0;
#endif
    auto iter = mutations.begin();
    for (; iter < mutations.end(); iter++) {
// check sorting in debug mode
#ifdef DETAIL_DEBUG_MUTATION_SORTED
        assert(lastPos < iter->get_position());
        lastPos = iter->get_position();
#endif
        if (iter->get_position() >= pos) {
            break;
        }
    }
    assert(iter == mutations.begin() || (iter - 1)->get_position() < pos);
    assert(iter == mutations.end() || (iter)->get_position() >= pos);
    return iter;
}


void Mutation_Annotated_Tree::Mutations_Collection::remove_invalid(){
    std::vector<Mutation_Annotated_Tree::Mutation> out;
    out.reserve(mutations.size());
    for (const auto& mut : mutations) {
        if (mut.is_valid()) {
            out.push_back(mut);
        }
    }
    mutations.swap(out);
}