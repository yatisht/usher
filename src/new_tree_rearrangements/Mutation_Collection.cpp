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
void Mutations_Collection::finalize(){
    std::lock_guard<mutex_type> lock(mutex);
    if(!is_dirty()) return;
    std::vector<Mutation> new_content;
    std::sort(new_inserts.begin(),new_inserts.end());
    new_content.reserve(mutations.size()+new_inserts.size());
    auto new_inserts_iter = new_inserts.begin();
    auto remove_or_replace_iter=remove_or_replace.begin();
#ifdef DETAIL_DEBUG_MUTATION_SORTED
    int last_pos_inserted = -1;
#endif
    // merge sort again
    for (auto& this_mutation : mutations) {
        while (new_inserts_iter != new_inserts.end() &&
               new_inserts_iter->get_position() < this_mutation.get_position()) {
            mutation_vector_check_order(new_inserts_iter->get_position());
            new_content.push_back(*new_inserts_iter);
            new_inserts_iter++;
        }
        if (new_inserts_iter == new_inserts.end() ||
            this_mutation.get_position() < new_inserts_iter->get_position()) {
            mutation_vector_check_order(this_mutation.get_position());
            switch (*remove_or_replace_iter) {
                case -1: break;
                case 0: new_content.push_back(this_mutation);
                default:
                    new_content.back().set_mut_one_hot((*remove_or_replace_iter)&0xff);
                    new_content.back().set_par_one_hot((*remove_or_replace_iter)>>4);
            }
        } else {
            assert(false);
        }
        new_inserts_iter++;
    }
    while (new_inserts_iter < new_inserts.end()) {
        mutation_vector_check_order(new_inserts_iter->get_position());
        new_content.push_back(*new_inserts_iter);
        new_inserts_iter++;
    }
    dirty_flag=0;
    new_inserts.clear();
    remove_or_replace.clear();
    mutations.swap(new_content);
}
/*
bool Mutations_Collection::dirty_set_difference(Mutations_Collection& common, Mutations_Collection& original){
    std::vector<Mutation> changes;
    common.reserve(mutations.size());
    original.reserve(mutations.size());
    changes.reserve(mutations.size()+new_inserts.size());
    auto new_inserts_iter = new_inserts.begin();
    auto remove_or_replace_iter=remove_or_replace.begin();
#ifdef DETAIL_DEBUG_MUTATION_SORTED
    int last_pos_inserted = -1;
#endif
    // merge sort again
    for (auto& this_mutation : mutations) {
        while (new_inserts_iter != new_inserts.end() &&
               new_inserts_iter->get_position() < this_mutation.get_position()) {
            mutation_vector_check_order(new_inserts_iter->get_position());
            changes.push_back(*new_inserts_iter);
            new_inserts_iter++;
        }
        if (new_inserts_iter == new_inserts.end() ||
            this_mutation.get_position() < new_inserts_iter->get_position()) {
            mutation_vector_check_order(this_mutation.get_position());
            switch (*remove_or_replace_iter) {
                case -1: original.mutations.push_back(this_mutation); break;
                case 0: common.mutations.push_back(this_mutation); break;
                default:
                    original.mutations.push_back(this_mutation);
                    changes.push_back(this_mutation);
                    changes.back().mut_nuc=(*remove_or_replace_iter)&0xff;
                    changes.back().par_nuc=(*remove_or_replace_iter)>>4;
            }
        } else {
            assert(false);
        }
        new_inserts_iter++;
    }
    while (new_inserts_iter < new_inserts.end()) {
        mutation_vector_check_order(new_inserts_iter->get_position());
        changes.push_back(*new_inserts_iter);
        new_inserts_iter++;
    }
    dirty_flag=0;
    new_inserts.clear();
    remove_or_replace.clear();
    mutations.swap(changes);
    return true;
}
*/
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
std::pair<bool,bool> Mutations_Collection::dirty_remove(int pos) {

    auto iter = find(pos);
    std::pair<bool,bool> result={false,false};
    if (iter == mutations.end()) {
        return result;
    }
    result.second=is_dirty();
    result.first=true;
    init_dirty();
    remove_or_replace[iter - mutations.begin()] = 0xff;
    return result;
}
std::pair<bool,bool> Mutations_Collection::dirty_insert(const Mutation &mut, char keep_self) {
    assert(mut.get_par_one_hot() != mut.get_mut_one_hot());
    std::pair<bool,bool> result={false,false};
    auto iter = find_next(mut.get_position());
    if (iter != mutations.end() && iter->get_position() == mut.get_position()) {
        assert(keep_self != -1);
        if (!keep_self) {
            if (*iter == mut) {
                return result;
            }
            result.second=is_dirty();
            result.first=true;
            init_dirty();
            remove_or_replace[iter - mutations.begin()] =
                ((mut.get_par_one_hot() << 4) | mut.get_mut_one_hot());
        }
        return result;
    }
    result.second=is_dirty();
    result.first=true;
    init_dirty();
    {
        std::lock_guard<mutex_type> lock(mutex);
        new_inserts.push_back(mut);
    }
    return result;
}
void Tree::finalize(){
    auto& dirty_nodes=this->dirty_nodes;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dirty_nodes.size()),[&dirty_nodes](tbb::blocked_range<size_t> r){
        for (size_t s=r.begin(); s<=r.end(); s++) {
            dirty_nodes[s]->finalize();
        }
    });
}

bool Node::dirty_remove(int pos) {
    auto result = mutations.dirty_remove(pos);
    if (result.first && (!result.second)) {
        tree->dirty_nodes.push_back(this);
    }
    return result.first;
}
bool Node::dirty_insert(const Mutation &mut, char keep_self) {
    auto result = mutations.dirty_insert(mut, keep_self);
    if (result.first && (!result.second)) {
        tree->dirty_nodes.push_back(this);
    }
    return result.first;
}
/*
void Node::finalize(){
    if (not_sample()) {
        mutations.finalize();
    }else {
        std::lock_guard<Mutations_Collection::mutex_type> lock(mutations.mutex);
        if(!mutations.is_dirty()) return;
        auto iter=std::find(parent->children.begin(),parent->children.end(),this);
        Node* common=new Node();
        Node* sample=new Node();
        *iter=common;
        common->children.push_back(this);
        common->children.push_back(sample);
        mutations.dirty_set_difference(common->mutations,sample->mutations);
    }
}
*/

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