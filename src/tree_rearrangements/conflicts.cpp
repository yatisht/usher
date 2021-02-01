#include "conflicts.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <vector>
bool check_mut_conflict(MAT::Node *node,
                        const std::vector<Fitch_Sankoff_Result *> &states) {
    auto iter = repeatedly_mutating_loci.find(node);
    if (iter == repeatedly_mutating_loci.end()) {
        return false;
    }
    const std::unordered_set<int> &subtree_mut_changed = iter->second;
    for (const Fitch_Sankoff_Result *e : states) {
        if (subtree_mut_changed.count(e->mutation.position)) {
            return true;
        }
    }
    return false;
}
bool check_loop_conflict(MAT::Node *src, MAT::Node *dst) {
    auto iter = potential_crosses.find(src);
    if (iter == potential_crosses.end()) {
        return false;
    }
    std::unordered_set<MAT::Node *> &nodes_moved_to_subtree = iter->second;
    while (dst) {
        if (nodes_moved_to_subtree.count(dst)) {
            return true;
        }
        dst = dst->parent;
    }
    return false;
}
void register_mut_conflict(MAT::Node *node,
                           const std::vector<Fitch_Sankoff_Result *> &states) {
    while (node) {
        auto iter =
            repeatedly_mutating_loci.insert({node, std::unordered_set<int>()});
        std::unordered_set<int> &to_insert = iter.first->second;
        for (const Fitch_Sankoff_Result *e : states) {
            to_insert.insert(e->mutation.position);
        }
        node = node->parent;
    }
}
void register_loop_conflict(MAT::Node *src, MAT::Node *dst) {
    while (dst) {
        auto iter =
            potential_crosses.emplace(dst, std::unordered_set<MAT::Node *>());
        iter.first->second.insert(src);
        dst = dst->parent;
    }
}
Profitable_Moves *
multiple_moves_resolver::operator()(tbb::flow_control fc) const {
    assert(valid);
    while (iter != multiple_optimal_deferred.end()) {
        Profitable_Moves *this_move = *iter;
        size_t min_level = 0xffffff;
        Move *min_leve_move = nullptr;
        MAT::Node *src = this_move->src;
        for (Move &m : this_move->moves) {
            if (m.dst->level < min_level) {
                if (!(check_mut_conflict(src, m.states) ||
                      check_mut_conflict(m.dst, m.states) ||
                      check_loop_conflict(src, m.dst))) {
                    min_leve_move = &m;
                    min_level = m.dst->level;
                }
            }
        }
        if(min_leve_move){
            std::vector<Move> temp({*min_leve_move});
            this_move->moves.swap(temp);
            return this_move;
        }else{
            postponed.push_back(src);
        }
    }
    fc.stop();
    return nullptr;
}