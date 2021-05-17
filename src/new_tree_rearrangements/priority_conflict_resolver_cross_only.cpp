#include "priority_conflict_resolver.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <mutex>
#include <utility>
#include <vector>

bool Conflict_Resolver::register_single_move_no_conflict(
    Profitable_Moves_ptr_t &candidate_move) {
    int best_score = 0;
    struct Conflicting_Move_t {
        size_t bfs_idx;
        std::vector<std::vector<Profitable_Moves_ptr_t>::iterator>
            conflicting_moves;
        Conflicting_Move_t(
            MAT::Node *node,
            const std::vector<Profitable_Moves_ptr_t>::iterator &iter)
            : bfs_idx(node->bfs_index), conflicting_moves({iter}) {}
    };
    std::vector<Conflicting_Move_t> conflicting_moves;
    for (auto node : candidate_move->dst_to_LCA) {
        auto &temp = potential_crosses[node->bfs_index].moves;

        for (std::vector<Profitable_Moves_ptr_t>::iterator iter = temp.begin();
             iter < temp.end(); iter++) {
            auto &move = *iter;
            if (move->LCA == candidate_move->LCA) {
                const auto &other_dst_to_LCA = move->dst_to_LCA;
                if (std::find(other_dst_to_LCA.begin(), other_dst_to_LCA.end(),
                              move->src) != other_dst_to_LCA.end()) {
                    best_score = std::min(best_score, move->score_change);
                    if (conflicting_moves.empty() ||
                        conflicting_moves.back().bfs_idx != node->bfs_index) {
                        conflicting_moves.emplace_back(node, iter);
                    } else {
                        conflicting_moves.back().conflicting_moves.push_back(
                            iter);
                    }
                }
            }
        }
    }
    if (candidate_move->score_change < best_score) {
        for (auto temp : conflicting_moves) {
            auto &to_shrink = potential_crosses[temp.bfs_idx].moves;
            auto exchange_iter = to_shrink.end() - 1;
            for (auto iter : temp.conflicting_moves) {
                std::iter_swap(iter, exchange_iter);
                exchange_iter--;
            }
            exchange_iter++;
            to_shrink.erase(exchange_iter, to_shrink.end());
        }
        potential_crosses[candidate_move->src->bfs_index].moves.push_back(
            candidate_move);
        return true;
    }
    return false;
}

char Conflict_Resolver::operator()(
    std::vector<Profitable_Moves_ptr_t> &candidate_move) {
    char ret = 0;
    Profitable_Moves_ptr_t selected_move = nullptr;
    for (Profitable_Moves_ptr_t &move : candidate_move) {
#ifndef SINGLE_THREAD_TEST
        std::lock_guard<std::mutex> lk(register_lock);
#endif
        register_single_move_no_conflict(move);
        ret = 1;
        selected_move = move;
        break;
    }

#ifndef NDEBUG
    if (selected_move) {
        // nodes_inside++;
    }
#endif
    for (Profitable_Moves_ptr_t move : candidate_move) {
        if (move != selected_move) {
#ifdef CHECK_LEAK
            move->destructed = true;
#else
            delete move;
#endif
        }
    }
    return ret;
}

void Conflict_Resolver::schedule_moves(
    std::vector<Profitable_Moves_ptr_t> &out) {
    for (const auto &node : potential_crosses) {
        out.insert(out.end(), node.moves.begin(), node.moves.end());
    }

}
