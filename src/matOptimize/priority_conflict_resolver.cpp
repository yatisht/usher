//If the path of 2 moves intersect, they are considered conflicting
#include "priority_conflict_resolver.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

bool Conflict_Resolver::check_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move)const {
    int best_score=0;
    //gather the minimium parsimony score of all the moves intersecting with this move
    candidate_move->apply_nodes([&best_score,this](MAT::Node* node) {
        if (potential_crosses[node->bfs_index]) {
            best_score=std::min(best_score,potential_crosses[node->bfs_index]->score_change);
        }
    });
    //only insert if its score change is the most negative among all conflicting moves
    if (candidate_move->score_change<best_score) {
        return true;
    }
    return false;
}
//Remove other ove from the set of move to apply, by clearing it from all nodes on its path
// also reset the minimum parsimony score change of all moves whose path include these nodes to 0
//, except for the "exclude" node, as it have been set for the new move.
static void remove_move(Cross_t &potential_crosses, const std::shared_ptr<Profitable_Moves>& other_move,
                        MAT::Node *exclude) {

    other_move->apply_nodes([&](MAT::Node *other_nodes_in_path) {
        if (other_nodes_in_path != exclude) {
            auto& other_node_conflict = potential_crosses[other_nodes_in_path->bfs_index];
            if (other_node_conflict) {
                other_node_conflict=nullptr;
            }
        }
    });
}
// Reggister "candidate_move" to apply
bool Conflict_Resolver::register_single_move_no_conflict (
    Profitable_Moves_ptr_t& candidate_move) const {
    candidate_move->apply_nodes([&](MAT::Node* node) {
        //for each node on the path
        auto& this_node_move = potential_crosses[node->bfs_index];
        //clear all the conflicting moves
        if (this_node_move) {
            remove_move(potential_crosses, this_node_move, node);
        }
        this_node_move=candidate_move;
        assert(potential_crosses[node->bfs_index]==this_node_move);
    });
    return true;
}

char Conflict_Resolver::operator()(std::vector<Profitable_Moves_ptr_t>* candidate_move_ptr) const {
    std::vector<Profitable_Moves_ptr_t>& candidate_move=*candidate_move_ptr;
    char ret=0;
    Profitable_Moves_ptr_t selected_move=nullptr;
    for (Profitable_Moves_ptr_t& move : candidate_move) {
        move->populate_involved_nodes();
        //don't need check-lock-check-set, as there is little contension on conflict resolver
        if (check_single_move_no_conflict(move)) {
            //fprintf(stderr, "registered move\n");
            register_single_move_no_conflict(move);
            ret =1;
            selected_move = move;
            break;
        }
    }

    if(!selected_move&&(!candidate_move.empty())) {
        for (Profitable_Moves_ptr_t move : candidate_move) {
            deferred_moves.emplace_back(move->src->identifier,move->dst->identifier);
        }
    }
    delete candidate_move_ptr;
    return ret;
}
//output all the nodes for apply
void schedule_moves(Cross_t& potential_crosses, std::vector<Profitable_Moves_ptr_t>& out) {
#ifdef CONFLICT_RESOLVER_DEBUG
    std::unordered_map<size_t,std::pair<size_t,Profitable_Moves_ptr_t>> pushed_moves;
#endif
    //go from leaf to root (reverse bfs order), so will get the terminal of the move (src and dst) first,
    //before most of the other
    for (auto& this_move:potential_crosses) {
        if (this_move) {
                out.push_back(this_move);
                //remove it to prevent it from enqueued again
                remove_move(potential_crosses, this_move, 0);
        }
    }
#ifndef NDEBUG
    for(const auto& a:potential_crosses) {
        assert(!a);
    }
#endif
}
