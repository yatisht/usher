#ifndef Priority_Conflict_Resolver
#define Priority_Conflict_Resolver
#include "mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include "tbb/concurrent_vector.h"
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>
//For recording which moves that have path crossing this node is pending to be applied, and its parsimony score improvement
typedef std::vector<Profitable_Moves_ptr_t> Cross_t;
typedef std::vector<std::pair<std::string,std::string>> Deferred_Move_t;
struct Conflict_Resolver {
    Cross_t& potential_crosses;
    //int& nodes_inside;
    Deferred_Move_t & deferred_moves;
    Conflict_Resolver(Cross_t& potential_crosses,Deferred_Move_t& deferred_moves):potential_crosses(potential_crosses),deferred_moves(deferred_moves) {}
    bool check_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move)const;
    bool register_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move) const;
    //enqueuing a move
    char operator()(std::vector<Profitable_Moves_ptr_t>* candidate_move) const;
    //output non-conflicting moves
};
void schedule_moves(Cross_t& potential_crosses,std::vector<Profitable_Moves_ptr_t>& out);

#endif