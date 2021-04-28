#ifndef Priority_Conflict_Resolver
#define Priority_Conflict_Resolver
#include "mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <tbb/queuing_rw_mutex.h>
#include <unordered_set>
#include <utility>
#include <vector>

struct Conflict_Set{
    std::atomic<int> parsimony_score_change;
    std::vector<Profitable_Moves_ptr_t> moves;
    Conflict_Set():parsimony_score_change(0){}
};
typedef std::vector<Conflict_Set> Cross_t;

struct Conflict_Resolver{
    Cross_t potential_crosses;
    //int& nodes_inside;
    std::mutex register_lock;
    Conflict_Resolver(size_t node_count):potential_crosses(node_count){}
    bool check_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move)const;
    void register_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move);
    char operator()(std::vector<Profitable_Moves_ptr_t>& candidate_move);
    void schedule_moves(std::vector<Profitable_Moves_ptr_t>& out);
};

#endif