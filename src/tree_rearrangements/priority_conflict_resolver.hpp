#ifndef Priority_Conflict_Resolver
#define Priority_Conflict_Resolver
#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>
struct Conflict_Set{
    int parsimony_score;
    std::vector<Profitable_Move*> moves;
};
typedef std::unordered_map<MAT::Node*,Conflict_Set,Node_Idx_Hash,Node_Idx_Eq> Cross_t;

struct Conflict_Resolver{
    Cross_t& potential_crosses;
    std::vector<MAT::Node*>& deferred_nodes;
    bool check_single_move_no_conflict(Profitable_Move* candidate_move)const;
    void register_single_move_no_conflict(Profitable_Move* candidate_move)const;
    char operator()(Profitable_Moves_From_One_Source*) const;
};
void schedule_moves(Cross_t& found_moves, std::vector<Profitable_Move*>& out);

#endif