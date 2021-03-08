#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"

typedef std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> Cross_t;
typedef std::unordered_map<MAT::Node*, std::unordered_map<int,Profitable_Move*>,Node_Idx_Hash,Node_Idx_Eq> Mut_t;
struct Conflict_Resolver{
    std::vector<Profitable_Move*>& non_conflicting_moves;
    Cross_t& potential_crosses;
    Mut_t& repeatedly_mutating_loci;
    std::vector<MAT::Node*>& deferred_nodes;
    bool check_single_move_no_conflict(Profitable_Move* candidate_move)const;
    void register_single_move_no_conflict(Profitable_Move* candidate_move)const;
    char operator()(Profitable_Moves_From_One_Source*) const;
};