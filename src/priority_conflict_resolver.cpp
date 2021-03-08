#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <vector>
struct Conflict_Set{
    int parsimony_score;
    std::vector<Profitable_Move*> moves;
};
typedef std::unordered_map<MAT::Node*,Conflict_Set,Node_Idx_Hash,Node_Idx_Eq> Cross_t;

struct Conflict_Resolver{
    std::vector<Profitable_Move*>& non_conflicting_moves;
    Cross_t& potential_crosses;
    Cross_t& LCA;
    std::vector<MAT::Node*>& deferred_nodes;
    bool check_single_move_no_conflict(Profitable_Move* candidate_move)const;
    void register_single_move_no_conflict(Profitable_Move* candidate_move)const;
    char operator()(Profitable_Moves_From_One_Source*) const;
};
bool Conflict_Resolver::check_single_move_no_conflict(Profitable_Move* candidate_move)const{
    int best_score=0;
    for(auto node:candidate_move->path){
        auto iter=potential_crosses.find(node);
        if (iter!=potential_crosses.end()) {
            best_score=std::min(best_score,iter->second.parsimony_score);
        }
    }
    if (candidate_move->score_change<best_score) {
        return true;
    }
    return false;
}


