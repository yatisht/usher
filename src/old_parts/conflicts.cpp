#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <vector>
#include "conflict.hpp"
typedef std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> Cross_t;
typedef std::unordered_map<MAT::Node*, std::unordered_map<int,Profitable_Move*>,Node_Idx_Hash,Node_Idx_Eq> Mut_t;

static Profitable_Move* check_mut_conflict(MAT::Node *node,const std::vector<Fitch_Sankoff_Result_Final> &states,const Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter = repeatedly_mutating_loci.find(node);
        if (iter != repeatedly_mutating_loci.end()) {
            const std::unordered_map<int,Profitable_Move*> &subtree_mut_changed = iter->second;
            for (const Fitch_Sankoff_Result_Final& e : states) {
                auto mut_iter=subtree_mut_changed.find(e.mutation.position);
                if (mut_iter!=subtree_mut_changed.end()) {
                    return mut_iter->second;
                }
            }
        }
        node=node->parent;
    }
    return nullptr;
}
static bool check_loop_conflict(const std::vector<MAT::Node*>& path,const Cross_t& potential_crosses) {
    for (size_t node_idx=0;node_idx<path.size();node_idx++) {
        if (potential_crosses.find(path[node_idx])!=potential_crosses.end()) {
            return true;
        }
    }
    return false;
}
static void register_mut_conflict(MAT::Node *node,
                           const std::vector<Fitch_Sankoff_Result_Final> &states,Profitable_Move* m,Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter =
            repeatedly_mutating_loci.insert({node, std::unordered_map<int,Profitable_Move*>()});
        std::unordered_map<int,Profitable_Move*> &to_insert = iter.first->second;
        for (const Fitch_Sankoff_Result_Final& e : states) {
            to_insert.emplace(e.mutation.position,m);
        }
        node = node->parent;
    }
}
static void register_loop_conflict(const std::vector<MAT::Node*>& path,Cross_t& potential_crosses) {
    for (size_t node_idx=0;node_idx<path.size();node_idx++) {
        potential_crosses.insert(path[node_idx]);
    }
}
bool Conflict_Resolver::check_single_move_no_conflict(Profitable_Move* candidate_move)const{
    const std::vector<MAT::Node*>& path=candidate_move->path;
    if(check_loop_conflict(path,potential_crosses)){
        delete candidate_move;
        return false;
    }
    if (check_mut_conflict(candidate_move->LCA,candidate_move->states,repeatedly_mutating_loci)){
        deferred_nodes.push_back(candidate_move->src);
        delete candidate_move;
        return false;
    }
    return true;
}
void Conflict_Resolver::register_single_move_no_conflict(Profitable_Move* candidate_move) const{
    const std::vector<MAT::Node*>& path=candidate_move->path;
    register_loop_conflict(path,potential_crosses);
    register_mut_conflict(candidate_move->LCA,candidate_move->states,candidate_move,repeatedly_mutating_loci);
    deferred_nodes.push_back(candidate_move->src->parent);
    non_conflicting_moves.push_back(candidate_move);
}

char Conflict_Resolver::operator()(Profitable_Moves_From_One_Source* candidate_move)const{
    if (candidate_move->profitable_moves.empty()) {
        return 0;
    }
    char ret=MOVE_FOUND_MASK;
    for(Profitable_Move* move:candidate_move->profitable_moves){
        if (ret&NONE_CONFLICT_MOVE_MASK) {
            delete move;
        }else {
            if (check_single_move_no_conflict(move)) {
                register_single_move_no_conflict(move);
                ret|=NONE_CONFLICT_MOVE_MASK;
            }
        }
    }
    return ret;
}