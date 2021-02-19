#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <vector>
typedef std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> Cross_t;
typedef std::unordered_map<MAT::Node*, std::unordered_map<int,Move*>,Node_Idx_Hash,Node_Idx_Eq> Mut_t;
struct Move_Comparator{
    bool operator()(Move* first,Move* second){
        return first->score_change<second->score_change;
    }
};
static Move* check_mut_conflict(MAT::Node *node,const std::vector<Fitch_Sankoff_Result *> &states,const Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter = repeatedly_mutating_loci.find(node);
        if (iter != repeatedly_mutating_loci.end()) {
            const std::unordered_map<int,Move*> &subtree_mut_changed = iter->second;
            for (const Fitch_Sankoff_Result *e : states) {
                auto mut_iter=subtree_mut_changed.find(e->mutation.position);
                if (mut_iter!=subtree_mut_changed.end()) {
                    return mut_iter->second;
                }
            }
        }
        node=node->parent;
    }
    return nullptr;
}
static bool check_loop_conflict(std::vector<MAT::Node*> path,const Cross_t& potential_crosses) {
    for (auto n:path) {
        if (potential_crosses.find(n)!=potential_crosses.end()) {
            return true;
        }
    }
    return false;
}
static void register_mut_conflict(MAT::Node *node,
                           const std::vector<Fitch_Sankoff_Result *> &states,Move* m,Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter =
            repeatedly_mutating_loci.insert({node, std::unordered_map<int,Move*>()});
        std::unordered_map<int,Move*> &to_insert = iter.first->second;
        for (const Fitch_Sankoff_Result *e : states) {
            to_insert.emplace(e->mutation.position,m);
        }
        node = node->parent;
    }
}
static void register_loop_conflict(std::vector<MAT::Node*> path,Cross_t& potential_crosses) {
    potential_crosses.insert(path.begin(),path.end());
}
void resolve_conflict(tbb::concurrent_vector<Move*>& candidate_moves, std::vector<Move*>& non_conflicting_moves, std::vector<MAT::Node*>& deferred_nodes){
    deferred_nodes.clear();
    Cross_t potential_crosses;
    Mut_t repeatedly_mutating_loci;
    std::sort(candidate_moves.begin(),candidate_moves.end(),Move_Comparator());
    for(auto m:candidate_moves){
        if (check_loop_conflict(m->path,potential_crosses)||check_mut_conflict(m->LCA,m->states,repeatedly_mutating_loci)) {
            deferred_nodes.push_back(m->src);
            for (size_t i=0; i<m->states.size(); i++) {
                delete m->states[i];
            }
            delete m;
        }else {
            register_loop_conflict(m->path,potential_crosses);
            register_mut_conflict(m->LCA,m->states,m,repeatedly_mutating_loci);
            non_conflicting_moves.push_back(m);
        }
    }
}