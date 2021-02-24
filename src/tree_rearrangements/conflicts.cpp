#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <vector>
typedef std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> Cross_t;
typedef std::unordered_map<MAT::Node*, std::unordered_map<int,Profitable_Move_Deserialized*>,Node_Idx_Hash,Node_Idx_Eq> Mut_t;

static Profitable_Move_Deserialized* check_mut_conflict(MAT::Node *node,const std::vector<Fitch_Sankoff_Result_Deserialized> &states,const Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter = repeatedly_mutating_loci.find(node);
        if (iter != repeatedly_mutating_loci.end()) {
            const std::unordered_map<int,Profitable_Move_Deserialized*> &subtree_mut_changed = iter->second;
            for (const Fitch_Sankoff_Result_Deserialized& e : states) {
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
static bool check_loop_conflict(MAT::Node** path,size_t n_hop,const Cross_t& potential_crosses) {
    for (size_t node_idx=0;node_idx<n_hop;node_idx++) {
        if (potential_crosses.find(path[node_idx])!=potential_crosses.end()) {
            return true;
        }
    }
    return false;
}
static void register_mut_conflict(MAT::Node *node,
                           const std::vector<Fitch_Sankoff_Result_Deserialized> &states,Profitable_Move_Deserialized* m,Mut_t& repeatedly_mutating_loci) {
    while (node) {
        auto iter =
            repeatedly_mutating_loci.insert({node, std::unordered_map<int,Profitable_Move_Deserialized*>()});
        std::unordered_map<int,Profitable_Move_Deserialized*> &to_insert = iter.first->second;
        for (const Fitch_Sankoff_Result_Deserialized& e : states) {
            to_insert.emplace(e.mutation.position,m);
        }
        node = node->parent;
    }
}
static void register_loop_conflict(MAT::Node** path,size_t n_hop,Cross_t& potential_crosses) {
    for (size_t node_idx=0;node_idx<n_hop;node_idx++) {
        potential_crosses.insert(path[node_idx]);
    }
}
void resolve_conflict(Profitable_Moves_Cacher& candidate_moves, std::vector<Profitable_Move_Deserialized*>& non_conflicting_moves, std::vector<MAT::Node*>& deferred_nodes){
    deferred_nodes.clear();
    Cross_t potential_crosses;
    Mut_t repeatedly_mutating_loci;
    while(!candidate_moves.eof()){
        MAT::Node** path;
        size_t n_hop=candidate_moves.get_path(&path);
        if(check_loop_conflict(path,n_hop,potential_crosses)){
            ++candidate_moves;
            continue;
        }
        Profitable_Move_Deserialized* m=*candidate_moves;
        if (check_mut_conflict(m->LCA,m->states,repeatedly_mutating_loci)) {
            deferred_nodes.push_back(m->src);
            delete m;
        }else {
            register_loop_conflict(path,n_hop,potential_crosses);
            register_mut_conflict(m->LCA,m->states,m,repeatedly_mutating_loci);
            deferred_nodes.push_back(m->src->parent);
            non_conflicting_moves.push_back(m);
        }
        ++candidate_moves;
    }
}