#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>
MAT::Node* Move_Executor::get_parent(MAT::Node* node) const{
        auto iter=new_parents_map.find(node);
        if(iter==new_parents_map.end()){
            return node->parent;
        }else{
            return (MAT::Node*)iter->second;
        }
}
static size_t get_parsimony_score(MAT::Node* root){
    int par_score=root->mutations.size();
    for(auto child:root->children){
        par_score+=get_parsimony_score(child);
    }
    return par_score;
}
int count_mutation(MAT::Node* node,int pos){
    unsigned int count=0;
    if (node) {
        count=std::count_if(node->mutations.begin(), node->mutations.end(),[pos](const MAT::Mutation& m){return m.position==pos;} );
    }
    assert(count<=1);
    return count;
}
/*
#ifndef NDEBUG
void check_original_states(Fitch_Sankoff::States_Type& original_states, Mutation_Annotated_Tree::Mutation& mutation,const Original_State_t& mut){
    for(Fitch_Sankoff::State_Type& state:original_states){
        auto iter=mut.find(state.node->identifier);
        if (iter!=mut.end()) {
            assert(state.node->is_leaf());
            auto mut_iter=iter->second.find(mutation);
            if(mut_iter!=iter->second.end()){
                assert(state.state==mut_iter->mut_nuc);
            }else {
                assert(state.state==mutation.ref_nuc);
            }
        }else{
            assert(!state.node->is_leaf());
        }
    }
}
#endif
*/
void Move_Executor::operator()(tbb::blocked_range<size_t> &r) const {
    for (auto i = r.begin(); i < r.end(); i++) {
        Profitable_Move *this_move = moves[i];

        MAT::Node* new_leaf=nullptr;
        #ifndef NDEBUG
        for(auto m:tree_edits){
            const auto& other_removed=m.second.removed;
            assert(std::find(other_removed.begin(),other_removed.end(),this_move->src)==other_removed.end());
        }
        #endif
        // Register Move
        ConfirmedMove temp;
        auto op_node = tree_edits.insert(std::make_pair(this_move->src->parent, temp));
        op_node.first->second.removed.push_back(this_move->src);
        assert(op_node.first->first==this_move->src->parent);
        op_node = tree_edits.insert(std::make_pair(this_move->dst, temp));
        std::vector<MAT::Node*>& to_add=op_node.first->second.added;
        if(this_move->dst->is_leaf()){
            if(to_add.empty()){
                new_leaf=new MAT::Node();
                new_leaf->parent=this_move->dst;
                new_leaf->identifier=this_move->dst->identifier;
                to_add.push_back(new_leaf);
            }else{
                new_leaf=to_add.front();
                assert(new_leaf->identifier==this_move->dst->identifier);
            }
        }
        to_add.push_back(this_move->src);

        assert(new_parents_map.insert(std::make_pair(this_move->src, this_move->dst)).second);
        #ifndef NDEBUG
        int ori_score=get_parsimony_score(this_move->LCA);
        #endif
        // conflicts checked at the last step
        for (size_t j=0;j<this_move->states.size();j++) {
            Fitch_Sankoff_Result_Final& m=this_move->states[j];
            Fitch_Sankoff::sankoff_forward_pass(
                this_move->range, dfs_ordered_nodes, m.mutation, ori,
                m.scores, m.LCA_parent_state, this_move->src, this_move->dst,new_leaf);
            #ifndef NDEBUG
            int this_mut_score=count_mutation(new_leaf, m.mutation.position);
            for (size_t node_idx=this_move->range.first; node_idx<this_move->range.second; node_idx++) {
                this_mut_score+=count_mutation(dfs_ordered_nodes[node_idx], m.mutation.position);
            }
            assert(this_mut_score==m.optimized_score);
            #endif
        }
#ifndef NDEBUG
        int new_score=get_parsimony_score(this_move->LCA);
        Original_State_t copy(ori);
        Mutation_Set parental;
        MAT::Node* ancestor=get_parent(this_move->LCA);
        while(ancestor){
            for(const auto& m:ancestor->mutations){
                    parental.insert(m);
            }
            ancestor=get_parent(ancestor);
        }
        for (auto iter = parental.begin(); iter != parental.end();) {
            if (iter->mut_nuc == iter->ref_nuc) {
                auto prev_iter = iter;
                iter++;
                parental.erase(prev_iter);
            } else {
                iter++;
            }
        }
        check_samples_worker_with_pending_moves(this_move->LCA, parental, copy,tree_edits);
        int new_leaf_size=new_leaf?new_leaf->mutations.size():0;
        if (this_move->score_change<(new_score+new_leaf_size-ori_score)) {
            std::vector<Profitable_Move*> prev;
            for (size_t idx=0; idx<i; idx++) {
                Profitable_Move* other_move=moves[idx];
                if (other_move->range.first<=this_move->range.first&&other_move->range.second>this_move->range.first) {
                    assert(other_move->range.second>=this_move->range.second);
                    prev.push_back(other_move);
                }
            }
        assert(false);
        }
#endif
        //delete this_move;
    }
}
