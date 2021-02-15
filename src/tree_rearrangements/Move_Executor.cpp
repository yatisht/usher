#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <utility>
MAT::Node* Move_Executor::get_parent(MAT::Node* node) const{
        auto iter=new_parents_map.find((void*)node);
        if(iter==new_parents_map.end()){
            return node->parent;
        }else{
            return (MAT::Node*)iter->second;
        }
}

void check_original_states(Fitch_Sankoff::States_Type original_states, Mutation_Annotated_Tree::Mutation& mutation,const Sample_Mut_Type& mut){
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

void Move_Executor::operator()(tbb::blocked_range<size_t> &r) const {
    for (auto i = r.begin(); i < r.end(); i++) {
        Move *this_move = moves[i];

        MAT::Node* new_leaf=nullptr;
        for(auto m:tree_edits){
            const auto& other_removed=m.second.removed;
            assert(std::find(other_removed.begin(),other_removed.end(),this_move->src)==other_removed.end());
        }
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

        assert(new_parents_map.insert(std::make_pair((void*)this_move->src, (void*)this_move->dst)).second);
        // conflicts checked at the last step
        for (size_t j=0;j<this_move->states.size();j++) {
            Fitch_Sankoff_Result* m=this_move->states[j];
            check_original_states(m->original_state, m->mutation, ori);
            Fitch_Sankoff::sankoff_forward_pass(
                m->range, dfs_ordered_nodes, m->mutation, m->original_state,
                m->scores, m->LCA_parent_state, this_move->src, this_move->dst,new_leaf);
            delete m;
        }

        Sample_Mut_Type copy(ori);
        Mutation_Set parental;
        MAT::Node* ancestor=get_parent(this_move->LCA);
        while(ancestor){
            parental.insert(ancestor->mutations.begin(),ancestor->mutations.end());
            ancestor=get_parent(ancestor);
        }
        check_samples_worker_with_pending_moves(this_move->LCA, parental, copy,tree_edits);

        delete this_move;
    }
}
