#include "tree_rearrangement_internal.hpp"
#include <utility>
void Move_Executor::operator()(Profitable_Moves* in)const{
    if (!in) {
        return;
    }
    if (in->moves.size()>1) {
        multiple_optimal_deferred.push_back(in);
        return;
    }
    
    Move& this_move=in->moves[0];
    //conflicts checked at the last step
    for(auto m:this_move.states){
        Fitch_Sankoff::sankoff_forward_pass(in->shared.range, m->states, dfs_ordered_nodes, m->mutation,m->original_state);
    }
    //Register Move
    ConfirmedMove temp;
    auto op_node=pending_moves.insert(std::make_pair(in->src->parent,temp));
    op_node.first->second.removed.push_back(in->src);
    op_node=pending_moves.insert(std::make_pair(this_move.dst,temp));
    if (this_move.merger.empty()) {
        op_node.first->second.added.push_back(in->src);
    }else {
        op_node.first->second.merges.push_back(std::make_pair(in->src,this_move.merger));
    }
    delete in->src_tip_fs_result;
    delete in;
}
