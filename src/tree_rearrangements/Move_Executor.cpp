#include "tree_rearrangement_internal.hpp"
#include <utility>
void Move_Executor::operator()(tbb::blocked_range<size_t> &r) const {
    for (auto i = r.begin(); i < r.end(); i++) {
        Move *this_move = moves[i];
        // conflicts checked at the last step
        for (size_t j=0;j<this_move->states.size();j++) {
            Fitch_Sankoff_Result* m=this_move->states[j];
            Fitch_Sankoff::sankoff_forward_pass(m->range, m->states,
                                                dfs_ordered_nodes, m->mutation,
                                                m->original_state,m->scores);
            //delete m;
        }
        // Register Move
        ConfirmedMove temp;
        auto op_node = tree_edits.insert(std::make_pair(this_move->src->parent, temp));
        op_node.first->second.removed.push_back(this_move->src);
        op_node = tree_edits.insert(std::make_pair(this_move->dst, temp));
        op_node.first->second.added.push_back(this_move->src);
        delete this_move;
    }
}
