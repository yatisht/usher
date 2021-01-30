#include "tree_rearrangement_internal.hpp"
#include <algorithm>
Movable_Node_Enumerator::Movable_Node_Enumerator(std::vector<MAT::Node *>& to_check,std::vector<MAT::Node *>& dfs_ordered_nodes):dfs_ordered_nodes(dfs_ordered_nodes){
    for(auto n:to_check){
        n=n->parent;
        if(n)
        this_round.push_back(n->index);
    }
    std::sort(this_round.begin(),this_round.end());
}

Mutation_Annotated_Tree::Node* Movable_Node_Enumerator::operator() (tbb::flow_control fc) const{
    if (iter==this_round.end()) {
        if (next_round.empty()){
            fc.stop();
            return nullptr;
        }
        this_round.swap(next_round);
        iter=this_round.begin();
    }
    while((iter+1)<this_round.end()){
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (!check_grand_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter])) {
            return dfs_ordered_nodes[*iter];
        }else{
            next_round.push_back(*iter);
        }
        iter++;
    }
    assert((iter+1)==this_round.end());
    return dfs_ordered_nodes[*iter];
}
