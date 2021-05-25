#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <vector>
char merge_new_node_mutations(
    const MAT::Mutations_Collection &new_node_mutations,
    const MAT::Mutations_Collection &sibling_node_mutations,
    MAT::Mutations_Collection &shared_node_mutations_out,
    MAT::Mutations_Collection &sibling_node_mutations_out,
    MAT::Mutations_Collection &new_node_mutations_out,
    MAT::Mutations_Collection &to_merge_if_children);

static void find_split(std::vector<MAT::Node*>& all_nodes,std::vector<MAT::Node*>& to_split){
    for(auto node:all_nodes){
        
    }
}