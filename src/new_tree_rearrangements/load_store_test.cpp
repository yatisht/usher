#include "src/new_tree_rearrangements/check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <vector>
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
int main(int argc, char** argv){
    Original_State_t o;
    MAT::Tree ori_tree=load_tree(argv[1],o);
    ori_tree.save_detailed_mutations("test_save.pb");

    MAT::Tree new_tree;
    new_tree.load_detatiled_mutations("test_save.pb");
    std::vector<MAT::Node*> old_nodes=ori_tree.depth_first_expansion();
    std::vector<MAT::Node*> new_nodes=new_tree.depth_first_expansion();
    for(size_t node_idx=0;node_idx<old_nodes.size();node_idx++){
        assert(old_nodes[node_idx]->identifier==new_nodes[node_idx]->identifier);
        assert((!node_idx)||old_nodes[node_idx]->parent->identifier==new_nodes[node_idx]->parent->identifier);
        for(size_t mut_idx=0;mut_idx<new_nodes[node_idx]->mutations.size();mut_idx++){
            assert(old_nodes[node_idx]->mutations[mut_idx]==new_nodes[node_idx]->mutations[mut_idx]);
        }
    }
}