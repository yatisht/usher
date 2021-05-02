#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <unordered_map>
struct Changed_Node{
    std::unordered_map<int,std::unordered_map<size_t,uint8_t>> altered_state;
    Profitable_Moves_ptr_t move;
    bool set;
    Changed_Node():move(nullptr),set(false){}
};
static void apply_change_in_range(const std::vector<MAT::Node*> dfs_ordered_nodes,std::vector<Changed_Node>& changed_nodes, size_t start, size_t end){
        for (int node_idx=end; node_idx>=start; node_idx--) {
        Changed_Node& this_node_change=changed_nodes[node_idx];
        MAT::Node* this_node=dfs_ordered_nodes[node_idx];
        if (!this_node_change.set) {
            Profitable_Moves_ptr_t move_on_this_node=
            if (this_node_change.move) {
                MAT::Node* other_node=
            }
        }
        this_node_change.set=true;
    }
}
void apply_moves(std::vector<Profitable_Moves_ptr_t> & all_moves, MAT::Tree &t, tbb::concurrent_vector<MAT::Node *> &to_filter){
    std::vector<MAT::Node*> dfs_ordered_nodes=t.depth_first_expansion();
    std::vector<Changed_Node> changed_nodes;
    for(Profitable_Moves_ptr_t move:all_moves){
        assert(!changed_nodes[move->get_src()->dfs_index].move);
        assert(!changed_nodes[move->get_dst()->dfs_index].move);
        changed_nodes[move->get_src()->dfs_index].move=move;
        changed_nodes[move->get_dst()->dfs_index].move=move;
    }

}