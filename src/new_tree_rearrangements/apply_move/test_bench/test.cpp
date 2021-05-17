#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include "src/new_tree_rearrangements/priority_conflict_resolver.hpp"
#include <cstddef>
#include <random>
#include <chrono>
#include <vector>
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
static Profitable_Moves_ptr_t make_move(MAT::Node* src,MAT::Node* dst){
    std::vector<MAT::Node*> src_to_LCA;
    std::vector<MAT::Node*> dst_to_LCA;
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src_to_LCA.push_back(src);
            src=src->parent;
        }
        else if (src->dfs_index<dst->dfs_index) {
            dst_to_LCA.push_back(dst);
            dst=dst->parent;
        }
    }
    dst_to_LCA.push_back(src);
    Profitable_Moves_ptr_t ret=new Profitable_Moves;
    ret->dst_to_LCA=std::move(dst_to_LCA);
    ret->src_to_LCA=std::move(src_to_LCA);
    ret->LCA=src;
    ret->src=ret->src_to_LCA[0];
    return ret;
}
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_PRIMARY_MOVE
                 ,
                 Original_State_t original_state
#endif
) ;
int main(int argc, char** argv){
    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t=load_tree(argv[1], origin_states);
    Profitable_Moves_ptr_t move=make_move(t.get_node("node_56_condensed_2_leaves"), t.get_node("117"));
    std::vector<Profitable_Moves_ptr_t> all_moves{move};
    std::vector<MAT::Node*> bfs_ordered_nodes=t.breadth_first_expansion();
    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes,origin_states);

    //From https://www.cplusplus.com/reference/random/mersenne_twister_engine/mersenne_twister_engine/
    //unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed1 = 50;
    std::mt19937 g1 (seed1);
    while(true){
    all_moves.clear();
    std::vector<MAT::Node*> dfs_ordered_nodes=t.depth_first_expansion();
    t.breadth_first_expansion();
    std::uniform_int_distribution<int> src_idx_dist(1,dfs_ordered_nodes.size()-1);
    std::uniform_int_distribution<int> score_change_idx_dist(-5,-1);
    Conflict_Resolver resolver(dfs_ordered_nodes.size());
    for(int i=0;i<1;i++){
        size_t src_idx=src_idx_dist(g1);
        MAT::Node* src_node=dfs_ordered_nodes[src_idx];
        size_t dst_range=dfs_ordered_nodes.size()-src_node->dfs_end_index+src_idx;
        size_t dst_idx_raw=g1()%dst_range;
        if (dst_idx_raw>=src_idx) {
            dst_idx_raw+=(src_node->dfs_end_index-src_idx);
        }
        if(dst_idx_raw==src_node->parent->dfs_index){
            i--;
            continue;
        }
        Profitable_Moves_ptr_t move=make_move(src_node, dfs_ordered_nodes[dst_idx_raw]);
        move->score_change=score_change_idx_dist(g1);
        std::vector<Profitable_Moves_ptr_t>temp{move};
        resolver(temp);
    }
    resolver.schedule_moves(all_moves);
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes,origin_states);}
}