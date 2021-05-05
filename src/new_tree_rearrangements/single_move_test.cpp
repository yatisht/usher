#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstdio>
#include <utility>
#include "tbb/parallel_for.h"
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst){
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src=src->parent;
        }
        else if (src->dfs_index<dst->dfs_index) {
            dst=dst->parent;
        }
    }
    return src;
}
int main(int argc,char** argv){
    MAT::Tree old_tree=Mutation_Annotated_Tree::load_mutation_annotated_tree(argv[1]);
    Original_State_t origin_states;
    check_samples(old_tree.root, origin_states, &old_tree);
    auto new_map=new std::unordered_map<std::string, nuc_one_hot>();
    for(const auto& sample:origin_states){
        for(const auto& mut:sample.second){
            auto res=mutated_positions.insert(std::make_pair(mut, new_map));
            if (res.second) {
                new_map=new std::unordered_map<std::string, nuc_one_hot>();
            }
            res.first->second->emplace(sample.first,mut.get_mut_one_hot());
        }
    }
    MAT::Tree t;
    t.load_detatiled_mutations("detailed_mutation_out.pb");
    auto bfs_ordered_nodes=t.breadth_first_expansion();
    fprintf(stderr,"%zu",t.get_parsimony_score());
                        output_t out;
    tbb::concurrent_vector<MAT::Node*> nodes_to_search;
        find_nodes_to_move(bfs_ordered_nodes, nodes_to_search);
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                [&nodes_to_search](tbb::blocked_range<size_t> r) {
                    for (size_t i = r.begin(); i < r.end(); i++) {
                        output_t out;
                        find_profitable_moves(nodes_to_search[i], out,10);
                        if (!out.moves.empty()) {
                            fprintf(stderr, "move %d\n",out.moves[0]->score_change);
                        }else {
                            fprintf(stderr,"%zu no change",i);
                        }
                    }
                });

    MAT::Node* src=t.get_node("Shanghai_SH0007_2020");
    MAT::Node* dst=t.get_node("Wuhan_HBCDC-HB-04_2020");
    MAT::Node* LCA=get_LCA(src, dst);
    individual_move(src,dst,LCA);
    /*individual_move(t.get_node("3034"),t.get_node("EPI_ISL_460417"),t.get_node("3036"));;
    individual_move(t.get_node("5971"),t.get_node("EPI_ISL_469247"),t.get_node("5972"));;
    individual_move(t.get_node("4775"),t.get_node("EPI_ISL_437278"),t.get_node("4777"));;
    individual_move(t.get_node("EPI_ISL_481587"),t.get_node("2022"),t.get_node("2022"));;
    individual_move(t.get_node("1731"),t.get_node("EPI_ISL_476961"),t.get_node("1846"));;
    individual_move(t.get_node("EPI_ISL_452364"),t.get_node("7215"),t.get_node("7215"));;
    individual_move(t.get_node("EPI_ISL_452795"),t.get_node("414"),t.get_node("414"));;
    individual_move(t.get_node("EPI_ISL_455369"),t.get_node("7312"),t.get_node("7312"));;
    individual_move(t.get_node("EPI_ISL_455370"),t.get_node("311"),t.get_node("311"));;*/
}
/*
change of -1 @ 24389 
change of -1 @ 24390 
change of 1 @ 28881 
change of 1 @ 28882 
change of 1 @ 28883 
*/