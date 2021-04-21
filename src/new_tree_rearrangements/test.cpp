#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "priority_conflict_resolver.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <cstdio>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/partitioner.h>
#include <unordered_set>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;

int main(int argc, char **argv) {
    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t =load_tree(argv[1],origin_states);
    t.save_detailed_mutations("detailed_mutation_out.pb");
    size_t score_before = t.get_parsimony_score();
    size_t new_score=score_before;
    fprintf(stderr, "after state reassignment:%zu\n", score_before);

    int improvement = -1;
    int stalled=0;
    while (improvement < 0) {
        tbb::concurrent_vector<MAT::Node*> nodes_to_search;
        std::vector<MAT::Node*> bfs_ordered_nodes;
        bfs_ordered_nodes = t.breadth_first_expansion();
        size_t inner_loop_score_before=score_before;
        find_nodes_to_move(bfs_ordered_nodes, nodes_to_search);
    //individual_move(bfs_ordered_nodes[1], bfs_ordered_nodes[63],bfs_ordered_nodes[0]);
    //individual_move(bfs_ordered_nodes[2], bfs_ordered_nodes[0],bfs_ordered_nodes[0]);
    //individual_move(bfs_ordered_nodes[3], bfs_ordered_nodes[0],bfs_ordered_nodes[0]);
        while (!nodes_to_search.empty()) {
            //*
            bfs_ordered_nodes = t.breadth_first_expansion();
            std::unordered_set<std::string> pushed_node;

            fprintf(stderr, "Node size: %zu\n",bfs_ordered_nodes.size());
            fprintf(stderr, "Internal node size %zu\n",t.curr_internal_node);
                /*tbb::parallel_for(
                tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                [&nodes_to_search,&bfs_ordered_nodes](tbb::blocked_range<size_t> r) {
                    for (size_t i = r.begin(); i < r.end(); i++) {
                        output_t out;
                        individual_move(nodes_to_search[i], bfs_ordered_nodes[0],bfs_ordered_nodes[0]);
                    }
                });*/
            tbb::concurrent_vector<MAT::Node *> deferred_nodes;
            Conflict_Resolver resolver(bfs_ordered_nodes.size());
                        output_t out;
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                [&nodes_to_search, &resolver,&deferred_nodes](tbb::blocked_range<size_t> r) {
                    for (size_t i = r.begin(); i < r.end(); i++) {
                        output_t out;
                        find_profitable_moves(nodes_to_search[i], out,5);
                        if (!out.moves.empty()) {
                            deferred_nodes.push_back(out.moves[0]->get_src());
                            resolver(out.moves);
                        }
                    }
                });
            fprintf(stderr,"%zu nodes deferred \n",deferred_nodes.size());
            std::vector<Profitable_Moves_ptr_t> all_moves;
            resolver.schedule_moves(all_moves);
            apply_moves(all_moves, t, bfs_ordered_nodes,deferred_nodes);
            check_samples(t.root, origin_states, &t);
            new_score = t.get_parsimony_score();
            fprintf(stderr, "after optimizing:%zu\n", new_score);
            t.save_detailed_mutations("detailed_mutation_out.pb");
            /*{
                MAT::Tree to_match;
                bfs_ordered_nodes=t.breadth_first_expansion();
                to_match.load_detatiled_mutations("detailed_mutation_out.pb");
                auto new_bfs_ordered_nodes=to_match.breadth_first_expansion();
                assert(new_bfs_ordered_nodes.size()==bfs_ordered_nodes.size());
                for(size_t i=0;i<bfs_ordered_nodes.size();i++){
                    MAT::Node* old_node=bfs_ordered_nodes[i];
                    MAT::Node* new_node=new_bfs_ordered_nodes[i];
                    assert(old_node->identifier==new_node->identifier);
                    assert(old_node->mutations.size()==new_node->mutations.size());
                    for(size_t mut_idx=0;mut_idx<old_node->mutations.size();mut_idx++){
                        assert(old_node->mutations[mut_idx]==new_node->mutations[mut_idx]);
                    }
                }
            }*/
            nodes_to_search = std::move(deferred_nodes);
            fprintf(stderr,"%zu nodes to search \n",nodes_to_search.size());
            if(new_score>=inner_loop_score_before){
                stalled++;
            }else{
                inner_loop_score_before=new_score;
                stalled=0;
            }
            if (stalled>=10) {
                break;
            }
        }
        improvement =  new_score-score_before;
        score_before=new_score;
    }
    std::vector<MAT::Node*> dfs=t.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dfs.size()),[& dfs](tbb::blocked_range<size_t> r){
        for (size_t i=r.begin(); i<r.end(); i++) {
            dfs[i]->mutations.remove_invalid();
        }
    });
    fix_condensed_nodes(&t);
    check_samples(t.root, origin_states, &t);
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, argv[2]);
}