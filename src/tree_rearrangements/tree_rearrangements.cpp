#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/pipeline.h>
#include "check_samples.hpp"
#include <tbb/flow_graph.h>
#include "Twice_Bloom_Filter.hpp"
#include "src/tree_rearrangement.hpp"
tbb::concurrent_vector<MAT::Node*> postponed;


static void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output){
    Twice_Bloom_Filter filter;
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            filter.insert(m.position);
        }
    }
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            if(filter.query(m.position)){
                output.push_back(n);
                break;
            }
        }
    }
}

static void feed_nodes(std::vector<MAT::Node *> &to_feed,
                         tbb::flow::buffer_node<MAT::Node*>& input_node,
                         std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;

    {
        this_round_idx.reserve(to_feed.size());
        for (auto node : to_feed) {
            auto parent = node->parent;
            if (parent&&parent->parent)
                this_round_idx.push_back(parent->index);
        }
    }
    while (!this_round_idx.empty()) {
    std::sort(this_round_idx.begin(), this_round_idx.end());
    for (auto iter = this_round_idx.begin(); iter < this_round_idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (!check_grand_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter])) {
            MAT::Node* this_node=dfs_ordered_nodes[*iter];
            input_node.try_put(this_node);
            this_node=this_node->parent;
            if(this_node->parent){
            next_round_idx.push_back(this_node->index);
            }
        } else {
            next_round_idx.push_back(*iter);
        }
    }
    input_node.try_put(dfs_ordered_nodes[this_round_idx.back()]);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    }
}


void Tree_Rearrangement::refine_trees(std::vector<MAT::Tree> &optimal_trees,int radius) {

    for (auto this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Sample_Mut_Type ori;
        std::vector<MAT::Node*>& to_optimize=this_tree.new_nodes;
        find_nodes_with_recurrent_mutations(dfs_ordered_nodes, to_optimize);
        check_samples(this_tree.root, ori);
        Pending_Moves_t pending_moves;
        tbb::concurrent_vector<Move*> profitable_moves;
        //building pipeline
        tbb::flow::graph search_graph;
        tbb::flow::buffer_node<MAT::Node*> input(search_graph);
        Neighbors_Finder_t neighbors_finder(search_graph,tbb::flow::unlimited,Neighbors_Finder{radius});
        tbb::flow::make_edge(input,neighbors_finder);
        tbb::flow::interface11::function_node<Possible_Move*> profitable_move_enumerator(search_graph,tbb::flow::unlimited,Profitable_Moves_Enumerator{dfs_ordered_nodes,profitable_moves});
        tbb::flow::make_edge(neighbors_finder,profitable_move_enumerator);


        while (!to_optimize.empty()) {
            profitable_moves.clear();
            pending_moves.clear();
            feed_nodes(to_optimize,input,dfs_ordered_nodes);
            search_graph.wait_for_all();
            std::vector<Move *> non_conflicting_moves;
            resolve_conflict(profitable_moves, non_conflicting_moves, to_optimize);
            if(!non_conflicting_moves.empty()){
                Pending_Moves_t tree_edits;
                // tbb::parallel_for(tbb::blocked_range<size_t>(0,
                // non_conflicting_moves.size()),
                // Move_Executor{dfs_ordered_nodes,this_tree,non_conflicting_moves,tree_edits});
                Move_Executor temp{dfs_ordered_nodes, this_tree,
                                   non_conflicting_moves, tree_edits, ori};
                tbb::blocked_range<size_t> range_temp(
                    0, non_conflicting_moves.size());
                temp(range_temp);
                // this_tree.finalize();
                /*
                tbb::parallel_for_each(
                    tree_edits.begin(), tree_edits.end(),
                    [&this_tree,
                     &ori](const std::pair<MAT::Node *, ConfirmedMove> &in) {
                   });
                */
                for (Pending_Moves_t::const_iterator iter=tree_edits.begin(); iter!=tree_edits.end(); iter++) {
                    finalize_children(reinterpret_cast<MAT::Node *>(iter->first),
                                      const_cast<ConfirmedMove &>(iter->second),
                                      &this_tree, ori);
                }
                dfs_ordered_nodes = this_tree.depth_first_expansion();
#ifndef NDEBUG
                Sample_Mut_Type copy(ori);
                check_samples(this_tree.root, copy);
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
#endif
            }
        }

        check_samples(this_tree.root, ori);
        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        this_tree.reassign_level();
    }
}
