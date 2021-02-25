#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/pipeline.h>
#include "check_samples.hpp"
#include <tbb/flow_graph.h>
#include <tbb/queuing_rw_mutex.h>
#include <utility>
#include <vector>
#include "Twice_Bloom_Filter.hpp"
#include "src/tree_rearrangement.hpp"
tbb::concurrent_vector<MAT::Node*> postponed;
extern uint32_t num_cores;
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
static void fix_condensed_nodes(MAT::Tree* tree){
    std::vector<MAT::Node*> nodes_to_fix;
    for(auto iter:tree->all_nodes){
        if (tree->condensed_nodes.count(iter.first)&&(!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for(auto node:nodes_to_fix){
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier, std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier,node);
    }
}
static void push_all_nodes(MAT::Tree* tree,std::vector<MAT::Node*>& nodes){
    nodes.clear();
    for(auto a:tree->all_nodes){
        nodes.push_back(a.second);
    }
}
static void feed_nodes(std::vector<MAT::Node *> &to_feed,
                         tbb::flow::buffer_node<MAT::Node*>& input_node,
                         std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;

#define output_node(this_node) \
input_node.try_put(this_node);\
this_node=this_node->parent;\
if(this_node->parent){\
    next_round_idx.push_back(this_node->index);\
}

    {
        this_round_idx.reserve(to_feed.size());
        for (auto node : to_feed) {
            auto parent = node->parent;
            if (parent)
                this_round_idx.push_back(node->index);
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
            output_node(this_node);
        } else {
            next_round_idx.push_back(*iter);
        }
    }
    MAT::Node* last_node=dfs_ordered_nodes[this_round_idx.back()];
    output_node(last_node);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    }
}


void Tree_Rearrangement::refine_trees(std::vector<MAT::Tree> &optimal_trees,int radius) {

    for (MAT::Tree& this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Original_State_t ori;
        std::vector<MAT::Node*>& to_optimize=this_tree.new_nodes;
        check_samples(this_tree.root, ori,&this_tree);
        Pending_Moves_t pending_moves;
        tbb::concurrent_vector<Profitable_Move*> profitable_moves;
        //building pipeline
        tbb::queuing_rw_mutex mutex;
        tbb::flow::graph search_graph;
        tbb::flow::buffer_node<MAT::Node*> input(search_graph);
        tbb::flow::function_node<MAT::Node*,Possible_Moves*> neighbors_finder(search_graph,tbb::flow::unlimited,Neighbors_Finder(radius));
        tbb::flow::make_edge(input,neighbors_finder);
        
        tbb::flow::function_node<Possible_Moves*,Candidate_Moves*> parsimony_score_calculator(search_graph,tbb::flow::unlimited,Parsimony_Score_Calculator{ori,dfs_ordered_nodes});
        tbb::flow::make_edge(neighbors_finder,parsimony_score_calculator);
        
        tbb::flow::function_node<Candidate_Moves*,tbb::flow::continue_msg,tbb::flow::rejecting> profitable_move_enumerator(search_graph,num_cores,Profitable_Moves_Enumerator{dfs_ordered_nodes,profitable_moves,mutex,ori});
        tbb::flow::make_edge(parsimony_score_calculator,profitable_move_enumerator);

        bool have_improvement=true;
        while (have_improvement) {
            have_improvement=false;
            find_nodes_with_recurrent_mutations(dfs_ordered_nodes, to_optimize);
            //push_all_nodes(&this_tree, to_optimize);
        while (!to_optimize.empty()) {
            Profitable_Moves_Cacher cacher(profitable_moves,mutex);
            cacher.run();
            profitable_moves.clear();
            pending_moves.clear();
            feed_nodes(to_optimize,input,dfs_ordered_nodes);
            search_graph.wait_for_all();
            std::vector<Profitable_Move_Deserialized *> non_conflicting_moves;
            cacher.finish();
            to_optimize.clear();
            if(!cacher.eof()){
                resolve_conflict(cacher, non_conflicting_moves, to_optimize);
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
                std::vector<std::pair<MAT::Node*, MAT::Node*>> deleted_map;
                for (Pending_Moves_t::const_iterator iter=tree_edits.begin(); iter!=tree_edits.end(); iter++) {
                    finalize_children(reinterpret_cast<MAT::Node *>(iter->first),
                                      const_cast<ConfirmedMove &>(iter->second),
                                      &this_tree, ori,deleted_map);
                }
                for(MAT::Node*& next_node:to_optimize){
                    for(const auto& deleted:deleted_map){
                        if(next_node==deleted.first){
                            next_node=deleted.second;
                        }
                    }
                }
                dfs_ordered_nodes = this_tree.depth_first_expansion();
#ifndef NDEBUG
                Original_State_t copy(ori);
                check_samples(this_tree.root, copy,&this_tree);
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
#endif
            have_improvement=true;
            }
        }
        }

        check_samples(this_tree.root, ori,&this_tree);
        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        fix_condensed_nodes(&this_tree);
        this_tree.reassign_level();
    }
}
