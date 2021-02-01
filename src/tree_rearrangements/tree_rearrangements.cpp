#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>
#include <tbb/pipeline.h>
#include "check_samples.hpp"
#include "conflicts.hpp"
tbb::concurrent_vector<MAT::Node*> postponed;

static void finalize_tree(Pending_Moves_t& pending_moves){
    tbb::parallel_for_each(pending_moves.begin(),pending_moves.end(),[](std::pair<MAT::Node*,ConfirmedMove>){
        
    });
}

void refine_trees(std::vector<MAT::Tree> &optimal_trees,int radius) {

    for (auto this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Sample_Mut_Type ori;
        std::vector<MAT::Node*>& new_nodes=this_tree.new_nodes;
        check_samples(this_tree.root, ori);
        Pending_Moves_t pending_moves;
        tbb::concurrent_vector<Profitable_Moves*> multiple_optimal_deferred;
        while (!new_nodes.empty()) {
            tbb::parallel_pipeline(10,
                tbb::make_filter<void,MAT::Node*>(tbb::filter::serial_in_order,Movable_Node_Enumerator(new_nodes,dfs_ordered_nodes))&
                tbb::make_filter<MAT::Node*,Possible_Moves*>(tbb::filter::parallel,Neighbors_Finder(radius))&
                tbb::make_filter<Possible_Moves*,Possible_Moves*>(tbb::filter::parallel,Parsimony_Score_Calculator({dfs_ordered_nodes}))&
                tbb::make_filter<Possible_Moves*,Profitable_Moves*>(tbb::filter::parallel,Profitable_Moves_Enumerator({dfs_ordered_nodes}))&
                tbb::make_filter<Profitable_Moves*,void>(tbb::filter::parallel,Move_Executor({pending_moves,multiple_optimal_deferred,dfs_ordered_nodes,this_tree}))
            );
            if(!multiple_optimal_deferred.empty()){
                tbb::concurrent_vector<Profitable_Moves*> ignored;
                tbb::parallel_pipeline(10,
                    tbb::make_filter<void,Profitable_Moves*> (tbb::filter::serial_in_order,multiple_moves_resolver(multiple_optimal_deferred))&
                    tbb::make_filter<Profitable_Moves*,void>(tbb::filter::parallel,Move_Executor({pending_moves,ignored,dfs_ordered_nodes,this_tree})));
                assert(ignored.empty());
            }
            if (!pending_moves.empty()) {
                //new_nodes=finalize_moves(moves);
                Sample_Mut_Type copy(ori);
                check_samples(this_tree.root, copy);
                dfs_ordered_nodes = this_tree.depth_first_expansion();
            }
            new_nodes.clear();
            new_nodes.insert(new_nodes.begin(),postponed.begin(),postponed.end());
        }

        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        this_tree.reassign_level();
    }
}
