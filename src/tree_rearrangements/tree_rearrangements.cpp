#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include "check_samples.hpp"

tbb::concurrent_hash_map<Mutation_Annotated_Tree::Node*, Mutation_Annotated_Tree::Node*> potential_crosses;

tbb::concurrent_hash_map<Mutation_Annotated_Tree::Node*, tbb::concurrent_vector<int>> repeatedly_mutating_loci;

tbb::concurrent_vector<MAT::Node> postponed;


void refine_trees(std::vector<MAT::Tree> &optimal_trees,int radius) {

    for (auto this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Sample_Mut_Type ori;
        std::vector<MAT::Node*>& new_nodes=this_tree.new_nodes;
        check_samples(this_tree.root, ori);
        std::vector<ConfirmedMove> moves;
        while (!new_nodes.empty()) {
            tbb::interface6::parallel_pipeline(10,
                tbb::make_filter<void,MAT::Node*>(tbb::filter::serial_in_order,Movable_Node_Enumerator(new_nodes,dfs_ordered_nodes))&
                tbb::make_filter<MAT::Node*,Possible_Moves*>(tbb::filter::parallel,Neighbors_Finder(radius))&
                tbb::make_filter<Possible_Moves*,Possible_Moves*>(tbb::filter::parallel,Parsimony_Score_Calculator())&
                tbb::make_filter<Possible_Moves*,Profitable_Moves*>(tbb::filter::parallel,Profitable_Moves_Enumerator())&
                tbb::make_filter<Profitable_Moves*,void>(tbb::filter::serial_in_order,Move_Executor(moves))
            );
            if (!moves.empty()) {
                new_nodes=finalize_moves(moves);
                Sample_Mut_Type copy(ori);
                check_samples(this_tree.root, copy);
                dfs_ordered_nodes = this_tree.depth_first_expansion();
            }
        }

        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        this_tree.reassign_level();
    }
}
