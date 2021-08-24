#ifdef CHECK_STATE_REASSIGN
#include "apply_move.hpp"
#include "src/matOptimize/Fitch_Sankoff.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include <vector>
MAT::Tree reassign_state_full(MAT::Tree &tree_in) {
    MAT::Tree new_tree;
    MAT::Node *new_ancestor = new Mutation_Annotated_Tree::Node(
        *tree_in.root, nullptr, &new_tree, false);
    new_tree.curr_internal_node = tree_in.curr_internal_node;
    new_tree.root = new_ancestor;
    std::vector<MAT::Node *> bfs_ordered_nodes =
        new_tree.breadth_first_expansion();
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
            output(bfs_ordered_nodes.size());
    tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&bfs_ordered_nodes, &output,&tree_in](
            const std::pair<MAT::Mutation,
            std::unordered_map<std::string, nuc_one_hot> *>
    &pos) {
        std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
        Fitch_Sankoff_Whole_Tree(bfs_ordered_nodes, pos.first, *mutated,
                                 output,&tree_in);
    });
    tbb::affinity_partitioner ap;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
    [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            const auto &to_refill = output[i];
            bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                         to_refill.size());
        }
    },
    ap);

    return new_tree;
}
#endif