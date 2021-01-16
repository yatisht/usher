#include "mutation_annotated_tree.hpp"
namespace Tree_Rearrangement{
    bool move_nearest(Mutation_Annotated_Tree::Node *this_node,
                  std::vector<Mutation_Annotated_Tree::Node *> &dfs_ordered_nodes,
                  Mutation_Annotated_Tree::Tree &tree);
}