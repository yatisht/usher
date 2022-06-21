#ifndef INIT_PIPELINE_HPP
#define INIT_PIPELINE_HPP
#include <src/mutation_annotated_tree.hpp>
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

void preorder_traversal(MAT::Tree &T);
int find_long_branches(MAT::Tree &T, uint32_t branch_len, uint32_t num_descendants);

#endif
