#include "common.hpp"

std::vector<std::string> mutation_paths(MAT::Tree* T, std::vector<std::string> samples);
std::vector<std::string> clade_paths(MAT::Tree* T);
std::vector<std::string> all_nodes_paths(MAT::Tree* T);
void print_node_stats(Mutation_Annotated_Tree::Node* node, size_t& leaf_count, size_t& mut_count, std::ofstream& outfile);
