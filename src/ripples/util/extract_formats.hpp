#ifndef EXTRACT_FORMATS_H
#define EXTRACT_FORMATS_H

#include "src/mutation_annotated_tree.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

void get_trios(MAT::Tree T, std::string filepath);

void get_parents(Mutation_Annotated_Tree::Tree *T,
                 std::unordered_set<std::string> &need_parents,
                 std::unordered_set<std::string> &all_nodes);

void generate_sample_paths(MAT::Tree &T);
std::vector<std::string> mutation_paths_no_label(MAT::Tree* T, std::vector<std::string> samples);
void leaves_per_node(MAT::Tree &T);

#endif

