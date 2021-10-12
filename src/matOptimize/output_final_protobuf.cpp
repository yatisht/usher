#include "mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
std::vector<std::string> changed_nodes;
int main(int argc, char** argv) {
    Mutation_Annotated_Tree::Tree tree;
    tree.load_detatiled_mutations(argv[1]);
    save_final_tree(tree, argv[2]);
}