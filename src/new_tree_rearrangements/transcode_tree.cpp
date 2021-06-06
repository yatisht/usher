#include "mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstdio>
int main (int argc,char** argv){
    Mutation_Annotated_Tree::Tree tree;
    tree.load_detatiled_mutations(argv[1]);
    printf("condensed: %zu\n",tree.condensed_nodes.size());
    Original_State_t ori;
    save_final_tree(tree,ori,argv[2]);
}