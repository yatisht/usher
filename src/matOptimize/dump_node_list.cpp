#include "mutation_annotated_tree.hpp"
#include <cstdio>
int main(int argc, char** argv) {
    FILE* fh=fopen("nodelist", "w");
    Mutation_Annotated_Tree::Tree tree;
    tree.load_detatiled_mutations(argv[1]);
    for (const auto & node : tree.all_nodes) {
        fprintf(fh, "%s\n",node.first.c_str());
    }
    fclose(fh);
}