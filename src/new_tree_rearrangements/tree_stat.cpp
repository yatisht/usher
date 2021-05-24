#include "mutation_annotated_tree.hpp"
using namespace Mutation_Annotated_Tree;
using namespace std;
int main(int argc, char** argv){
    Tree t=load_mutation_annotated_tree(argv[1]);
    vector<Node*> nodes=t.breadth_first_expansion();
    auto src_node=t.get_node("2273");
    auto dst_node=t.get_node("2086");

}