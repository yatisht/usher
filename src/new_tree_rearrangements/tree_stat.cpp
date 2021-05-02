#include "mutation_annotated_tree.hpp"
using namespace Mutation_Annotated_Tree;
using namespace std;
int main(int argc, char** argv){
    Tree t=load_mutation_annotated_tree(argv[1]);
    vector<Node*> nodes=t.breadth_first_expansion();
    for(auto n:nodes){
        //if (n->children.size()>1) {
            printf("%zu\n",n->mutations.size());
        //}
    }
}