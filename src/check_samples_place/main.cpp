#include "src/matOptimize/mutation_annotated_tree.hpp"
namespace MAT = Mutation_Annotated_Tree;
#include "src/matOptimize/check_samples.hpp"
void Sample_Input(const char *name, Original_State_t& ori_state,
                  MAT::Tree &tree);
int main (int argc,char** argv){
    MAT::Tree ori_tree=MAT::load_mutation_annotated_tree(argv[1]);
    ori_tree.uncondense_leaves();
    Original_State_t ori_state;
    check_samples(ori_tree.root, ori_state, &ori_tree);
    Sample_Input(argv[2],ori_state,ori_tree);
    MAT::Tree new_tree=MAT::load_mutation_annotated_tree(argv[3]);
    new_tree.uncondense_leaves();
    check_samples(new_tree.root, ori_state, &new_tree);
}