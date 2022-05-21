#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#include <string>
namespace MAT = Mutation_Annotated_Tree;
#include "src/matOptimize/check_samples.hpp"
void Sample_Input(const char *name, Original_State_t& ori_state,
                  MAT::Tree &tree);
int main (int argc,char** argv){
    bool ignore_missings=false;
    if (argc==5) {
        ignore_missings=true;
        fprintf(stderr, "Ignoring missing samples\n");
    }
    MAT::Tree new_tree;
    MAT::load_mutation_annotated_tree(argv[3],new_tree);
    std::unordered_set<std::string> claimed_leaves;
    for (const auto leaf:new_tree.get_leaves()){
	auto name=new_tree.get_node_name(leaf->node_id);
    	auto ins_result=claimed_leaves.insert(name);
	if(!ins_result.second){
		fprintf(stderr,"%s duplicated\n", name.c_str());
	}
    }
    new_tree.uncondense_leaves();
    claimed_leaves.clear();
    for (const auto leaf:new_tree.get_leaves()){
	auto name=new_tree.get_node_name(leaf->node_id);
    	auto ins_result=claimed_leaves.insert(name);
	if(!ins_result.second){
		fprintf(stderr,"%s duplicated\n", name.c_str());
	}
    }
    MAT::Tree ori_tree;
    MAT::load_mutation_annotated_tree(argv[1],ori_tree);
    ori_tree.uncondense_leaves();
    Original_State_t ori_state;
    check_samples(ori_tree.root, ori_state, &ori_tree);
    Sample_Input(argv[2],ori_state,ori_tree);
    Original_State_t new_state;
    for (const auto& samp_pair : ori_state) {
        auto new_idx=new_tree.get_node(ori_tree.get_node_name(samp_pair.first));
        if (new_idx) {
            new_state.emplace(new_idx->node_id,std::move(samp_pair.second));            
        }else if (!ignore_missings) {
            fprintf(stderr, "Sample %s missing\n",ori_tree.get_node_name(samp_pair.first).c_str());
        }
    }
    check_samples(new_tree.root, new_state, &new_tree,ignore_missings);
}
