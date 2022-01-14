#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cstdint>
#include "src/matOptimize/check_samples.hpp"
#include <vector>
#define LEVEL_T uint8_t
#define MAX_LEVLEL UINT8_MAX
#define IDX_TREE_IDX_T uint32_t
#define EMPTY_POS UINT32_MAX
namespace MAT = Mutation_Annotated_Tree;
struct Sampled_Tree_Mutation{
    int position;
    uint8_t chrom_idx;
    uint8_t mut_nuc;
    uint8_t descendent_possible_nuc;
    uint8_t par_nuc;
    inline bool operator< (const Sampled_Tree_Mutation& m) const {
        return ((*this).position < m.position);
    }
};
struct Sampled_Tree_Node{
    Sampled_Tree_Node* parent;
    const MAT::Node* corresponding_main_tree_node;
    std::vector<Sampled_Tree_Mutation> mutations;
    std::vector<Sampled_Tree_Node*> children;
    int dfs_idx;
};
struct To_Place_Sample_Mutation{
    int position;
    uint8_t chrom_idx;
    uint8_t mut_nuc;
    union{
    uint8_t par_nuc;
    uint16_t range;
    };
};
Sampled_Tree_Node *sample_tree(MAT::Tree &in, int threshold);

struct Sample_Muts{
    std::string sample_name;
    std::vector<Sampled_Tree_Mutation> muts;
};
void Sample_Input(const char *name, std::vector<Sample_Muts> &sample_mutations,
                  MAT::Tree &tree);
#ifndef NDEBUG
Mutation_Set get_mutations(const MAT::Node *main_tree_node);
void check_sampled_tree(MAT::Tree &main_tree,
                    std::vector<Sampled_Tree_Node *>& sampled_tree_dfs,
                    int distance);
#endif

void place_sample(Sample_Muts &&sample_to_place,
                  Sampled_Tree_Node *sampled_tree_root, MAT::Tree &main_tree,
                  int sampling_radius
#ifndef NDEBUG
                  ,
                  Original_State_t &ori_state
#endif
);
void sample_tree_dfs(Sampled_Tree_Node *sampled_tree_root,std::vector<Sampled_Tree_Node *>& output);
void
set_parent_muts(std::vector<Sampled_Tree_Mutation> &mutations_to_set,
                const MAT::Node *node);
void fix_parent(Mutation_Annotated_Tree::Node *root);
void
check_sampled_main_correspondence(const Sampled_Tree_Node *sampled_tree_node);