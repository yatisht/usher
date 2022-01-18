#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cassert>
#include <cstdint>
#include "src/matOptimize/check_samples.hpp"
#include <vector>
#define LEVEL_T uint8_t
#define MAX_LEVLEL UINT8_MAX
#define IDX_TREE_IDX_T uint32_t
#define EMPTY_POS UINT32_MAX
namespace MAT = Mutation_Annotated_Tree;
struct To_Place_Sample_Mutation{
    int position;
    uint8_t chrom_idx;
    uint8_t mut_nuc;
    union{
    struct{
        uint8_t par_nuc;
        uint8_t descendent_possible_nuc;
    };
    uint16_t range;
    };
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc), range(0) {
        assert(mut_nuc == 0xf);
    }
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc,
                             uint8_t par_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc),par_nuc(par_nuc)
        ,descendent_possible_nuc(0xf) {
        assert(mut_nuc != 0xf);
    }
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc,
                             uint8_t par_nuc,uint8_t descendant_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc),par_nuc(par_nuc)
        ,descendent_possible_nuc(descendant_nuc) {
        assert(mut_nuc != 0xf);
    }
    int get_end_range() const{
        if (mut_nuc==0xf) {
            return position+range;
        }else {
            return position;
        }
    }
};
struct Sample_Muts{
    std::string sample_name;
    std::vector<MAT::Mutation> muts;
};
void Sample_Input(const char *name, std::vector<Sample_Muts> &sample_mutations,
                  MAT::Tree &tree);
#ifndef NDEBUG
Mutation_Set get_mutations(const MAT::Node *main_tree_node);
void check_descendant_nuc(const MAT::Node* node);
#endif

void place_sample(Sample_Muts &&sample_to_place, MAT::Tree &main_tree
#ifndef NDEBUG
                  ,
                  Original_State_t &ori_state
#endif
);
void fix_parent(Mutation_Annotated_Tree::Node *root);
void convert_mut_type(std::vector<MAT::Mutation> &in,
                      std::vector<To_Place_Sample_Mutation> &out);
void assign_descendant_muts(MAT::Tree &in);
void assign_levels(MAT::Node* root);