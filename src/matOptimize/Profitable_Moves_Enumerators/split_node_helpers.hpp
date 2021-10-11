#include "Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
template<typename T>
static inline void add_node_split(const T& pos,nuc_one_hot ori_major_allele,nuc_one_hot allele1,nuc_one_hot allele2, Mutation_Count_Change_Collection& out,int& par_score) {
    auto major_allele=allele1&allele2;
    if (!major_allele) {
        par_score++;
        major_allele=allele1|allele2;
    }
    auto incremented=major_allele&(~ori_major_allele);
    auto decremented=ori_major_allele&(~major_allele);
    if (ori_major_allele!=major_allele) {
        out.emplace_back(pos,decremented,incremented);
    }
}
//at current node, but not in added node
static inline void add_node_split(const MAT::Mutation& mut,Mutation_Count_Change_Collection& out,int& par_score) {
    return add_node_split(mut,mut.get_all_major_allele(),mut.get_all_major_allele(),mut.get_par_one_hot(),out,par_score);
}
//at added node, but not at current node
static inline void add_node_split(const Mutation_Count_Change& mut,Mutation_Count_Change_Collection& out,int& par_score) {
    return add_node_split(mut,mut.get_par_state(),mut.get_par_state(),mut.get_incremented(),out,par_score);
}
//at both
static inline void add_node_split(const MAT::Mutation& mut, nuc_one_hot allele1, nuc_one_hot allele2,Mutation_Count_Change_Collection& out,int& par_score) {
    return add_node_split(mut,mut.get_all_major_allele(),allele2,allele1,out,par_score);
}