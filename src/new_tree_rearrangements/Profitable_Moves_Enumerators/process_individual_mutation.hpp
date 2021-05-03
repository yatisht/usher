#ifndef PROCESS_INDIVIDUAL_MUTATION
#define PROCESS_INDIVIDUAL_MUTATION

#include "Profitable_Moves_Enumerators.hpp"
nuc_one_hot decrement_mutation_count(Mutation_Count_Change_Collection &out,
                         const MAT::Mutation &parent_mut,const Mutation_Count_Change& change_in,
                         int& score_change);
nuc_one_hot increment_mutation_count(Mutation_Count_Change_Collection &out,
                         const MAT::Mutation &parent_mut,
                         const Mutation_Count_Change& change_in,int& score_change) ;

nuc_one_hot decrement_increment_mutation_count(
    const MAT::Mutation &parent_mutation, Mutation_Count_Change change_in,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &score_change);


nuc_one_hot dbl_inc_dec_mutations(
    const Mutation_Count_Change &dst_add_count,
    const Mutation_Count_Change &src_count_change,
    const MAT::Mutation &this_mut, int &score_change,
    Mutation_Count_Change_Collection &parent_mutation_change_out);

int get_new_major_allele_binary_node(nuc_one_hot left_child,nuc_one_hot right_child,nuc_one_hot& major_allele_out);
#endif