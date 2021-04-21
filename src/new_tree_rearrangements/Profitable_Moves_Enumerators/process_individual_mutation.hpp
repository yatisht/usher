#include "Profitable_Moves_Enumerators.hpp"
nuc_one_hot decrement_mutation_count(Mutation_Count_Change_Collection &out,
                         const MAT::Mutation &parent_mut,const Mutation_Count_Change& change_in,
                         int& score_change,bool is_terminal);
nuc_one_hot increment_mutation_count(Mutation_Count_Change_Collection &out,
                         const MAT::Mutation &parent_mut,
                         const Mutation_Count_Change& change_in,int& score_change,bool is_terminal) ;

nuc_one_hot decrement_increment_mutation_count(
    const MAT::Mutation &parent_mutation, Mutation_Count_Change change_in,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &score_change);


nuc_one_hot dbl_inc_dec_mutations(
    const Mutation_Count_Change &dst_add_count, bool is_dst_terminal,
    const Mutation_Count_Change &src_count_change, bool is_src_terminal,
    const MAT::Mutation &this_mut, int &score_change,
    Mutation_Count_Change_Collection &parent_mutation_change_out);
