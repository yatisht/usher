#include "usher.hpp"
#include <vector>
#pragma once
template <typename Target_Type> struct Output {
    std::mutex mutex;
    int best_par_score;
    std::vector<Target_Type> targets;
};
struct Main_Tree_Target {
    MAT::Node *target_node;
    MAT::Node *parent_node;
    MAT::Mutations_Collection splited_mutations;
    std::vector<To_Place_Sample_Mutation> sample_mutations;
    MAT::Mutations_Collection shared_mutations;
    int distance_left;
};
std::tuple<std::vector<Main_Tree_Target>, int>
place_main_tree(const std::vector<To_Place_Sample_Mutation> &mutations,
                MAT::Tree &main_tree
#ifdef DETAILED_MERGER_CHECK
                ,
                Mutation_Set &sample_mutations
#endif
               ) ;
#ifndef NDEBUG
void check_mutations(Mutation_Set ref,const Main_Tree_Target& target_to_check);
void check_continuation(const MAT::Node* parent_node,Mutation_Set ref,const std::vector<To_Place_Sample_Mutation> &decendent_mutations);
#endif
void
set_parent_muts(std::vector<To_Place_Sample_Mutation> &mutations_to_set,const MAT::Node *node);