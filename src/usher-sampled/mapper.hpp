#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <vector>
struct Sampled_Place_Target {
    const Sampled_Tree_Node *target_node;
    std::vector<To_Place_Sample_Mutation> muts;
    Sampled_Place_Target() = default;
    Sampled_Place_Target(const Sampled_Tree_Node *target_node,
                         std::vector<To_Place_Sample_Mutation> &muts)
        : target_node(target_node), muts(muts) {}
};
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
std::vector<Sampled_Place_Target>
place_on_sampled_tree(Sampled_Tree_Node *sampled_tree_root,
                      std::vector<To_Place_Sample_Mutation> &&sample_mutations,
                      int& parsimony_score
#ifndef NDEBUG
                      ,
                      Mutation_Set &sample_mutations_set
#endif
) ;
std::tuple<Main_Tree_Target, int,int>
place_main_tree(std::vector<Sampled_Place_Target> &sampled_output,
                MAT::Tree &main_tree, int sampling_radius
#ifndef NDEBUG
                ,
                Mutation_Set &sample_mutations
#endif
) ;
#ifndef NDEBUG
void check_mutations(Mutation_Set ref,const Main_Tree_Target& target_to_check);
void optimality_check(Mutation_Set &sample_mutations, int parsimony,
                      MAT::Tree& main_tree, int sampling_radius,
                      Sampled_Tree_Node *sample_tree_root,
                      const std::vector<Sampled_Place_Target> &sampled_out,
                      std::vector<To_Place_Sample_Mutation>& samples_to_place);
void check_sampled_mutations(Mutation_Set ref,const Sampled_Place_Target& target_to_check);
void check_continuation(const MAT::Node* parent_node,Mutation_Set ref,const std::vector<To_Place_Sample_Mutation> &decendent_mutations);
#endif
void
set_parent_muts(std::vector<To_Place_Sample_Mutation> &mutations_to_set,const MAT::Node *node);