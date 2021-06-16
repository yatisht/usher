#include "Profitable_Moves_Enumerators.hpp"
//Tag for merges that don't care about data unique to arg2 
struct ignore_T2 {};
//Tag for merges that care about data unique to arg2 
struct use_T2 {};
//Template for merge
template <typename Functor> class merge_func {
    typedef typename Functor::T1 T1;
    typedef typename Functor::T2 T2;
    void process_rest_of_t2(range<T2> &T2_iter, Functor functor,
                            ignore_T2 tag) {}
    void process_rest_of_t2(range<T2> &T2_iter, Functor functor, use_T2 tag) {
        while (T2_iter) {
            functor.T2_only(*T2_iter);
            T2_iter++;
        }
    }

  public:
    void operator()(const T1 &arg1, const T2 &arg2, Functor &functor) {
        range<T2> T2_iter(arg2);
        for (const auto &ele : arg1) {
            while (T2_iter && T2_iter->get_position() < ele.get_position()) {
                functor.T2_only(*T2_iter);
                T2_iter++;
            }
            if (T2_iter && T2_iter->get_position() == ele.get_position()) {
                functor.T1_T2_match(ele, *T2_iter);
                T2_iter++;
            } else {
                functor.T1_only(ele);
            }
        }
        process_rest_of_t2(T2_iter, functor, typename Functor::T2_useful());
    }
};

void get_intermediate_nodes_mutations(
    const MAT::Node *node,const MAT::Node* changed_child,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change);
bool get_parsimony_score_change_from_add(
    MAT::Node *parent_node,
    const Mutation_Count_Change_Collection &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations,
    int &parsimony_score_change);
void get_parent_altered_remove(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    const MAT::Node *src, int &parent_parsimony_score_change);

int register_change_from_new_state(Mutation_Count_Change_Collection &out,int new_mut_count,
                            const MAT::Mutation &pos,nuc_one_hot new_state);

int get_new_major_allele_binary_node(nuc_one_hot left_child,nuc_one_hot right_child,nuc_one_hot& major_allele_out);

void get_child_removed_binary_node(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    const MAT::Node *src, int &parent_parsimony_score_change);

void get_LCA_mutation(
    const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change);