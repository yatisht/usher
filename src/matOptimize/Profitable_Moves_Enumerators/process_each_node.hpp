#include "Profitable_Moves_Enumerators.hpp"
//Tag for merges that don't care about data unique to arg2
struct ignore_T2 {};
//Tag for merges that care about data unique to arg2
struct use_T2 {};
//Template for merge
template <typename Functor> class merge_func {
    typedef typename Functor::T1 T1;
    typedef typename Functor::T2 T2;
    //Process arg2 elements after arg1 elements are exhausted
    void process_rest_of_t2(range<typename T2::value_type> &T2_iter, Functor& functor,
                            ignore_T2 tag) {}
    void process_rest_of_t2(range<typename T2::value_type> &T2_iter, Functor& functor, use_T2 tag) {
        while (T2_iter) {
            functor.T2_only(*T2_iter);
            T2_iter++;
        }
    }

    void process_arg1(const typename T1::value_type& ele,range<typename T2::value_type>& T2_iter, Functor &functor) {
        //consume arg2 element smaller than arg1
        while (T2_iter && T2_iter->get_position() < ele.get_position()) {
            functor.T2_only(*T2_iter);
            T2_iter++;
        }
        //equal
        if (T2_iter && T2_iter->get_position() == ele.get_position()) {
            functor.T1_T2_match(ele, *T2_iter);
            T2_iter++;
        } else {
            //arg2>arg1
            functor.T1_only(ele);
        }
    }
  public:
    void operator()(const T1 &arg1, const T2 &arg2, Functor &functor) {
        range<typename T2::value_type> T2_iter(arg2);
        for (const auto &ele : arg1) {
            process_arg1(ele, T2_iter, functor);
        }
        process_rest_of_t2(T2_iter, functor, typename Functor::T2_useful());
    }
    void operator()(range<typename T1::value_type> arg1, const T2 &arg2, Functor &functor) {
        range<typename T2::value_type> T2_iter(arg2);
        for (; arg1; arg1++) {
            process_arg1(*arg1, T2_iter, functor);
        }
        process_rest_of_t2(T2_iter, functor, typename Functor::T2_useful());
    }
};

/**
 * @brief Caculate the fitch set change of nodes that is neither dst, parent of src, or LCA
 * @param[in] node the intermediate node to be caculated
 * @param[in] changed_child the children of node whose fitch set have changed
 * @param[in] this_node_mutation_count_change changes to fitch set of that child
 * @param[out] parent_node_mutation_count_change changes to Fitch set of this node
 * @param[out] parent_parsimony_score_change apply change in parsimony score due to
    change in number of mutations on edge between node and its direct children
 */
void get_intermediate_nodes_mutations(
    const MAT::Node *node,const MAT::Node* changed_child,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change);
/**
 * @brief Calculate the fitch set difference between fitch set of
    the new node formed by splitting the branch between node and that of parent of node to insert src
    (always assumed branch splitting in caculate parsimony score change) and parsimony score
    increase due to mutations unique to src
 * @param[in] node dst node
 * @param[in] children_added_mutations mutations needed to maintain the state of src when move it as sibling of dst
 * @param[out] parent_added_mutations The fitch set difference output, equivalent to the fitch set change of the child node
    on the branch dst node is as seen by its parent node
 * @param[out] parsimony_score_change parsimony score change due to this split
 * @return whether there are mutations unique to node compare to mutations in children_added_mutations, if nono, we
    can insert at children of dst as well
 */
bool get_parsimony_score_change_from_add(
    MAT::Node *node,
    const range<Mutation_Count_Change> &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations,
    int &parsimony_score_change);
//Similar to the functions above, but for effect of removing src node on fitch set of its parent,
// and the parsimony score change among children of parent of src
void get_parent_altered_remove(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    const MAT::Node *src, int &parent_parsimony_score_change);

static int get_new_major_allele_binary_node(nuc_one_hot left_child,
        nuc_one_hot right_child,
        nuc_one_hot &major_allele_out) {
    major_allele_out = left_child & right_child;
    if (!major_allele_out) {
        major_allele_out = left_child | right_child;
        return 1;
    }
    return 0;
}
//Helper function for adding Fitch set change from original state of a loci and changed state of a loci
static void register_change_from_new_state(Mutation_Count_Change_Collection &out,
        int new_mut_count, const MAT::Mutation &pos,
        nuc_one_hot new_state) {
    if (new_state != pos.get_all_major_allele()) {
        nuc_one_hot incremented_allele =
            new_state & (~pos.get_all_major_allele());
        nuc_one_hot decremented_allele =
            pos.get_all_major_allele() & (~new_state);
        out.emplace_back(pos,decremented_allele, incremented_allele);
    }
}


/**
 * @brief Caculate the fitch set change of LCA, if dst is not LCA
 * @param LCA the LCA node
 * @param src_branch Direct child of LCA that is on the path from src to LCA
 * @param is_src_terminal whether the src_branch_node is src node itself,
    if true, src_branch_node is removed fromm children of LCA
 * @param from_src_remove Fitch set changes of src_branch_node, or don't care if is_src_terminal is true
 * @param from_dst_add Fitch set changes of dst_branch_node, the child of LCA that is on path from dst to LCA
 * @param LCA_parent_mutation_count_change_out Fitch set change of LCA node
 * @param parsimony_score_change parsimony score change on the edges from LCA to its direct children
 * @return (void)
 */
void get_LCA_mutation(
    const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change);
void check_parsimony_score_change_above_LCA(MAT::Node *LCA, int &parsimony_score_change,
        Mutation_Count_Change_Collection &parent_added,
        const std::vector<MAT::Node *> &node_stack_from_src,
        std::vector<MAT::Node *> &node_stack_above_LCA,
        Mutation_Count_Change_Collection &parent_of_parent_added,
        MAT::Node *ancestor);