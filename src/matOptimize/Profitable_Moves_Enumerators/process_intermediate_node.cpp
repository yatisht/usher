#include "process_each_node.hpp"
#include "process_individual_mutation.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cstdio>
#include <vector>
//Functor for general intermediate node
struct intermediate_mut_functor {
    typedef Mutation_Count_Change_Collection T1;
    typedef MAT::Mutations_Collection T2;
    typedef ignore_T2 T2_useful;
    Mutation_Count_Change_Collection &parent_node_mutation_count_change;
    int &parent_parsimony_score_change;
    intermediate_mut_functor(
        Mutation_Count_Change_Collection &parent_node_mutation_count_change,
        int &parent_parsimony_score_change)
        : parent_node_mutation_count_change(parent_node_mutation_count_change),
          parent_parsimony_score_change(parent_parsimony_score_change) {}
    // There is a change in state count (by at most 1) at its children, but the current node is not
    //sensitive at that loci (no allele have allele count difference less than 1 with major allele count),
    //so the state won't change, and increment parsimony score if it changed to a non-major allele.
    void T1_only(const T1::value_type &this_mut) {
        parent_parsimony_score_change += this_mut.get_default_change_internal();
    }
    void T2_only(const T2::value_type &) {}
    void T1_T2_match(const T1::value_type &this_mut,
                     const T2::value_type &LCA_mut) {
        //Match with a sensitive allele, use generic count change, as change_in can be in either direction
        decrement_increment_mutation_count(LCA_mut, this_mut,
                                           parent_node_mutation_count_change,
                                           parent_parsimony_score_change);
    }
};
struct intermediate_mut_functor_one_children {
    typedef Mutation_Count_Change_Collection T1;
    typedef MAT::Mutations_Collection T2;
    typedef ignore_T2 T2_useful;
    Mutation_Count_Change_Collection &parent_node_mutation_count_change;
    int &parent_parsimony_score_change;
    intermediate_mut_functor_one_children(
        Mutation_Count_Change_Collection &parent_node_mutation_count_change,
        int &parent_parsimony_score_change)
        : parent_node_mutation_count_change(parent_node_mutation_count_change),
          parent_parsimony_score_change(parent_parsimony_score_change) {}
    void T1_only(const T1::value_type &this_mut) {
        //specialization for nodes with one children, add the boundary one alleles (which is just all
        // the alleles other than major allele), they are not recorded in the tree explicitly, otherwise for leaves,
        //it is sensitive to all mutations, which unreasonably increase memory usage.
        MAT::Mutation this_node_mutation(this_mut.get_position());
        nuc_one_hot parent_state = this_mut.get_par_state();
        this_node_mutation.set_par_mut(parent_state, parent_state);
        this_node_mutation.set_auxillary(parent_state, (~parent_state) & 0xf);
        decrement_increment_mutation_count(this_node_mutation, this_mut,
                                           parent_node_mutation_count_change,
                                           parent_parsimony_score_change);
    }
    void T2_only(const T2::value_type &) {}
    void T1_T2_match(const T1::value_type &this_mut,
                     const T2::value_type &LCA_mut) {
        decrement_increment_mutation_count(LCA_mut, this_mut,
                                           parent_node_mutation_count_change,
                                           parent_parsimony_score_change);
    }
};
//Specialization for nodes with two children, they are most common kind of nodes

void get_intermediate_nodes_mutations(
    const MAT::Node *node,const MAT::Node* changed_child,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change) {
    //Dispatcher to use different specialization based on node count
    if (node->children.size() == 1) {
        intermediate_mut_functor_one_children functor(
            parent_node_mutation_count_change, parent_parsimony_score_change);
        merge_func<intermediate_mut_functor_one_children>()(
            this_node_mutation_count_change, node->mutations, functor);
    }else {
        intermediate_mut_functor functor(parent_node_mutation_count_change,
                                         parent_parsimony_score_change);
        merge_func<intermediate_mut_functor>()(this_node_mutation_count_change,
                                               node->mutations, functor);
    }
}