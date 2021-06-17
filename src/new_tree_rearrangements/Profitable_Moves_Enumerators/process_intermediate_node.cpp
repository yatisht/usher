#include "process_each_node.hpp"
#include "process_individual_mutation.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
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
int get_new_major_allele_binary_node(nuc_one_hot left_child,
                                     nuc_one_hot right_child,
                                     nuc_one_hot &major_allele_out) {
    major_allele_out = left_child & right_child;
    if (!major_allele_out) {
        major_allele_out = left_child | right_child;
        return 1;
    }
    return 0;
}

int register_change_from_new_state(Mutation_Count_Change_Collection &out,
                                   int new_mut_count, const MAT::Mutation &pos,
                                   nuc_one_hot new_state) {
    if (new_state != pos.get_all_major_allele()) {
        nuc_one_hot incremented_allele =
            new_state & (~pos.get_all_major_allele());
        nuc_one_hot decremented_allele =
            pos.get_all_major_allele() & (~new_state);
        out.emplace_back(pos,decremented_allele, incremented_allele,
                              new_state);
    }
    int old_mutation_count =
        !(pos.get_left_child_state() & pos.get_right_child_state());
    return new_mut_count - old_mutation_count;
}
//Specialization for nodes with two children, they are most common kind of nodes
struct get_two_child_intermediate_node_mutations {
    Mutation_Count_Change_Collection &parent_node_mutation_count_change;
    int &parent_parsimony_score_change;
    bool changed_first_children;
    typedef Mutation_Count_Change_Collection T1;
    typedef MAT::Mutations_Collection T2;
    typedef ignore_T2 T2_useful;
    get_two_child_intermediate_node_mutations(
        Mutation_Count_Change_Collection &parent_node_mutation_count_change,
        int &parent_parsimony_score_change, const MAT::Node *this_node,
        const MAT::Node *changed_child)
        : parent_node_mutation_count_change(parent_node_mutation_count_change),
          parent_parsimony_score_change(parent_parsimony_score_change) {
        if (changed_child == this_node->children[0]) {
            changed_first_children = true;
        } else {
            changed_first_children = false;
            assert(changed_child == this_node->children[1]);
        }
    }
    void T1_only(const T1::value_type &this_mut) {
        parent_parsimony_score_change += this_mut.get_default_change_internal();
    }
    void T2_only(const T2::value_type &) {}
    void T1_T2_match(const T1::value_type &this_mut,
                     const T2::value_type &LCA_mut) {
        //Basically substitute the changed node and recalculate majjor allele from scratch
        nuc_one_hot major_allele,left_child_state,right_child_state;
        if (changed_first_children) {
            left_child_state = this_mut.get_new_state();
            right_child_state = LCA_mut.get_right_child_state();
        } else {
            left_child_state = LCA_mut.get_left_child_state();
            right_child_state = this_mut.get_new_state();
        }
        int new_count = get_new_major_allele_binary_node(
            left_child_state, right_child_state, major_allele);
        parent_parsimony_score_change +=
            register_change_from_new_state(parent_node_mutation_count_change,
                                           new_count, LCA_mut, major_allele);
    }
};

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
    } else if(node->children.size()==2){
        struct get_two_child_intermediate_node_mutations functor(
        parent_node_mutation_count_change, parent_parsimony_score_change,
        node, changed_child);
    merge_func<struct get_two_child_intermediate_node_mutations>()(
        this_node_mutation_count_change, node->mutations, functor);
    }else{
        intermediate_mut_functor functor(parent_node_mutation_count_change,
                                         parent_parsimony_score_change);
        merge_func<intermediate_mut_functor>()(this_node_mutation_count_change,
                                               node->mutations, functor);
    }
}