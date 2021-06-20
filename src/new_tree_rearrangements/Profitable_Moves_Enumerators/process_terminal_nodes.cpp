#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "process_individual_mutation.hpp"
//functor for adding src to split the edge between dst and its parent, same as geting major allele of a binary node
struct get_parsimony_score_from_add {
    typedef Mutation_Count_Change_Collection T1;
    typedef MAT::Mutations_Collection T2;
    typedef use_T2 T2_useful;
    Mutation_Count_Change_Collection &parent_added_mutations;
    int &parsimony_score_change;
    bool &have_not_shared;
    get_parsimony_score_from_add(
        Mutation_Count_Change_Collection &parent_node_mutation_count_change,
        int &parent_parsimony_score_change, bool &have_not_shared)
        : parent_added_mutations(parent_node_mutation_count_change),
          parsimony_score_change(parent_parsimony_score_change),
          have_not_shared(have_not_shared) {}
    void T1_only(const T1::value_type &added_child_mutation) {
        //mutation unique to moved src, so the dst node have par_state
        nuc_one_hot major_allele;
        int new_score = get_new_major_allele_binary_node(
            added_child_mutation.get_par_state(),
            added_child_mutation.get_incremented(), major_allele);
        have_not_shared = have_not_shared || (new_score);
        parsimony_score_change += new_score;
        if (major_allele != added_child_mutation.get_par_state()) {
            parent_added_mutations.emplace_back(added_child_mutation,
                                                major_allele);
        }
    }
    void T2_only(const T2::value_type &sibling_addable) {
        // Adding a children with parent state
        nuc_one_hot major_alleles = sibling_addable.get_all_major_allele() &
                                    sibling_addable.get_par_one_hot();
        if (!major_alleles) {
            major_alleles = sibling_addable.get_all_major_allele() |
                            sibling_addable.get_par_one_hot();
            have_not_shared = true;
            parsimony_score_change++;
        }
        register_change_from_new_state(parent_added_mutations, 0,
                                       sibling_addable, major_alleles);
    }
    void T1_T2_match(const T1::value_type &added_child_mutation,
                     const T2::value_type &sibling_addable) {
        nuc_one_hot major_alleles;
        int new_score = get_new_major_allele_binary_node(
            sibling_addable.get_all_major_allele(),
            added_child_mutation.get_incremented(), major_alleles);
        have_not_shared = have_not_shared || (new_score);
        parsimony_score_change += new_score;
        register_change_from_new_state(parent_added_mutations, new_score,
                                       sibling_addable, major_alleles);
    }
};

bool get_parsimony_score_change_from_add(
    MAT::Node *node,
    const range<Mutation_Count_Change> &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations,
    int &parsimony_score_change) {
    bool have_not_shared = false;
    get_parsimony_score_from_add functor(parent_added_mutations,parsimony_score_change,have_not_shared);
    merge_func<get_parsimony_score_from_add>()(children_added_mutations,node->mutations,functor);
    return have_not_shared;
}
//Functor for removing src from a non-binary node
struct get_parent_altered_remove {
    typedef MAT::Mutations_Collection T1;
    typedef MAT::Mutations_Collection T2;
    typedef use_T2 T2_useful;
    Mutation_Count_Change_Collection &parent_mutation_count_change_out;
    int &parent_parsimony_score_change;
    get_parent_altered_remove(
        Mutation_Count_Change_Collection &parent_node_mutation_count_change,
        int &parent_parsimony_score_change)
        : parent_mutation_count_change_out(parent_node_mutation_count_change),
          parent_parsimony_score_change(parent_parsimony_score_change) {}
    void T1_only(const T1::value_type &src_mut) {
            //only src have this mutation, so decrement parsimony score if removing a valid mutation
            if (src_mut.is_valid()) {
                parent_parsimony_score_change--;
            }
    }
    void T2_only(const T2::value_type &parent_variable) {
            Mutation_Count_Change temp(parent_variable,parent_variable.get_mut_one_hot(),0,0,true);
            decrement_mutation_count(parent_mutation_count_change_out, parent_variable, temp,parent_parsimony_score_change);
            //decrement_mutation_count manipulate parsimony score change only on how such change change the number of children able to have major allele state (as if src branch node is not removed), so need to decrement artificailly to account for src->parent have one less children
            parent_parsimony_score_change--;
    }
    void T1_T2_match(const T1::value_type &src_mut,
                     const T2::value_type &parent_variable) {
            Mutation_Count_Change temp(src_mut,src_mut.get_all_major_allele(), 0, 0,true);
            decrement_mutation_count(
                parent_mutation_count_change_out, parent_variable, temp,
                parent_parsimony_score_change);
            //same as above
            parent_parsimony_score_change--;
    }
};
/*
parent of parent (not involved in this function, but the output
parent_mutation_count_change is used by it to adjust its mutation list)
|
parent (parent_mutation_count_change reports what happens to be mutation count
at each loci for this node)
|
src (being removed, amount to decrement count of all its mutations)
*/
void get_parent_altered_remove(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    const MAT::Node *src, int &parent_parsimony_score_change
) {
    struct get_parent_altered_remove functor(parent_mutation_count_change_out,parent_parsimony_score_change);
    merge_func<struct get_parent_altered_remove>()(src->mutations,src->parent->mutations,functor);
}

//specialization when src->parent is binary, then it just have the same major allele set as its remaining children 
void get_child_removed_binary_node(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    const MAT::Node *src, int &parent_parsimony_score_change) {
    MAT::Node *parent_node = src->parent;
    bool removed_left_child = parent_node->children[0] == src;
    assert(removed_left_child || parent_node->children[1] == src);
    for (const auto &mutation : parent_node->mutations) {
        nuc_one_hot new_major_allele = removed_left_child
                                           ? mutation.get_right_child_state()
                                           : mutation.get_left_child_state();
        int score_change = register_change_from_new_state(
            parent_mutation_count_change_out, 0, mutation, new_major_allele);
        parent_parsimony_score_change += score_change;
    }
}
