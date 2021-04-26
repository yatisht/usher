#include "Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
static void register_change(Mutation_Count_Change_Collection &out,
                            const MAT::Mutation &pos,
                            nuc_one_hot decremented_allele,
                            nuc_one_hot incremented_allele,nuc_one_hot new_state) {
    out.emplace_back(pos);
    Mutation_Count_Change &last_added = out.back();
    last_added.set_change(decremented_allele, incremented_allele,new_state);
}
nuc_one_hot decrement_mutation_count(Mutation_Count_Change_Collection &out,
                                     const MAT::Mutation &parent_mut,
                                     const Mutation_Count_Change &change_in,
                                     int &score_change, bool is_terminal) {
    // Score_change assumes we are operating on internal nodes, so it is
    // removing a possibility of the state a child may take, not the whole child
    // altogether
    nuc_one_hot decremented_allele = change_in.get_decremented();
    nuc_one_hot all_nuc =
        parent_mut.get_all_major_allele();
    nuc_one_hot effectively_decremented_alleles = all_nuc & decremented_allele;
    if (effectively_decremented_alleles) {
        // Decrementing tied alleles, inform parent of parent an allele is no
        // longer tied in parent, thus decrement it count
        nuc_one_hot not_decremented =  all_nuc & (~decremented_allele);
        nuc_one_hot boundary1_not_decremented =
            parent_mut.get_boundary1_one_hot() & (~decremented_allele);
        if (not_decremented) {
            // Some originally major alleles is not decremented, so they become
            // new major alleles (no change in parent of parent allele count
            // needed)
            // Decremented major alleles are no longer tie, need to be removed
            // from parent of parent alleles count
            register_change(out, parent_mut, effectively_decremented_alleles,
                            0,not_decremented);
            // Parsimony score change: If the allele is originally tied with
            // another allele, then they have same number of children that can
            // only be that state, the removed mutation is not one of those,
            // otherwise it is removing a child altogether. Then the parsimony
            // score remain the number of children whose state can only be that
            // of the other major allele, which stay the same
            
            if(is_terminal){
                score_change--;
            }
            return not_decremented;
        } else {
            // One additional child can only have state of minor alleles, if it
            // is not complete removal
            if (!is_terminal) {
                score_change++;
            }

            // all alleles originally tied got decremented
            if (boundary1_not_decremented) {
                // Then everything become tie with the old boundary1 allele that
                // aren't decremented
                nuc_one_hot major_allele=boundary1_not_decremented|all_nuc;
                register_change(out, parent_mut, 0, boundary1_not_decremented,major_allele);
                return major_allele;
            }
            // everything got decremented, major alleles stay the same, boundary
            // 1 tie with boundary2, but no impact on parent of parent
            return all_nuc;
        }
    }
    // Decrementing minor alleles that is an alternative have no impact if it
    // doesn't remove the child
    if (is_terminal&&(!(change_in.get_par_state()&change_in.get_decremented()))) {
        score_change--;
    }
    return all_nuc;
}

nuc_one_hot increment_mutation_count(Mutation_Count_Change_Collection &out,
                                     const MAT::Mutation &parent_mut,
                                     const Mutation_Count_Change &change_in,
                                     int &score_change, bool is_terminal) {
    // Break ties as before
    nuc_one_hot incremented_allele = change_in.get_incremented();
    nuc_one_hot all_nuc =
        parent_mut.get_all_major_allele();
    nuc_one_hot effectively_incremented_alleles = all_nuc & incremented_allele;
    if (is_terminal&&(!(effectively_incremented_alleles))) {
        score_change++;
    }
    if (effectively_incremented_alleles) {
        // Not incremented alleles no longer in the tie
        nuc_one_hot not_incremented = effectively_incremented_alleles ^ all_nuc;
        if ((!is_terminal)) {
            // originally a mutation to some other allele, now can follow major
            // (potentially new )allele
            score_change--;
        }
        if (not_incremented) {
            // Incremented major alleles are still major alleles
            // Not incremented major alleles was originally counted in count of
            // children mutations by parent of parent, they do not change the
            // number of children whose state that can only be minor alleles
            register_change(out, parent_mut, not_incremented, 0,effectively_incremented_alleles);
        }
        // it just increment everything, no change in tie

        // boundary alleles will never tie with major allele, if any major/tied
        // alleles are incremented
        return effectively_incremented_alleles;
    } else if (parent_mut.get_boundary1_one_hot() & incremented_allele) {
        // If none of the originally tieing major alleles are incremented,
        // boundary 1 alleles may join ties
        nuc_one_hot tieing_alleles =
            parent_mut.get_all_major_allele() |
            (parent_mut.get_boundary1_one_hot() & incremented_allele);
        register_change(out, parent_mut, 0,
                        parent_mut.get_boundary1_one_hot() &
                            incremented_allele,tieing_alleles);
        // Either the parent is original allele, or the new major allele, number
        // of children who cant follow parent state is the same, and may
        // increment if this is a new child
        return tieing_alleles;
    }
    // Adding a minor allele, if adding new child, increment
    return parent_mut.get_all_major_allele();
}

nuc_one_hot decrement_increment_mutation_count(
    const MAT::Mutation &parent_mutation, Mutation_Count_Change change_in,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &score_change) {
    nuc_one_hot incremented = change_in.get_incremented();
    nuc_one_hot decremented = change_in.get_decremented();
    assert(!(decremented & incremented));

    nuc_one_hot all_major_alleles = parent_mutation.get_all_major_allele();
    nuc_one_hot major_alleles = all_major_alleles;

    nuc_one_hot boundary1_state = parent_mutation.get_boundary1_one_hot();

    // handle decrementing major allele, incrementing minor allele
    nuc_one_hot major_allele_decremented = all_major_alleles & decremented;
    nuc_one_hot major_allele_incremented = all_major_alleles & incremented;
    nuc_one_hot major_allele_not_decremented =
        all_major_alleles & (~decremented);

    if (major_allele_incremented) {
        // Regardless of what happens to minor alleles or other major alleles,
        // the incremented major alleles will be the new major allele
        // (incremented minor alleles is at least 1 less than incremented major
        // allele) It is also the same in reporting change to parent of parent,
        // whether a major allele is decremented directly or not, it will report
        // a decrement if it is not incremented, as they are no longer in tie
        // with incremented major allele either way.
        major_alleles = increment_mutation_count(
            parent_node_mutation_count_change, parent_mutation, change_in,
            score_change, false);
        change_in.set_change(0,major_allele_incremented,major_alleles);
        return major_alleles;
    } else if (major_allele_decremented) {
        // this is adding to a minor allele while removing a major
        // allele
        if (boundary1_state & incremented || major_allele_not_decremented) {
            assert(!major_allele_incremented);
            // tied alleles not decremented (no action needed)
            // and boundary 1 allele incremented (need to be incremented in
            // parent of parent counts) become new major allele
            // Tied major alleles decremented joins minor alleles (need to be
            // decremented in parent of parent counts)
            major_alleles =
                major_allele_not_decremented | (boundary1_state & incremented);
            register_change(parent_node_mutation_count_change, parent_mutation,
                major_allele_decremented,
                boundary1_state & incremented,major_alleles);
            // The number of children that can follow any of the major alleles
            // didn't change, as there are still same numbe of children that can
            // follow major_alleles
            return major_alleles;
        } else if (parent_mutation.get_boundary_by_2() & incremented) {
            assert(!(major_allele_not_decremented));
            assert(!(boundary1_state & incremented));
            // All major allele are decremented to tie with boundary 1 alleles,
            // no boundary1 alleles are incremented,some boundary2 alleles get
            // incremented. They become the new major alleles
            // so no change in mutations (it just follow parent),just change in
            // children mutation count in parent of parent;
            nuc_one_hot new_major_alleles =
                (boundary1_state & (!decremented)) |
                (parent_mutation.get_boundary_by_2() & incremented);

            major_alleles = new_major_alleles | all_major_alleles;

            register_change(parent_node_mutation_count_change, parent_mutation,
                            0, new_major_alleles,major_alleles);
            // One less node that follow parent
            score_change++;
            return major_alleles;
        } else {
            // No thing interesting have been incremented
            major_alleles=decrement_mutation_count(parent_node_mutation_count_change,
                                            parent_mutation, change_in,
                                            score_change, false);
            change_in.set_change(decremented,0,major_alleles);
            return major_alleles;
        }
    } else {
        // Nothing happens to major alleles, for minor alleles, only increment
        // is interesting
        if (incremented) {
            return increment_mutation_count(parent_node_mutation_count_change,
                                            parent_mutation, change_in,
                                            score_change, false);
        }
    }
    return major_alleles;
}

static void score_change_internal_node(nuc_one_hot original_state,const Mutation_Count_Change & change, int& score_change ){
    if (original_state&change.get_incremented()) {
        score_change--;
    }
    if(original_state&change.get_decremented()){
        score_change++;
    }
}

nuc_one_hot dbl_inc_dec_mutations(
    const Mutation_Count_Change &dst_add_count, bool is_dst_terminal,
    const Mutation_Count_Change &src_count_change, bool is_src_terminal,
    const MAT::Mutation &this_mut, int &score_change,
    Mutation_Count_Change_Collection &parent_mutation_change_out) {
    nuc_one_hot dst_increment=dst_add_count.get_incremented();
    nuc_one_hot dst_decrement=dst_add_count.get_decremented();
    nuc_one_hot src_increment=src_count_change.get_incremented();
    nuc_one_hot src_decrement=src_count_change.get_decremented();
    assert((!is_dst_terminal) || (!is_src_terminal));

    nuc_one_hot dbl_increment_raw=dst_increment&src_increment;
    nuc_one_hot dbl_decrement_raw=dst_decrement&src_decrement;
    nuc_one_hot singe_decrement_raw=(dst_decrement|src_decrement)&(~dbl_decrement_raw);
    nuc_one_hot singe_increment_raw=(dst_increment|src_increment)&(~dbl_increment_raw);

    nuc_one_hot dbl_inc_w_singe=dbl_increment_raw&(~dbl_decrement_raw);
    nuc_one_hot dbl_dec_w_singe=dbl_decrement_raw&(~dbl_increment_raw);

    nuc_one_hot double_increment_canceled_by_single=dbl_inc_w_singe&singe_decrement_raw;
    nuc_one_hot double_increment=dbl_inc_w_singe&(~double_increment_canceled_by_single);
    nuc_one_hot remaining_single_decrement=singe_decrement_raw&(~double_increment_canceled_by_single);

    nuc_one_hot double_decrement_canceled_by_single=dbl_dec_w_singe&singe_increment_raw;
    nuc_one_hot double_decrement=dbl_dec_w_singe&(~double_decrement_canceled_by_single);
    nuc_one_hot remaining_single_increment=singe_increment_raw&(~double_decrement_canceled_by_single);

    nuc_one_hot single_canceled=remaining_single_decrement&remaining_single_increment;
    nuc_one_hot single_increment=(remaining_single_increment&(~single_canceled))|double_increment_canceled_by_single;
    nuc_one_hot single_decrement=(remaining_single_decrement&(~single_canceled))|double_decrement_canceled_by_single;


    nuc_one_hot any_increment = single_increment | double_increment;
    nuc_one_hot any_decrement = single_decrement | double_decrement;

    nuc_one_hot all_mut =this_mut.get_all_major_allele();
    nuc_one_hot major_allele = all_mut;

    nuc_one_hot boundary1_dbl_inc =
        double_increment & this_mut.get_boundary1_one_hot();

    nuc_one_hot boundary2_dbl_inc =
        double_increment & this_mut.get_boundary2_one_hot();

    if (double_increment) {
        assert(!is_src_terminal);
        // Double increment
        if (all_mut & double_increment) {

            major_allele = all_mut & double_increment;
            // Nor incremented one are no longer tie
            nuc_one_hot major_allele_not_incremented =
                all_mut & (!double_increment);
            if (major_allele_not_incremented) {
                register_change(parent_mutation_change_out, this_mut,
                                major_allele_not_incremented, 0,major_allele);
            }
            // it is incrementing the same majority alleles, so their count
            // increases by 2, minor allele count shrink by 2
            if (src_count_change) {
                score_change--;
            }
            if (dst_add_count && (!is_dst_terminal)) {
                score_change--;
            }

            return major_allele;
        } else if (boundary1_dbl_inc) {

            nuc_one_hot major_allele_inc = all_mut & any_increment;
            nuc_one_hot major_allele_not_inc = all_mut & (!any_increment);

            // Majority allele incremented and twice incremented boundary 1
            // become tie
            major_allele = major_allele_inc | boundary1_dbl_inc;
            register_change(parent_mutation_change_out, this_mut,
                            major_allele_not_inc, boundary1_dbl_inc,major_allele);
            // Major allele count only increment by 1

            if (src_count_change || (dst_add_count && !is_dst_terminal)) {
                score_change--;
            }
            return major_allele;
        } else if (boundary2_dbl_inc) {
            nuc_one_hot major_allele_inc = all_mut & any_increment;
            if (major_allele_inc) {
                major_allele = major_allele_inc;
                nuc_one_hot major_allele_not_inc = all_mut & (!any_increment);
                register_change(parent_mutation_change_out, this_mut,
                                major_allele_not_inc, 0,major_allele);
                if (src_count_change &&
                    src_increment & major_allele_inc) {
                    score_change--;
                } else if (dst_add_count &&
                           dst_increment & major_allele_inc) {
                    score_change--;
                }
            } else {
                nuc_one_hot all_mut_decremented = all_mut & any_decrement;
                nuc_one_hot boundary1_allele_inc =
                    this_mut.get_boundary1_one_hot() & any_increment;
                major_allele = (all_mut ^ (~any_decrement)) |
                               boundary1_allele_inc | boundary2_dbl_inc;
                register_change(parent_mutation_change_out, this_mut,
                                all_mut_decremented,
                                boundary1_allele_inc | boundary2_dbl_inc,major_allele);
            }

            return major_allele;
        }
    }
    Mutation_Count_Change temp(src_count_change);
    temp.set_change(any_decrement, any_increment,0);
    int score_change_temp=0;
    major_allele = decrement_increment_mutation_count(
        this_mut, temp, parent_mutation_change_out, score_change_temp);
    // Fix up double decrement
    if(!parent_mutation_change_out.empty()){
    auto &last_change = parent_mutation_change_out.back();
    if (last_change.get_position() == this_mut.get_position()) {
        last_change.get_incremented() =
            last_change.get_incremented() & (~double_decrement);
        last_change.get_decremented() =
            last_change.get_decremented() | (double_decrement & all_mut);
    }
    }
    major_allele = major_allele & (~double_decrement);
    if (!(is_dst_terminal ^ is_src_terminal)) {
        score_change+=score_change_temp;
        // Total number of nodes didn't change, equivalent to swapping a child
        if ((!(all_mut & (~double_decrement)) &&
             (!(this_mut.get_boundary1_one_hot() & (~any_decrement)))) &&
            (!((all_mut | this_mut.get_boundary1_one_hot() |
                this_mut.get_boundary2_one_hot()) &
               any_increment))) {
            score_change++;
        }
    } else if (is_dst_terminal) {
        if (dst_increment!=major_allele) {
            score_change++;
        }
        score_change_internal_node(this_mut.get_mut_one_hot(), src_count_change, score_change);
        assert(!is_src_terminal);
    } else {
        assert(is_src_terminal);
        score_change_internal_node(this_mut.get_mut_one_hot(), dst_add_count, score_change);
        if(!(src_decrement&this_mut.get_mut_one_hot())){score_change--;}
    }
    return major_allele;
}
