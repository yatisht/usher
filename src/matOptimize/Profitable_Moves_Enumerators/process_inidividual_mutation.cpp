#include "Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
static void register_change(Mutation_Count_Change_Collection &out,
                            const MAT::Mutation &pos,
                            nuc_one_hot decremented_allele,
                            nuc_one_hot incremented_allele) {
    out.emplace_back(pos,decremented_allele, incremented_allele);
}
nuc_one_hot decrement_mutation_count(Mutation_Count_Change_Collection &out,
                                     const MAT::Mutation &parent_mut,
                                     nuc_one_hot decremented_allele,
                                     int &score_change) {
    // Score_change assumes we are operating on internal nodes, so it is
    // removing a possibility of the state a child may take, not the whole child
    // altogether
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
                            0);
            // Parsimony score change: If the allele is originally tied with
            // another allele, then they have same number of children that can
            // only be that state, the removed mutation is not one of those,
            // otherwise it is removing a child altogether. Then the parsimony
            // score remain the number of children whose state can only be that
            // of the other major allele, which stay the same
            return not_decremented;
        } else {
            // One additional child can only have state of minor alleles, if it
            // is not complete removal
            score_change++;

            // all alleles originally tied got decremented
            if (boundary1_not_decremented) {
                // Then everything become tie with the old boundary1 allele that
                // aren't decremented
                nuc_one_hot major_allele=boundary1_not_decremented|all_nuc;
                register_change(out, parent_mut, 0, boundary1_not_decremented);
                return major_allele;
            }
            // everything got decremented, major alleles stay the same, boundary
            // 1 tie with boundary2, but no impact on parent of parent
            return all_nuc;
        }
    }
    // Decrementing minor alleles that is an alternative have no impact if it
    // doesn't remove the child
    /*if (is_terminal&&(!(change_in.get_par_state()&change_in.get_decremented()))) {
        score_change--;
    }*/
    return all_nuc;
}

nuc_one_hot increment_mutation_count(Mutation_Count_Change_Collection &out,
                                     const MAT::Mutation &parent_mut,
                                     nuc_one_hot incremented_allele,
                                     int &score_change) {
    // Break ties as before
    nuc_one_hot all_nuc =
        parent_mut.get_all_major_allele();
    nuc_one_hot effectively_incremented_alleles = all_nuc & incremented_allele;
    if (effectively_incremented_alleles) {
        // Not incremented alleles no longer in the tie
        nuc_one_hot not_incremented = effectively_incremented_alleles ^ all_nuc;

        // originally a mutation to some other allele, now can follow major
        // (potentially new )allele
        score_change--;

        if (not_incremented) {
            // Incremented major alleles are still major alleles
            // Not incremented major alleles was originally counted in count of
            // children mutations by parent of parent, they do not change the
            // number of children whose state that can only be minor alleles
            register_change(out, parent_mut, not_incremented, 0);
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
                        incremented_allele);
        // Either the parent is original allele, or the new major allele, number
        // of children who cant follow parent state is the same, and may
        // increment if this is a new child
        return tieing_alleles;
    }
    // Adding a minor allele, if adding new child, increment
    return parent_mut.get_all_major_allele();
}

nuc_one_hot decrement_increment_mutation_count(
    const MAT::Mutation &parent_mutation,nuc_one_hot decremented, nuc_one_hot incremented,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &score_change) {
    assert(!(decremented & incremented));
    //specialization when change_in is pure increment/decrement, which is the case for nodes under LCA
    //Remove node (base case)|| all state change is decrementing allele (induction) ->
    //  some of decremented allele decremented (if some major allele not decremented) OR
    //  some of non-decremented allele (allele count diff with major allele=1 ) incremented

    //Add node (base case)|| all state change is incrementing allele (induction) ->
    //  some of incremented allele incremented (none of incremented allele was originally major allele ) OR
    //  some of non-incremented allele (original major allele not incremented ) decremented

    if (incremented&&(!decremented)) {
        return increment_mutation_count(parent_node_mutation_count_change,parent_mutation,incremented,score_change);
    } else if (decremented&&(!incremented)) {
        return decrement_mutation_count(parent_node_mutation_count_change,parent_mutation,decremented,score_change);
    }

    //For mixed increment/decrement (rare case)
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
                            parent_node_mutation_count_change, parent_mutation, incremented,
                            score_change);
        return major_alleles;
    } else if (major_allele_decremented) {
        // this is adding to a minor allele while removing a major
        // allele
        if (boundary1_state & incremented || major_allele_not_decremented) {
            //assert(!major_allele_incremented);
            // tied alleles not decremented (no action needed)
            // and boundary 1 allele incremented (need to be incremented in
            // parent of parent counts) become new major allele
            // Tied major alleles decremented joins minor alleles (need to be
            // decremented in parent of parent counts)
            major_alleles =
                major_allele_not_decremented | (boundary1_state & incremented);
            register_change(parent_node_mutation_count_change, parent_mutation,
                            major_allele_decremented,
                            boundary1_state & incremented);
            // The number of children that can follow any of the major alleles
            // didn't change, as there are still same numbe of children that can
            // follow major_alleles
            return major_alleles;
        } else {
            // No thing interesting have been incremented
            major_alleles=decrement_mutation_count(parent_node_mutation_count_change,
                                                   parent_mutation, decremented,
                                                   score_change);
            return major_alleles;
        }
    } else {
        // Nothing happens to major alleles, for minor alleles, only increment
        // is interesting
        if (incremented) {
            return increment_mutation_count(parent_node_mutation_count_change,
                                            parent_mutation, incremented,
                                            score_change);
        }
    }
    return major_alleles;
}

static void score_change_internal_node(nuc_one_hot original_state,const Mutation_Count_Change & change, int& score_change ) {
    if (original_state&change.get_incremented()) {
        score_change--;
    }
    if(original_state&change.get_decremented()) {
        score_change++;
    }
}

nuc_one_hot dbl_inc_dec_mutations(
    const Mutation_Count_Change &dst_add_count,
    const Mutation_Count_Change &src_count_change,
    const MAT::Mutation &this_mut, int &score_change,
    Mutation_Count_Change_Collection &parent_mutation_change_out) {
    nuc_one_hot dst_increment=dst_add_count.get_incremented();
    nuc_one_hot dst_decrement=dst_add_count.get_decremented();
    nuc_one_hot src_increment=src_count_change.get_incremented();
    nuc_one_hot src_decrement=src_count_change.get_decremented();
    //Merging allele change toghether
    //Here we assume no allele will change in the same direction on both src branch and dst branch
    //Because on src branch, where a node is removed, some of its allele are decremented or some of the complement of
    //its allele are incremented at src_branch_node.
    //When it is added back at dst branch, some of its allele are incremented, some of the complement of its alleles are
    // decremented at dst branch node, so any individual allele cannot change in the same direction.
    nuc_one_hot any_increment = (src_increment&(~dst_decrement))|(dst_increment&(~src_decrement));
    nuc_one_hot any_decrement = (src_decrement&(~dst_increment))|(dst_decrement&(~src_increment));

    nuc_one_hot all_mut =this_mut.get_all_major_allele();
    nuc_one_hot major_allele = all_mut;

    major_allele = decrement_increment_mutation_count(
                       this_mut, any_decrement,any_increment, parent_mutation_change_out, score_change);
    return major_allele;
}
