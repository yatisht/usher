#include "process_each_node.hpp"
#include "process_individual_mutation.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <climits>
#include <cstdio>
struct ignore_remaining_LCA {};
struct use_remaining_LCA {};
// src_branch_node is not src, and LCA have 2 children, reconstitute that 2
// children with their new state to get new major allele and parsimony score
// from scratch
// More than 2 children with src node being direct children of LCA,
// so all alleles cared by src are decremented,
// calling respective process_individual_mutation layer functions
struct process_LCA_more_than_two_src_terminal {
    typedef use_remaining_LCA remaining_LCA_useful;
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out;
    int &parsimony_score_change;
    process_LCA_more_than_two_src_terminal(
        const MAT::Node *LCA,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change)
        : LCA_parent_mutation_count_change_out(
              LCA_parent_mutation_count_change_out),
          parsimony_score_change(parsimony_score_change) {}
    void LCA_no_match(const MAT::Mutation &LCA_mut) {
        Mutation_Count_Change temp(LCA_mut, LCA_mut.get_mut_one_hot(), 0);
        decrement_mutation_count(LCA_parent_mutation_count_change_out, LCA_mut,
                                 temp, parsimony_score_change);
        parsimony_score_change--;
    }

    void LCA_dst_match(const MAT::Mutation &LCA_matching_mutation,
                       const Mutation_Count_Change &dst_add_count) {
        nuc_one_hot major_allele;
        Mutation_Count_Change temp(LCA_matching_mutation,
                                   LCA_matching_mutation.get_mut_one_hot(), 0);
        major_allele = dbl_inc_dec_mutations(
                           dst_add_count, temp, LCA_matching_mutation, parsimony_score_change,
                           LCA_parent_mutation_count_change_out);
        parsimony_score_change--;
    }
    void LCA_src_dst_match(const MAT::Mutation &LCA_mutation,
                           const Mutation_Count_Change &src_count_change,
                           const Mutation_Count_Change &dst_count_change) {
        dbl_inc_dec_mutations(dst_count_change, src_count_change, LCA_mutation,
                              parsimony_score_change,
                              LCA_parent_mutation_count_change_out);
        parsimony_score_change--;
    }
    void LCA_src_match(const MAT::Mutation &LCA_mutation,
                       const Mutation_Count_Change &src_count_change) {
        nuc_one_hot major_allele;
        decrement_mutation_count(LCA_parent_mutation_count_change_out,
                                 LCA_mutation, src_count_change,
                                 parsimony_score_change);
        parsimony_score_change--;
    }
};
// Most normal case of LCA have 2 children, and src_branch_node is not being
// removed, basically just
// a fascade to process_individual_mutation level functions.
struct process_LCA_more_than_two_src_not_terminal {
    typedef ignore_remaining_LCA remaining_LCA_useful;
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out;
    int &parsimony_score_change;
    process_LCA_more_than_two_src_not_terminal(
        const MAT::Node *LCA,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change)
        : LCA_parent_mutation_count_change_out(
              LCA_parent_mutation_count_change_out),
          parsimony_score_change(parsimony_score_change) {}
    void LCA_no_match(const MAT::Mutation &LCA_mut) {}

    void LCA_dst_match(const MAT::Mutation &LCA_matching_mutation,
                       const Mutation_Count_Change &dst_add_count) {
        nuc_one_hot major_allele;
        major_allele = decrement_increment_mutation_count(
                           LCA_matching_mutation, dst_add_count,
                           LCA_parent_mutation_count_change_out, parsimony_score_change);
    }
    void LCA_src_dst_match(const MAT::Mutation &LCA_mutation,
                           const Mutation_Count_Change &src_count_change,
                           const Mutation_Count_Change &dst_count_change) {
        dbl_inc_dec_mutations(dst_count_change, src_count_change, LCA_mutation,
                              parsimony_score_change,
                              LCA_parent_mutation_count_change_out);
    }
    void LCA_src_match(const MAT::Mutation &LCA_mutation,
                       const Mutation_Count_Change &src_count_change) {
        nuc_one_hot major_allele;
        decrement_increment_mutation_count(LCA_mutation, src_count_change,
                                           LCA_parent_mutation_count_change_out,
                                           parsimony_score_change);
    }
};
// This template does the real work of merging LCA sensitive loci, src chamged
// alleles, and dst changed alleles toghether to get Fitch set change and
// parsimony score change
template <typename Functor> class process_LCA {
    // Sensitive loci at LCA node, but have no explicit change in src or dst,
    // only useful for removing src, when the mut_nuc (a major allele) is
    // decremented
    void LCA_no_match(int end_pos, range<MAT::Mutation> &LCA_mut,
                      Functor &functor) {
        while (LCA_mut && LCA_mut->get_position() < end_pos) {
            functor.LCA_no_match(*LCA_mut);
            LCA_mut++;
        }
    }
    // Basically the same as above, but add tag for translating to noop if it is
    // not removing src_branch_node
    void LCA_no_match_remaining(range<MAT::Mutation> &LCA_mut,
                                Functor &functor, ignore_remaining_LCA tag) {}
    void LCA_no_match_remaining(range<MAT::Mutation> &LCA_mut,
                                Functor &functor, use_remaining_LCA tag) {
        while (LCA_mut) {
            functor.LCA_no_match(*LCA_mut);
            LCA_mut++;
        }
    }
    // process change on dst branch node before a change in src branch node is
    // encountered
    void LCA_dst(range<MAT::Mutation> &LCA_mut,
                 range<Mutation_Count_Change> &dst_add_count_iter,
                 int end_pos, Functor &functor) {

        while (dst_add_count_iter &&
                dst_add_count_iter->get_position() < end_pos) {
            LCA_no_match(dst_add_count_iter->get_position(), LCA_mut, functor);
            if (LCA_mut &&
                    LCA_mut->get_position() == dst_add_count_iter->get_position()) {
                functor.LCA_dst_match(*LCA_mut, *dst_add_count_iter);
                LCA_mut++;
            } else {
                functor.parsimony_score_change +=
                    dst_add_count_iter->get_default_change_internal();
            }
            dst_add_count_iter++;
        }
    }
    // process src_branch_node fitch set change that coincide with a sensitive
    // loci at LCA
    void
    LCA_src_match(range<Mutation_Count_Change> &dst_add_count_iter,
                  const Mutation_Count_Change &src_count_change,
                  const MAT::Mutation &LCA_mutation, Functor &functor) {
        if (dst_add_count_iter && dst_add_count_iter->get_position() ==
                src_count_change.get_position()) {
            functor.LCA_src_dst_match(LCA_mutation, src_count_change,
                                      *dst_add_count_iter);
            dst_add_count_iter++;
        } else {
            functor.LCA_src_match(LCA_mutation, src_count_change);
        }
    }
    // Specialization for handing changes on src_branch_node that do not
    // coincide with sensitive loci at LCA, depending on whether it is removed
    void src_unmatch_default(const Mutation_Count_Change &src_count_change,
                             Functor &functor, use_remaining_LCA tag) {
        functor.parsimony_score_change +=
            src_count_change.get_default_change_terminal();
    }
    void src_unmatch_default(const Mutation_Count_Change &src_count_change,
                             Functor &functor, ignore_remaining_LCA tag) {
        functor.parsimony_score_change +=
            src_count_change.get_default_change_internal();
    }
    // Iterate other iterators (LCA sensitive loci, changes on dst branch node)
    // to process src branch change
    template <typename Tag_t>
    void proc_src(const Mutation_Count_Change &src_count_change,
                  range<MAT::Mutation> &LCA_mut,
                  range<Mutation_Count_Change> &dst_add_count_iter,
                  Functor &functor, Tag_t tag) {
        // process fitch set change on dst branch node before this changed
        // allele from src branch node
        LCA_dst(LCA_mut, dst_add_count_iter, src_count_change.get_position(),
                functor);
        // same for smaller LCA sensitive loci
        LCA_no_match(src_count_change.get_position(), LCA_mut, functor);
        if (LCA_mut &&
                LCA_mut->get_position() == src_count_change.get_position()) {
            LCA_src_match(dst_add_count_iter, src_count_change, *LCA_mut,
                          functor);
            LCA_mut++;
        } else {
            src_unmatch_default(src_count_change, functor, tag);
            // also a mismatch for dst if src is mismatch, and they won't
            // interact to change the major allele, as the allele change tend
            // toward cancelling
            if (dst_add_count_iter && dst_add_count_iter->get_position() ==
                    src_count_change.get_position()) {
                functor.parsimony_score_change +=
                    dst_add_count_iter->get_default_change_internal();
                dst_add_count_iter++;
            }
        }
    }
    // loop for processing src branch node change when src branch node is
    // removed,where all its major alleles are removed
    template <typename Tag_t>
    void iter_src(const Mutation_Count_Change_Collection &from_src_remove,
                  range<MAT::Mutation> &LCA_mut,
                  range<Mutation_Count_Change> &dst_add_count_iter,
                  Functor &functor, Tag_t tag) {
        for (const Mutation_Count_Change &src_count_change : from_src_remove) {
            proc_src(src_count_change, LCA_mut, dst_add_count_iter, functor,
                     tag);
        }
    }
    // for when src branch node is not removed
    template <typename Tag_t>
    void iter_src(const MAT::Mutations_Collection &from_src_remove,
                  range<MAT::Mutation> &LCA_mut,
                  range<Mutation_Count_Change> &dst_add_count_iter,
                  Functor &functor, Tag_t tag) {
        for (const MAT::Mutation &src_mut : from_src_remove) {
            Mutation_Count_Change src_count_change(
                src_mut, src_mut.get_all_major_allele(), 0);
            proc_src(src_count_change, LCA_mut, dst_add_count_iter, functor,
                     tag);
        }
    }

  public:
    template <typename src_t>
    void operator()(
        const MAT::Node *LCA,
        const src_t &from_src_remove,
        const Mutation_Count_Change_Collection &from_dst_add,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change) {
        Functor functor(LCA, LCA_parent_mutation_count_change_out,
                        parsimony_score_change);
        range<MAT::Mutation> LCA_mut(LCA->mutations);
        range<Mutation_Count_Change> dst_add_count_iter(
            from_dst_add);
        iter_src(from_src_remove, LCA_mut, dst_add_count_iter, functor,
                 typename Functor::remaining_LCA_useful());
        LCA_dst(LCA_mut, dst_add_count_iter, INT_MAX, functor);
        LCA_no_match_remaining(LCA_mut, functor,
                               typename Functor::remaining_LCA_useful());
    }
};
// This just dispatch the appropriate function based on number of children on
// LCA, and whether src_branch_node is removed
void get_LCA_mutation(
    const MAT::Node *LCA, const MAT::Node *src,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change) {
    assert(LCA->children.size()>1);
    if (src->parent==LCA) {
        process_LCA<process_LCA_more_than_two_src_terminal>()(
            LCA,  src->mutations, from_dst_add,
            LCA_parent_mutation_count_change_out, parsimony_score_change);
    } else {
        process_LCA<process_LCA_more_than_two_src_not_terminal>()(
            LCA,  from_src_remove, from_dst_add,
            LCA_parent_mutation_count_change_out, parsimony_score_change);
    }

}