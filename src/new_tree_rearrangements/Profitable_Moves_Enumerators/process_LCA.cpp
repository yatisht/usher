#include "process_each_node.hpp"
#include "process_individual_mutation.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <climits>
#include <cstdio>
struct process_LCA_two_children {
    bool is_src_left;
    bool is_src_terminal;
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out;
    int &parsimony_score_change;
    process_LCA_two_children(
        const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change)
        : is_src_left(src_branch == LCA->children[0]),
          is_src_terminal(is_src_terminal),
          LCA_parent_mutation_count_change_out(
              LCA_parent_mutation_count_change_out),
          parsimony_score_change(parsimony_score_change) {}
    void LCA_no_match(const MAT::Mutation &LCA_mut) {
        if (is_src_terminal) {
            nuc_one_hot major_allele = is_src_left
                                           ? LCA_mut.get_right_child_state()
                                           : LCA_mut.get_left_child_state();
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, 0, LCA_mut, major_allele);
        }
    }
    void LCA_dst_match(const MAT::Mutation &LCA_matching_mutation,
                       const Mutation_Count_Change &dst_add_count) {
        nuc_one_hot major_allele;
        if (is_src_terminal) {
            major_allele = dst_add_count.get_new_state();
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, 0, LCA_matching_mutation,
                major_allele);
        } else {
            int new_par_score = get_new_major_allele_binary_node(
                is_src_left ? LCA_matching_mutation.get_left_child_state()
                            : LCA_matching_mutation.get_right_child_state(),
                dst_add_count.get_new_state(), major_allele);
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, new_par_score,
                LCA_matching_mutation, major_allele);
        }
    }
    void LCA_src_dst_match(const MAT::Mutation &LCA_mutation,
                           const Mutation_Count_Change &src_count_change,
                           const Mutation_Count_Change &dst_count_change) {
        nuc_one_hot major_allele;
        if (is_src_terminal) {
            major_allele = dst_count_change.get_new_state();
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, 0, LCA_mutation,
                major_allele);
        } else {
            int new_mut_count = get_new_major_allele_binary_node(
                src_count_change.get_new_state(),
                dst_count_change.get_new_state(), major_allele);
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, new_mut_count,
                LCA_mutation, major_allele);
        }
    }
    void LCA_src_match(const MAT::Mutation &LCA_mutation,
                       const Mutation_Count_Change &src_count_change) {
        nuc_one_hot major_allele;
        if (is_src_terminal) {
            major_allele = is_src_left ? LCA_mutation.get_right_child_state()
                                       : LCA_mutation.get_left_child_state();
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, 0, LCA_mutation,
                major_allele);
        } else {
            int new_mut_count = get_new_major_allele_binary_node(
                src_count_change.get_new_state(),
                is_src_left ? LCA_mutation.get_right_child_state()
                            : LCA_mutation.get_left_child_state(),
                major_allele);
            parsimony_score_change += register_change_from_new_state(
                LCA_parent_mutation_count_change_out, new_mut_count,
                LCA_mutation, major_allele);
        }
    }
};

struct process_LCA_more_than_two {
    bool is_src_terminal;
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out;
    int &parsimony_score_change;
    process_LCA_more_than_two(
        const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change)
        : is_src_terminal(is_src_terminal),
          LCA_parent_mutation_count_change_out(
              LCA_parent_mutation_count_change_out),
          parsimony_score_change(parsimony_score_change) {}
    void LCA_no_match(const MAT::Mutation &LCA_mut) {
        if (is_src_terminal) {
            Mutation_Count_Change temp(LCA_mut,LCA_mut.get_mut_one_hot(),0,0,true);
            decrement_mutation_count(LCA_parent_mutation_count_change_out,
                                     LCA_mut, temp, parsimony_score_change);
            parsimony_score_change--;
        }
    }

    void LCA_dst_match(const MAT::Mutation &LCA_matching_mutation,
                       const Mutation_Count_Change &dst_add_count) {
        nuc_one_hot major_allele;
        if (is_src_terminal) {
            Mutation_Count_Change temp(LCA_matching_mutation,LCA_matching_mutation.get_mut_one_hot(), 0, 0,
                            true);
            major_allele = dbl_inc_dec_mutations(
                dst_add_count, temp, LCA_matching_mutation,
                parsimony_score_change, LCA_parent_mutation_count_change_out);
            parsimony_score_change--;
        } else {
            major_allele = decrement_increment_mutation_count(
                LCA_matching_mutation, dst_add_count,
                LCA_parent_mutation_count_change_out, parsimony_score_change);
        }
    }
    void LCA_src_dst_match(const MAT::Mutation &LCA_mutation,
                           const Mutation_Count_Change &src_count_change,
                           const Mutation_Count_Change &dst_count_change) {
        dbl_inc_dec_mutations(dst_count_change, src_count_change, LCA_mutation,
                              parsimony_score_change,
                              LCA_parent_mutation_count_change_out);
        parsimony_score_change += (0 - is_src_terminal);
    }
    void LCA_src_match(const MAT::Mutation &LCA_mutation,
                       const Mutation_Count_Change &src_count_change) {
        nuc_one_hot major_allele;
        if (is_src_terminal) {
            decrement_mutation_count(LCA_parent_mutation_count_change_out,
                                     LCA_mutation, src_count_change,
                                     parsimony_score_change);
            parsimony_score_change--;
        } else {
            decrement_increment_mutation_count(
                LCA_mutation, src_count_change,
                LCA_parent_mutation_count_change_out, parsimony_score_change);
        }
    }
};

template <typename Functor> class process_LCA {
    void LCA_no_match(int end_pos, range<MAT::Mutations_Collection> &LCA_mut,
                      Functor &functor) {
        while (LCA_mut && LCA_mut->get_position() < end_pos) {
            functor.LCA_no_match(*LCA_mut);
            LCA_mut++;
        }
    }
    void LCA_dst(range<MAT::Mutations_Collection> &LCA_mut,
                 range<Mutation_Count_Change_Collection> &dst_add_count_iter,
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
    void
    LCA_src_match(range<Mutation_Count_Change_Collection> &dst_add_count_iter,
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
    void proc_src(const Mutation_Count_Change &src_count_change,
                  range<MAT::Mutations_Collection> &LCA_mut,
                  range<Mutation_Count_Change_Collection> &dst_add_count_iter,
                  Functor &functor) {
        LCA_dst(LCA_mut, dst_add_count_iter, src_count_change.get_position(),
                functor);
        LCA_no_match(src_count_change.get_position(), LCA_mut, functor);
        if (LCA_mut &&
            LCA_mut->get_position() == src_count_change.get_position()) {
            LCA_src_match(dst_add_count_iter, src_count_change, *LCA_mut,
                          functor);
            LCA_mut++;
        } else {
            functor.parsimony_score_change +=
                functor.is_src_terminal
                    ? src_count_change.get_default_change_terminal()
                    : src_count_change.get_default_change_internal();
            if (dst_add_count_iter && dst_add_count_iter->get_position() ==
                                          src_count_change.get_position()) {
                functor.parsimony_score_change +=
                    dst_add_count_iter->get_default_change_internal();
                dst_add_count_iter++;
            }
        }
    }

    void iter_src(const Mutation_Count_Change_Collection &from_src_remove,
                  range<MAT::Mutations_Collection> &LCA_mut,
                  range<Mutation_Count_Change_Collection> &dst_add_count_iter,
                  Functor &functor) {
        for (const Mutation_Count_Change &src_count_change : from_src_remove) {
            proc_src(src_count_change, LCA_mut, dst_add_count_iter, functor);
        }
    }
    void iter_src(const MAT::Mutations_Collection &from_src_remove,
                  range<MAT::Mutations_Collection> &LCA_mut,
                  range<Mutation_Count_Change_Collection> &dst_add_count_iter,
                  Functor &functor) {
        for (const MAT::Mutation  &src_mut : from_src_remove) {
            Mutation_Count_Change src_count_change(src_mut,src_mut.get_all_major_allele(),0,0,true);
            proc_src(src_count_change, LCA_mut, dst_add_count_iter, functor);
        }
    }
  public:
    template<typename src_t>
    void operator()(
        const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
        const src_t &from_src_remove,
        const Mutation_Count_Change_Collection &from_dst_add,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change) {
        Functor functor(LCA, src_branch, is_src_terminal,
                        LCA_parent_mutation_count_change_out,
                        parsimony_score_change);
        range<MAT::Mutations_Collection> LCA_mut(LCA->mutations);
        range<Mutation_Count_Change_Collection> dst_add_count_iter(
            from_dst_add);
        iter_src(from_src_remove,LCA_mut,dst_add_count_iter,functor);
        LCA_dst(LCA_mut, dst_add_count_iter, INT_MAX, functor);
        LCA_no_match(INT_MAX, LCA_mut, functor);
    }
};
void get_LCA_mutation(
    const MAT::Node *LCA, const MAT::Node *src_branch, bool is_src_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change) {
    if (LCA->children.size() == 2) {
        process_LCA<process_LCA_two_children>()(
            LCA, src_branch, is_src_terminal, from_src_remove, from_dst_add,
            LCA_parent_mutation_count_change_out, parsimony_score_change);
    } else if (LCA->children.size() > 2) {
        if (is_src_terminal) {
            process_LCA<process_LCA_more_than_two>()(
                LCA, src_branch, is_src_terminal, src_branch->mutations, from_dst_add,
                LCA_parent_mutation_count_change_out, parsimony_score_change);
        }else {
            process_LCA<process_LCA_more_than_two>()(
                LCA, src_branch, is_src_terminal, from_src_remove, from_dst_add,
                LCA_parent_mutation_count_change_out, parsimony_score_change);
        }
    } else {
        assert(false);
    }
}