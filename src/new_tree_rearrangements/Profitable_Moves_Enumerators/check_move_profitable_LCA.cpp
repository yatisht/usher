#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"

static void LCA_place_mezzanine_update_debug(
    std::vector<state_change_hist_dbg>::const_iterator &dst_debug_in,
    std::vector<state_change_hist_dbg>::const_iterator &dst_debug_end,
    std::vector<state_change_hist_dbg> &debug_above_LCA,
    Mutation_Count_Change_Collection &out, int position,
    nuc_one_hot major_allele, nuc_one_hot par_nuc, int score_change) {
    while (dst_debug_in != dst_debug_end && dst_debug_in->position < position) {
        debug_above_LCA.emplace_back(
            dst_debug_in->position, dst_debug_in->par_nuc,
            dst_debug_in->par_nuc, 0, Mutation_Count_Change());
        dst_debug_in++;
    }
    if (position == dst_debug_in->position) {
        dst_debug_in++;
    }
    debug_above_LCA.emplace_back(
        dst_debug_in->position, par_nuc, major_allele, score_change,
        (!out.empty() && out.back().get_position() == position)
            ? out.back()
            : Mutation_Count_Change());
}

static void LCA_place_mezzanine(
    MAT::Node *LCA,
    const Mutation_Count_Change_Collection &dst_unique_mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    Mutation_Count_Change_Collection &out, int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> &src_dbg_states_in,
    std::vector<state_change_hist_dbg> &debug_above_LCA
#endif
) {
    auto LCA_change_iter = root_mutations_altered.begin();
    auto LCA_change_end = root_mutations_altered.end();

    auto added_iter = dst_unique_mutations.begin();
    auto added_end = dst_unique_mutations.end();

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    auto src_debug_in_iter = src_dbg_states_in.begin();
    auto src_debug_end = src_dbg_states_in.end();
#endif

    for (const MAT::Mutation &m : LCA->mutations) {
        while (added_iter != added_end &&
               added_iter->get_position() < m.get_position()) {
            parsimony_score_change++;
            out.emplace_back(*added_iter);
            nuc_one_hot major_allele =
                added_iter->get_par_state() | added_iter->get_incremented();
            out.back().set_change(0, added_iter->get_incremented(),
                                  major_allele);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            LCA_place_mezzanine_update_debug(
                src_debug_in_iter, src_debug_end, debug_above_LCA, out,
                added_iter->get_position(), major_allele,
                added_iter->get_par_state(), 1);
#endif
        }
        nuc_one_hot major_allele = m.get_all_major_allele();
        nuc_one_hot new_major_allele = major_allele;
        int score_change = 0;
        bool overriden = false;
        if (LCA_change_iter != LCA_change_end &&
            LCA_change_iter->get_position() == m.get_position()) {
            major_allele = LCA_change_iter->get_new_state();
            LCA_change_iter++;
            overriden = true;
        }
        if (added_iter != added_end &&
            m.get_position() == added_iter->get_position()) {
            score_change = get_new_major_allele_binary_node(
                major_allele, added_iter->get_incremented(), new_major_allele);
            register_change_from_new_state(out, score_change, m,
                                           new_major_allele);
        } else if (overriden) {
            score_change = get_new_major_allele_binary_node(
                LCA_change_iter->get_ori_state(),
                LCA_change_iter->get_new_state(), new_major_allele);
            register_change_from_new_state(out, score_change, m,
                                           new_major_allele);
        }
        parsimony_score_change += score_change;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        LCA_place_mezzanine_update_debug(src_debug_in_iter, src_debug_end,
                                         debug_above_LCA, out, m.get_position(),
                                         new_major_allele, m.get_par_one_hot(),
                                         score_change);
#endif
    }
}

static void
LCA_no_change_update_debug(const MAT::Node *LCA,
                           std::vector<state_change_hist_dbg> &debug_above_LCA,
                           std::vector<LCA_merged_states> &merged_states) {
    auto iter = LCA->mutations.begin();
    auto end = LCA->mutations.end();
    for (auto &in : merged_states) {
        while (iter != end && iter->get_position() < in.position) {
            iter++;
        }
        if (iter != end && iter->get_position() == in.position) {
            debug_above_LCA.emplace_back(
                iter->get_position(), iter->get_par_one_hot(),
                iter->get_all_major_allele(), 0, Mutation_Count_Change());
        } else {
            debug_above_LCA.emplace_back(in.position, in.par_allele,
                                         in.par_allele, 0,
                                         Mutation_Count_Change());
        }
    }
}
MAT::Node *check_move_profitable_LCA(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change,
    const std::vector<MAT::Node *> &node_stack_from_dst,
    Mutation_Count_Change_Collection &parent_added,
    const std::vector<MAT::Node *>& node_stack_from_src, bool have_shared,
    const Mutation_Count_Change_Collection &dst_unique_mutations,
    std::vector<MAT::Node *> node_stack_above_LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> debug_from_src,
    const std::vector<state_change_hist_dbg> debug_from_dst,
    std::vector<state_change_hist_dbg> debug_above_LCA
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(debug_from_dst.empty() ||
           debug_from_dst[0].mutation_score_change.size() ==
               node_stack_from_dst.size());
    std::vector<LCA_merged_states> merged_states;
    prep_LCA_checker(debug_from_src, debug_from_dst, merged_states);
#endif
    MAT::Node *ancestor = LCA->parent;
    Mutation_Count_Change_Collection parent_of_parent_added;
    if (dst == LCA && have_shared) {
        LCA_place_mezzanine(LCA, dst_unique_mutations, root_mutations_altered, parent_of_parent_added, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        , debug_from_src, debug_above_LCA
#endif
        );
    } else {
        bool is_src_terminal = src->parent == LCA;
        bool is_dst_terminal = dst == LCA;
        if ((!(root_mutations_altered.empty() && parent_added.empty())) ||
            is_src_terminal || is_dst_terminal) {
            get_LCA_mutation(LCA, is_src_terminal ? src : 0,
                             is_dst_terminal ? dst : 0, root_mutations_altered,
                             is_dst_terminal ? mutations : parent_added,
                             parent_of_parent_added, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                             ,
                             debug_above_LCA, merged_states,
                             node_stack_above_LCA,
                             is_src_terminal ? src : node_stack_from_src.back(),
                             is_dst_terminal ? dst : node_stack_from_dst.back()
#endif
            );
        }
    }
    node_stack_above_LCA.push_back(LCA);
    ancestor = LCA->parent;
    parent_added = std::move(parent_of_parent_added);
    while (ancestor && (!parent_added.empty())) {
        get_intermediate_nodes_mutations(ancestor, parent_added,
                                         parent_of_parent_added,
                                         parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                         ,
                                         debug_above_LCA, node_stack_above_LCA
#endif
        );
        node_stack_above_LCA.push_back(ancestor);
        parent_added = std::move(parent_of_parent_added);
        ancestor = ancestor->parent;
    }
    if (!ancestor) {
        for (auto &a : parent_added) {
            parsimony_score_change += a.get_default_change_internal();
        }
    }

    else {
        node_stack_above_LCA.push_back(LCA);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        LCA_no_change_update_debug(LCA, debug_above_LCA, merged_states);
#endif
    }
    return ancestor;
}