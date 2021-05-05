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
            dst_debug_in->position, dst_debug_in->par_nuc.back(),
            dst_debug_in->par_nuc.back(), 0, Mutation_Count_Change());
        dst_debug_in++;
    }
    if (dst_debug_in!=dst_debug_end&&position == dst_debug_in->position) {
        dst_debug_in++;
    }
    debug_above_LCA.emplace_back(
        position, par_nuc, major_allele, score_change,
        (!out.empty() && out.back().get_position() == position)
            ? out.back()
            : Mutation_Count_Change());
}

static bool LCA_place_mezzanine(
    MAT::Node *src_branch_node,
    const Mutation_Count_Change_Collection &dst_mutations,
    const Mutation_Count_Change_Collection &src_branch_node_mutations_altered,
    Mutation_Count_Change_Collection &out, int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> &src_dbg_states_in,
    std::vector<state_change_hist_dbg> &debug_above_LCA
#endif
) {
    bool have_shared=false;
    bool have_not_shared=false;

    auto src_branch_node_change_iter = src_branch_node_mutations_altered.begin();
    auto src_branch_node_change_end = src_branch_node_mutations_altered.end();

    auto added_iter = dst_mutations.begin();
    auto added_end = dst_mutations.end();

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    auto src_debug_in_iter = src_dbg_states_in.begin();
    auto src_debug_end = src_dbg_states_in.end();
#endif

    for (const MAT::Mutation &m : src_branch_node->mutations) {
        /*if (m.get_position()==13270) {
            fputc('ab',stderr);
        }*/
        while (added_iter != added_end &&
               added_iter->get_position() < m.get_position()) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    int old_parsimony_score=parsimony_score_change;
#endif
        /*if (added_iter->get_position()==13270) {
            fputc('ab',stderr);
        }*/
            assert(added_iter->get_par_state()==get_parent_state(src_branch_node, added_iter->get_position()));
            nuc_one_hot major_allele =
                added_iter->get_par_state() & added_iter->get_incremented();
            if(!major_allele){
                parsimony_score_change++;
                major_allele = added_iter->get_par_state() | added_iter->get_incremented();
                out.emplace_back(*added_iter);
                out.back().set_ori_state(added_iter->get_par_state());
                out.back().set_change(0, added_iter->get_incremented(),
                                  major_allele);
            }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            LCA_place_mezzanine_update_debug(
                src_debug_in_iter, src_debug_end, debug_above_LCA, out,
                added_iter->get_position(), major_allele,
                added_iter->get_par_state(), parsimony_score_change-old_parsimony_score);
#endif
            added_iter++;
        }
        nuc_one_hot major_allele = m.get_all_major_allele();
        nuc_one_hot new_major_allele = major_allele;
        int score_change = 0;
        if (src_branch_node_change_iter != src_branch_node_change_end &&
            src_branch_node_change_iter->get_position() == m.get_position()) {
            major_allele = src_branch_node_change_iter->get_new_state();
            src_branch_node_change_iter++;
        }
        if (added_iter != added_end &&
            m.get_position() == added_iter->get_position()) {
            //This is a newly added mutation while traversing to LCA
                if (added_iter->get_incremented()&major_allele) {
                    have_shared=true;
                    new_major_allele=added_iter->get_incremented()&major_allele;
                }else {
                    have_not_shared=true;
                    new_major_allele=added_iter->get_incremented()|major_allele;
                    score_change++;
                }
            added_iter++;
        }else {
            //The newly added node follow LCA
            have_not_shared=true;
            score_change+=get_new_major_allele_binary_node(major_allele, m.get_par_one_hot(), new_major_allele);
        }
        register_change_from_new_state(out, 0, m, new_major_allele);
        parsimony_score_change += score_change;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        LCA_place_mezzanine_update_debug(src_debug_in_iter, src_debug_end,
                                         debug_above_LCA, out, m.get_position(),
                                         new_major_allele, m.get_par_one_hot(),
                                         score_change);
#endif
    }
        return have_not_shared;
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
    int& parsimony_score_change,
    const std::vector<MAT::Node *> &node_stack_from_dst,
    Mutation_Count_Change_Collection &parent_added,
    const std::vector<MAT::Node *>& node_stack_from_src,
    std::vector<MAT::Node *>& node_stack_above_LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg>& debug_from_src,
    const std::vector<state_change_hist_dbg>& debug_from_dst,
    std::vector<state_change_hist_dbg>& debug_above_LCA
#endif
) {
    Mutation_Count_Change_Collection parent_of_parent_added;
    MAT::Node *ancestor = LCA;
    if (dst==LCA) {
        if(!LCA_place_mezzanine(node_stack_from_src.back(), mutations, root_mutations_altered, parent_of_parent_added, parsimony_score_change, debug_from_src, debug_above_LCA)){
            return nullptr;
        }
    }else{
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(debug_from_dst.empty() ||
           debug_from_dst[0].mutation_score_change.size() ==
               node_stack_from_dst.size());
    std::vector<LCA_merged_states> merged_states;
    prep_LCA_checker(debug_from_src, debug_from_dst, merged_states);
#endif
        bool is_src_terminal = src->parent == LCA;
        bool is_dst_terminal = dst == LCA;
        if ((!(root_mutations_altered.empty() && parent_added.empty())) ||
            is_src_terminal || is_dst_terminal) {
            if (LCA->children.size()==2) {
            get_LCA_mutation_binary_node(LCA, is_src_terminal ? src : node_stack_from_src.back(),is_src_terminal ,is_dst_terminal ? dst : node_stack_from_dst.back(),
                             is_dst_terminal, root_mutations_altered,
                             is_dst_terminal ? mutations : parent_added,
                             parent_of_parent_added, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                             ,
                             debug_above_LCA, merged_states,
                             node_stack_above_LCA
#endif
            );
            }else{
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
            );}
        }else {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        LCA_no_change_update_debug(LCA, debug_above_LCA, merged_states);
#endif
        }
        ancestor = LCA->parent;
    }
    node_stack_above_LCA.push_back(LCA);
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
        ancestor=node_stack_above_LCA.back();
    }

    return ancestor;
}