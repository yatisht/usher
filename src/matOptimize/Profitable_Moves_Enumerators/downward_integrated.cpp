#include "Bounded_Search.hpp"
#include "Profitable_Moves_Enumerators.hpp"
#include "process_each_node.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
typedef Bounded_Mut_Change_Collection::const_iterator Bounded_Mut_Iter;
static void
add_remaining_dst_to_LCA_nodes(MAT::Node *cur, const MAT::Node *LCA,
                               std::vector<MAT::Node *> &dst_stack) {
    while (cur != LCA) {
        // assert(dst_stack.empty()||dst_stack.back()!=cur);
        dst_stack.push_back(cur);
        cur = cur->parent;
    }
}

void output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,
                   const std::vector<MAT::Node *> &node_stack_from_src,
                   std::vector<MAT::Node *> &node_stack_from_dst,
                   std::vector<MAT::Node *> &node_stack_above_LCA,
                   int radius_left);
static void output_not_LCA(Mutation_Count_Change_Collection &parent_added,
                    MAT::Node *dst_node, int parsimony_score_change,
                    int lower_bound, const src_side_info &src_side,
                    int radius) {
#ifndef CHECK_BOUND
    if (use_bound&&lower_bound > src_side.output.score_change) {
        return;
    }
#else
    Mutation_Count_Change_Collection parent_backup(parent_added);
#endif
    std::vector<MAT::Node *> node_stack_from_dst;
    Mutation_Count_Change_Collection parent_of_parent_added;
    parent_of_parent_added.reserve(parent_added.size());
    node_stack_from_dst.push_back(dst_node);
    auto this_node = dst_node->parent;
    if (dst_node->dfs_index==58) {
        fputc('a', stderr);
    }
    while (this_node != src_side.LCA) {
        parent_of_parent_added.clear();
        get_intermediate_nodes_mutations(this_node, node_stack_from_dst.back(),
                                         parent_added, parent_of_parent_added,
                                         parsimony_score_change);
        node_stack_from_dst.push_back(this_node);
        parent_added.swap(parent_of_parent_added);
        if (parent_added.empty()) {
            add_remaining_dst_to_LCA_nodes(this_node->parent, src_side.LCA,
                                           node_stack_from_dst);
            break;
        }
        this_node = this_node->parent;
    }
    std::vector<MAT::Node *> node_stack_above_LCA;
    parent_of_parent_added.reserve(
        parent_added.size() + src_side.allele_count_change_from_src.size());
    // Adjust LCA node and above
    parent_of_parent_added.clear();
    bool is_src_terminal = src_side.src->parent == src_side.LCA;
    if ((!(src_side.allele_count_change_from_src.empty() &&
           parent_added.empty())) ||
        is_src_terminal) {
        get_LCA_mutation(src_side.LCA,
                         is_src_terminal ? src_side.src
                                         : src_side.node_stack_from_src.back(),
                         is_src_terminal, src_side.allele_count_change_from_src,
                         parent_added, parent_of_parent_added,
                         parsimony_score_change);
    }
    node_stack_above_LCA.push_back(src_side.LCA);
    parent_added.swap(parent_of_parent_added);
    parent_of_parent_added.clear();
    check_parsimony_score_change_above_LCA(
        src_side.LCA, parsimony_score_change, parent_added,
        src_side.node_stack_from_src, node_stack_above_LCA,
        parent_of_parent_added, src_side.LCA->parent);
    #ifdef CHECK_PAR_MAIN
    output_t temp;
    auto refout=individual_move(src_side.src, dst_node, src_side.LCA, temp);
    if (refout>=0) {
        assert(parsimony_score_change>=0);
    }else {
        assert(parsimony_score_change==refout);
    }
    #endif
    #ifdef CHECK_BOUND
    assert((!use_bound)||(parsimony_score_change >= lower_bound)||(parsimony_score_change>=0));
#endif
    output_result(src_side.src, dst_node, src_side.LCA, parsimony_score_change,
                  src_side.out, src_side.node_stack_from_src,
                  node_stack_from_dst, node_stack_above_LCA, radius);
}

static void add_mut(const Mutation_Count_Change_W_Lower_Bound &mut, int radius_left,
             MAT::Node *node, int &descendant_lower_bound,
             Bounded_Mut_Change_Collection &mut_out) {
    mut_out.push_back(mut);
    if (use_bound) {
    if (!mut_out.back().to_decendent(node, radius_left)) {
        descendant_lower_bound++;
    }
    }
}

static int downward_integrated(MAT::Node *node, int radius_left,
                        const Bounded_Mut_Change_Collection &from_parent,
                        Bounded_Mut_Change_Collection &mut_out,
                        const src_side_info &src_side
#ifdef CHECK_BOUND
                        ,
                        int prev_lower_bound
#endif

) {
    Bounded_Mut_Iter iter = from_parent.begin();
    Bounded_Mut_Iter end = from_parent.end();
    if (node->dfs_index==58) {
        fputc('a', stderr);
    }
    Mutation_Count_Change_Collection split_allele_cnt_change;
    split_allele_cnt_change.reserve(
        std::min(from_parent.size(), node->mutations.size()));
    int par_score_from_split = src_side.par_score_change_from_src_remove;
    int par_score_from_split_lower_bound = src_side.src_par_score_lower_bound;

    int descendant_lower_bound = src_side.src_par_score_lower_bound;

    for (const auto &mut : node->mutations) {
        while (iter->get_position() < mut.get_position()) {
            auto have_not_shared = add_node_split(
                *iter, split_allele_cnt_change, par_score_from_split);
            if (have_not_shared && (!iter->offsetable)) {
                par_score_from_split_lower_bound++;
            }
            add_mut(*iter, radius_left, node, descendant_lower_bound, mut_out);
            mut_out.back().offsetable=false;
            iter++;
        }
        if (iter->get_position() == mut.get_position()) {
            auto have_not_shared = add_node_split(
                mut, mut.get_all_major_allele(), iter->get_incremented(),
                split_allele_cnt_change, par_score_from_split);
            if (have_not_shared && (!iter->offsetable)) {
                par_score_from_split_lower_bound++;
            }

            //add_mut(*iter, radius_left, node, descendant_lower_bound, mut_out);
            mut_out.push_back(*iter);
            mut_out.back().set_par_nuc(mut.get_mut_one_hot());
            if (!mut_out.back().to_decendent(node, radius_left)) {
                descendant_lower_bound++;
            }
            mut_out.back().offsetable=mut.get_sensitive_increment()&iter->get_incremented();
            //assert((!mut_out.back().offsetable)||mut_out.back().have_content);
            iter++;
        } else {
            if (mut.is_valid()) {
                mut_out.emplace_back(mut, 0, mut.get_par_one_hot());
                mut_out.back().set_par_nuc(mut.get_mut_one_hot());
                mut_out.back().offsetable=mut.get_sensitive_increment()&mut.get_par_one_hot();
                if (!mut_out.back().init(node, radius_left)) {
                    descendant_lower_bound++;
                }
                //assert((!mut_out.back().offsetable)||mut_out.back().have_content);
            }
            add_node_split(mut, split_allele_cnt_change, par_score_from_split);
        }
    }
    while (iter ->get_position()!=INT_MAX) {
        auto have_not_shared = add_node_split(*iter, split_allele_cnt_change,
                                              par_score_from_split);
        if (have_not_shared & (!iter->offsetable)) {
            par_score_from_split_lower_bound++;
        }
        add_mut(*iter, radius_left, node, descendant_lower_bound, mut_out);
        mut_out.back().offsetable=false;
        iter++;
    }
    mut_out.emplace_back();
#ifdef CHECK_BOUND
    assert((!use_bound)||par_score_from_split_lower_bound >= prev_lower_bound);
#endif
    output_not_LCA(split_allele_cnt_change, node, par_score_from_split,
                   par_score_from_split_lower_bound, src_side, radius_left);
    return descendant_lower_bound;
}
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                    int radius_left,
                    const Bounded_Mut_Change_Collection &par_muts,
                    int lower_bound) {
    Bounded_Mut_Change_Collection muts;
    lower_bound = downward_integrated(node, radius_left, par_muts, muts,
                                      src_side, lower_bound);
    if (!radius_left) {
        return;
    }
#ifndef CHECK_BOUND
    if (use_bound&&lower_bound > src_side.out.score_change) {
        return;
    }
#endif
    for (auto child : node->children) {
        search_subtree_bounded(child, src_side, radius_left-1, muts, lower_bound);
    }
}
