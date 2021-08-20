#include "Bounded_Search.hpp"
#include "process_individual_mutation.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <tuple>
typedef Bounded_Mut_Change_Collection::iterator Bounded_Mut_Iter;
void check_parsimony_score_change_above_LCA(
    MAT::Node *LCA, int &parsimony_score_change,
    Mutation_Count_Change_Collection &parent_added,
    const std::vector<MAT::Node *> &node_stack_from_src,
    std::vector<MAT::Node *> &node_stack_above_LCA,
    Mutation_Count_Change_Collection &parent_of_parent_added,
    MAT::Node *ancestor);
void output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,
                   const std::vector<MAT::Node *> &node_stack_from_src,
                   std::vector<MAT::Node *> &node_stack_from_dst,
                   std::vector<MAT::Node *> &node_stack_above_LCA,
                   int radius_left);

    void output_LCA(Mutation_Count_Change_Collection
                        &allele_count_change_from_splitting_LCA,
                    int par_score_change, const src_side_info& src_side,int radius_left) {
        Mutation_Count_Change_Collection parent_of_parent_added;
        parent_of_parent_added.reserve(
            allele_count_change_from_splitting_LCA.size());
        std::vector<MAT::Node *> node_stack_above_LCA;
        check_parsimony_score_change_above_LCA(
            src_side.LCA, par_score_change, allele_count_change_from_splitting_LCA,
            src_side.node_stack_from_src, node_stack_above_LCA, parent_of_parent_added,
            src_side.LCA);
        std::vector<MAT::Node *> ignored;
        output_result(src_side.src, src_side.LCA, src_side.LCA, par_score_change, src_side.out,
                      src_side.node_stack_from_src, ignored, node_stack_above_LCA,
                      radius_left);
    }

static void add_mut(
    // IN
    const Mutation_Count_Change_W_Lower_Bound &in, const MAT::Node *parent_node,
    int radius_left, bool offsetable,
    // OUT
    // Go to sibling
    Bounded_Mut_Change_Collection &left_mut, int &left_lower_bound,
    Bounded_Mut_Change_Collection &right_mut, int &right_lower_bound,
    // Go to ancestor
    Bounded_Mut_Change_Collection &mutations_out) {
    mutations_out.push_back(in);
    mutations_out.back().offsetable = offsetable;
    auto res = mutations_out.back().to_parent(parent_node, radius_left);
    left_mut.push_back(std::get<0>(res));
    if (std::get<2>(res)) {
        left_lower_bound++;
    }
    right_mut.push_back(std::get<1>(res));
    if (std::get<3>(res)) {
        right_lower_bound++;
    }
}

void search_subtree_first_level(MAT::Node *node, MAT::Node *to_exclude,
                                const src_side_info &src_side, int radius_left,
                                const Bounded_Mut_Change_Collection &left,
                                int left_lower_bound,
                                const Bounded_Mut_Change_Collection &right,
                                int right_lower_bound) {
    for (auto child : node->children) {
        if (child == to_exclude) {
            continue;
        }
        if (child->dfs_end_index < to_exclude->dfs_index) {
#ifndef CHECK_BOUND
            if (left_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, left,
                                   left_lower_bound);
        } else {
            assert(child->dfs_index > to_exclude->dfs_end_index);
#ifndef CHECK_BOUND
            if (right_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, right,
                                   right_lower_bound);
        }
    }
}
/*
1. accumulate mutation vector of src to what is needed if placed as children of
parent
2. whether it is profitable to just place src as children of parent node
3. change to major allele set at this node
                parent
            /           \
            ? (2.split)   1. mutation vector if placed here
          /
        node (3. change to major allele set, also parsimony score change)
        /  \
    ...     also search other children
    /
    src
*/
static void upward_integrated(src_side_info &src_side,
                              Bounded_Mut_Change_Collection &src_mut_in,
                              int radius_left,
                              Bounded_Mut_Change_Collection &mut_out) {

    MAT::Node *node = src_side.LCA->parent;
    if (!node) {
        return;
    }

    // IN: alelle cnt change from this node (that will change the major allele
    // at parent node)
    auto src_allele_cnt_change_iter =
        src_side.allele_count_change_from_src.begin();
    auto src_allele_cnt_change_end =
        src_side.allele_count_change_from_src.end();

    // IN: Mutations on src node if placed as children of this node to keep its
    // and its descendant unchanged
    Bounded_Mut_Iter src_mut_iter = src_mut_in.begin();
    Bounded_Mut_Iter src_mut_end = src_mut_in.end();

    // OUT: Mutations on src node if placed as children of parent node, have two
    // copies, with same mutations, but different range Left side: nodes whose
    // DFS index is less than this node
    auto max_mut_size = node->mutations.size() + src_mut_in.size();
    Bounded_Mut_Change_Collection left_mut;
    left_mut.reserve(max_mut_size);
    int left_lower_bound = src_side.src_par_score_lower_bound;
    // Right side: nodes whose DFS idx > all descendants of this node
    Bounded_Mut_Change_Collection right_mut;
    right_mut.reserve(max_mut_size);
    int right_lower_bound = src_side.src_par_score_lower_bound;
    // mut_out is for going up to parent, with range covering the whole subtree
    // of parent node
    mut_out.reserve(max_mut_size);

    // allele cnt change from parent node
    Mutation_Count_Change_Collection allele_change_out;
    allele_change_out.reserve(
        std::min(node->mutations.size(),
                 src_side.allele_count_change_from_src.size()) +
        1);
    // par socre mchange from src node removal on src->parent node (LCA) branch
    int next_src_par_score = src_side.par_score_change_from_src_remove;
    // Less increments from old_LCA that can improve parsimony at ancestors of
    // parent node as lower bound
    int next_src_par_score_lower_bound_diff = 0;

    // Only include increment form adding src between parent node and its
    // parent, par score change on the src branch node is next_src_par_score
    // calculated in the same loop
    int par_score_change_split_LCA = 0;
    Mutation_Count_Change_Collection split_allele_count_change_out;

    for (const auto &mut : node->mutations) {
        // Calculate allele count change for parent (and new major allele at
        // parent node for splitting between parent node and its parent )
        nuc_one_hot major_allele = mut.get_all_major_allele();
        while (src_allele_cnt_change_iter->get_position() <
               mut.get_position()) {
            auto change =
                src_allele_cnt_change_iter->get_default_change_internal();
            next_src_par_score += change;
            src_allele_cnt_change_iter++;
        }
        if (src_allele_cnt_change_iter->get_position() == mut.get_position()) {
            major_allele = decrement_increment_mutation_count(
                mut, (Mutation_Count_Change)*src_allele_cnt_change_iter,
                allele_change_out, next_src_par_score);
            // this node incrementing some allele that can reduce parsimony
            // score among ancestors of  parent node
            if (src_allele_cnt_change_iter->get_incremented() &
                mut.get_sensitive_increment()) {
                next_src_par_score_lower_bound_diff--;
            }
        }

        // Get mutations if placed as children of parent node, and whether it is
        // profitable to place between parent node and parent of parent node
        while (src_mut_iter->get_position() < mut.get_position()) {
            add_mut(*src_mut_iter, node, radius_left, false, left_mut,
                    left_lower_bound, right_mut, right_lower_bound, mut_out);
            add_node_split(*src_mut_iter, split_allele_count_change_out,
                           par_score_change_split_LCA);
            src_mut_iter++;
        }
        if (src_mut_iter->get_position() == mut.get_position()) {
            // accumulate mutation for placement
            auto new_par_nuc = mut.get_par_one_hot();
            if (new_par_nuc != src_mut_iter->get_incremented()) {
                add_mut(*src_mut_iter, node, radius_left,
                        mut.get_sensitive_increment() &
                            src_mut_iter->get_incremented(),
                        left_mut, left_lower_bound, right_mut,
                        right_lower_bound, mut_out);
                mut_out.back().set_par_nuc(new_par_nuc);
            }
            // split LCA
            add_node_split(mut, src_mut_iter->get_incremented(), major_allele,
                           split_allele_count_change_out,
                           par_score_change_split_LCA);
            src_mut_iter++;
        } else {
            if (mut.is_valid()) {
                mut_out.emplace_back(mut, 0,mut.get_mut_one_hot());
                // get range
                mut_out.back().init(node);
            }
            add_node_split(mut, split_allele_count_change_out,
                           par_score_change_split_LCA);
        }
    }
    // Finish off allele count change
    while (src_allele_cnt_change_iter < src_allele_cnt_change_end) {
        auto change = src_allele_cnt_change_iter->get_default_change_internal();
        next_src_par_score += change;
        src_allele_cnt_change_iter++;
    }
    // Add sentinel
    allele_change_out.emplace_back();
    while (src_mut_iter != src_mut_end) {
        add_mut(*src_mut_iter, node, radius_left, false, left_mut,
                left_lower_bound, right_mut, right_lower_bound, mut_out);
        add_node_split(*src_mut_iter, split_allele_count_change_out,
                       par_score_change_split_LCA);
        src_mut_iter++;
    }

    search_subtree_first_level(node, src_side.LCA, src_side, radius_left,
                               left_mut, left_lower_bound, right_mut,
                               right_lower_bound);

    src_side.par_score_change_from_src_remove=next_src_par_score;
    src_side.src_par_score_lower_bound=next_src_par_score+next_src_par_score_lower_bound_diff;
    src_side.LCA=node;
    src_side.allele_count_change_from_src=allele_change_out;

    output_LCA(split_allele_count_change_out, par_score_change_split_LCA,
               src_side,radius_left);
    src_side.node_stack_from_src.push_back(node);

}