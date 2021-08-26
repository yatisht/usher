#include "Bounded_Search.hpp"
#include "process_individual_mutation.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
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

void output_LCA(
    Mutation_Count_Change_Collection &allele_count_change_from_splitting_LCA,
    int par_score_change, const src_side_info &src_side, int radius_left) {
    Mutation_Count_Change_Collection parent_of_parent_added;
    parent_of_parent_added.reserve(
        allele_count_change_from_splitting_LCA.size());
    auto actual_LCA=src_side.LCA->parent;
    if (!actual_LCA) {
        return;
    }
    std::vector<MAT::Node *> node_stack_above_LCA{actual_LCA};
    check_parsimony_score_change_above_LCA(
        actual_LCA, par_score_change, allele_count_change_from_splitting_LCA,
        src_side.node_stack_from_src, node_stack_above_LCA,
        parent_of_parent_added, actual_LCA);
    std::vector<MAT::Node *> ignored;
    #ifdef CHECK_PAR_MAIN
    output_t temp;
    auto ref_out=individual_move(src_side.src, actual_LCA, actual_LCA, temp);
    if (ref_out>=0) {
        if (par_score_change<0) {
            fprintf(stderr, "%s\n",src_side.src->identifier.c_str());        
        }
        assert(par_score_change>=0);
    }else {
        if (ref_out!=par_score_change) {
            fprintf(stderr, "%s\n",src_side.src->identifier.c_str());
        }
        assert(ref_out==par_score_change);
    }
    #endif
    output_result(src_side.src, actual_LCA, actual_LCA, par_score_change,
                  src_side.out, src_side.node_stack_from_src, ignored,
                  node_stack_above_LCA, radius_left);
}

static void split_children(const MAT::Node *&parent_node, int &radius_left,
               Bounded_Mut_Change_Collection &left_mut, int &left_lower_bound,
               Bounded_Mut_Change_Collection &right_mut, int &right_lower_bound,
               Mutation_Count_Change_W_Lower_Bound &in) {
    if (use_bound) {
        auto res = in.to_parent(parent_node, radius_left);
        left_mut.push_back(std::get<0>(res));
        bool count=!(in.get_par_state()&in.get_incremented());
        if ((!std::get<0>(res).have_content)&&count) {
            left_lower_bound++;
        }
        right_mut.push_back(std::get<1>(res));
        if ((!std::get<1>(res).have_content)&&count) {
            right_lower_bound++;
        }
    }
}
static void add_mut(
    // IN
    const Mutation_Count_Change_W_Lower_Bound &in, const MAT::Node *parent_node,
    int radius_left,
    // OUT
    // Go to sibling
    Bounded_Mut_Change_Collection &left_mut, int &left_lower_bound,
    Bounded_Mut_Change_Collection &right_mut, int &right_lower_bound,
    // Go to ancestor
    Bounded_Mut_Change_Collection &mutations_out) {
    mutations_out.push_back(in);
    split_children(parent_node, radius_left, left_mut, left_lower_bound, right_mut,
              right_lower_bound, mutations_out.back());
}
static void add_mut_not_parent(
    // IN
    const Mutation_Count_Change_W_Lower_Bound &in, const MAT::Node *parent_node,
    int radius_left,
    // OUT
    // Go to sibling
    Bounded_Mut_Change_Collection &left_mut, int &left_lower_bound,
    Bounded_Mut_Change_Collection &right_mut, int &right_lower_bound) {
    Mutation_Count_Change_W_Lower_Bound temp(in);
    split_children(parent_node, radius_left, left_mut, left_lower_bound, right_mut,
              right_lower_bound, temp);
}

void search_subtree_first_level(MAT::Node *node, MAT::Node *to_exclude,
                                const src_side_info &src_side, int radius_left,
                                const Bounded_Mut_Change_Collection &left,
                                int left_lower_bound,
                                const Bounded_Mut_Change_Collection &right,
                                int right_lower_bound,int split_lower_bound) {
    for (auto child : node->children) {
        if (child == to_exclude) {
            continue;
        }
        if (child->dfs_end_index < to_exclude->dfs_index) {
#ifndef CHECK_BOUND
            if (use_bound&&left_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, left,
                                   split_lower_bound,left_lower_bound);
        } else {
            assert(child->dfs_index > to_exclude->dfs_end_index);
#ifndef CHECK_BOUND
            if (use_bound&&right_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, right,
                                   split_lower_bound,right_lower_bound);
        }
    }
}
void search_subtree_first_level(MAT::Node *node, MAT::Node *to_exclude,
                                const src_side_info &src_side, int radius_left,
                                const Bounded_Mut_Change_Collection &either) {
    for (auto child : node->children) {
        if (child == to_exclude) {
            continue;
        }
        if (child->dfs_end_index < to_exclude->dfs_index) {
#ifndef CHECK_BOUND
            if (use_bound&&left_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, either,0,
                                   0);
        } else {
            assert(child->dfs_index > to_exclude->dfs_end_index);
#ifndef CHECK_BOUND
            if (use_bound&&right_lower_bound > src_side.out.score_change) {
                continue;
            }
#endif
            search_subtree_bounded(child, src_side, radius_left - 1, either,0,
                                   0);
        }
    }
}
static void allele_count_reserve(
    Mutation_Count_Change_Collection &allele_change_out, size_t par_cnt,
    const Mutation_Count_Change_Collection &src_branch_change) {
    allele_change_out.reserve(par_cnt + 1);
}
static void
allele_count_reserve(Mutation_Count_Change_Collection &allele_change_out,
                     size_t par_cnt,
                     const MAT::Mutations_Collection &src_branch_change) {
    allele_change_out.reserve(std::min(par_cnt, src_branch_change.size()) + 1);
}
static nuc_one_hot allele_cnt_change_middle(
    Mutation_Count_Change_Collection::const_iterator
        &src_allele_cnt_change_iter,
    Mutation_Count_Change_Collection::const_iterator &src_allele_cnt_change_end,
    const MAT::Mutation &mut,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &next_src_par_score_lower_bound_diff) {
    nuc_one_hot major_allele = mut.get_all_major_allele();
    while (src_allele_cnt_change_iter->get_position() < mut.get_position()) {
        auto change = src_allele_cnt_change_iter->get_default_change_internal();
        next_src_par_score += change;
        src_allele_cnt_change_iter++;
    }
    if (src_allele_cnt_change_iter->get_position() == mut.get_position()) {
        int score_change = 0;
        major_allele = decrement_increment_mutation_count(
            mut, (Mutation_Count_Change)*src_allele_cnt_change_iter,
            allele_change_out, score_change);
        next_src_par_score += score_change;
        // this node incrementing some allele that can reduce parsimony
        // score among ancestors of  parent node
        if (score_change >= 0 && src_allele_cnt_change_iter->get_incremented() &
                                     mut.get_sensitive_increment()) {
            next_src_par_score_lower_bound_diff--;
        }
        src_allele_cnt_change_iter++;
    }
    return major_allele;
}
static void allele_cnt_change_end(
    Mutation_Count_Change_Collection::const_iterator
        &src_allele_cnt_change_iter,
    Mutation_Count_Change_Collection::const_iterator &src_allele_cnt_change_end,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &next_src_par_score_lower_bound_diff) {
    while (src_allele_cnt_change_iter != src_allele_cnt_change_end) {
        auto change = src_allele_cnt_change_iter->get_default_change_internal();
        next_src_par_score += change;
        src_allele_cnt_change_iter++;
    }
}

static nuc_one_hot allele_cnt_change_middle(
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_iter,
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_end,
    const MAT::Mutation &mut,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &next_src_par_score_lower_bound_diff) {
    nuc_one_hot major_allele = mut.get_all_major_allele();
    while (src_allele_cnt_change_iter!=src_allele_cnt_change_end&&src_allele_cnt_change_iter->get_position() < mut.get_position()) {
        if (src_allele_cnt_change_iter->is_valid()) {
            next_src_par_score--;
        }
        src_allele_cnt_change_iter++;
    }
    if (src_allele_cnt_change_iter!=src_allele_cnt_change_end&&src_allele_cnt_change_iter->get_position() == mut.get_position()) {
        int score_change = -1;
        major_allele = decrement_mutation_count(
            allele_change_out, mut,
            src_allele_cnt_change_iter->get_all_major_allele(), score_change);
        // this node incrementing some allele that can reduce parsimony
        // score among ancestors of  parent node
        next_src_par_score += score_change;
        if (!(src_allele_cnt_change_iter->get_all_major_allele() &
              mut.get_sensitive_decrement())) {
            next_src_par_score_lower_bound_diff--;
        }
        src_allele_cnt_change_iter++;
    } else {
        int score_change = -1;
        major_allele = decrement_mutation_count(
            allele_change_out, mut, mut.get_mut_one_hot(), score_change);
        // this node incrementing some allele that can reduce parsimony
        // score among ancestors of  parent node
        next_src_par_score += score_change;
        if (!(mut.get_par_one_hot() & mut.get_sensitive_decrement())) {
            next_src_par_score_lower_bound_diff--;
        }
    }
    return major_allele;
}
static void allele_cnt_change_end(
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_iter,
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_end,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &next_src_par_score_lower_bound_diff) {
    while (src_allele_cnt_change_iter != src_allele_cnt_change_end) {
        if (src_allele_cnt_change_iter->is_valid()) {
            next_src_par_score--;
        }
        src_allele_cnt_change_iter++;
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
template <typename T>
static bool
upward_integrated(src_side_info &src_side,
                  Bounded_Mut_Change_Collection &src_mut_in, int radius_left,
                  Bounded_Mut_Change_Collection &mut_out, const T &src_branch) {
                      auto old_LCA=src_side.LCA;
    MAT::Node *node = old_LCA->parent;
    if (!node) {
        return false;
    }
    if (node->dfs_index==2330) {
        //fputc('a', stderr);
    }
    // IN: alelle cnt change from this node (that will change the major allele
    // at parent node)
    auto src_allele_cnt_change_iter = src_branch.begin();
    auto src_allele_cnt_change_end = src_branch.end();

    // IN: Mutations on src node if placed as children of this node to keep its
    // and its descendant unchanged
    Bounded_Mut_Iter src_mut_iter = src_mut_in.begin();

    // OUT: Mutations on src node if placed as children of parent node, have two
    // copies, with same mutations, but different range Left side: nodes whose
    // DFS index is less than this node
    auto max_mut_size = node->mutations.size() + src_mut_in.size();
    Bounded_Mut_Change_Collection left_mut;
    int left_lower_bound = 0;
    // Right side: nodes whose DFS idx > all descendants of this node
    Bounded_Mut_Change_Collection right_mut;
    int right_lower_bound = 0;
    // mut_out is for going up to parent, with range covering the whole subtree
    // of parent node
    if (use_bound) {
    left_mut.reserve(max_mut_size);
    right_mut.reserve(max_mut_size);
    }
    mut_out.reserve(max_mut_size);

    // allele cnt change from parent node
    Mutation_Count_Change_Collection allele_change_out;
    allele_count_reserve(allele_change_out, node->mutations.size(), src_branch);
    // par socre mchange from src node removal on src->parent node (LCA) branch
    int next_src_par_score = src_side.par_score_change_from_src_remove;
    // Less increments from old_LCA that can improve parsimony at ancestors of
    // parent node as lower bound
    int next_src_par_score_lower_bound_diff = 0;
    int next_split_par_score=0;

    // Only include increment form adding src between parent node and its
    // parent, par score change on the src branch node is next_src_par_score
    // calculated in the same loop
    int par_score_change_split_LCA = 0;
    Mutation_Count_Change_Collection split_allele_count_change_out;

    for (const auto &mut : node->mutations) {
        if (mut.get_position()==3106) {
        fputc('a', stderr);
        }
        // Calculate allele count change for parent (and new major allele at
        // parent node for splitting between parent node and its parent )
        nuc_one_hot major_allele = allele_cnt_change_middle(
            src_allele_cnt_change_iter, src_allele_cnt_change_end, mut,
            allele_change_out, next_src_par_score,
            next_src_par_score_lower_bound_diff);
        // Get mutations if placed as children of parent node, and whether it is
        // profitable to place between parent node and parent of parent node
        while (src_mut_iter->get_position() < mut.get_position()) {
            add_mut(*src_mut_iter, node, radius_left, left_mut,
                    left_lower_bound, right_mut, right_lower_bound, mut_out);
            if (!(src_mut_iter->get_incremented()&src_mut_iter->get_par_state())) {
                next_split_par_score++;                
            }
            add_node_split(*src_mut_iter, split_allele_count_change_out,
                           par_score_change_split_LCA);
            src_mut_iter++;
        }
        if (src_mut_iter->get_position() == mut.get_position()) {
            // accumulate mutation for placement
            auto new_par_nuc = mut.get_par_one_hot();
            if (new_par_nuc != src_mut_iter->get_incremented()) {
                add_mut(*src_mut_iter, node, radius_left,
                        left_mut, left_lower_bound, right_mut,
                        right_lower_bound, mut_out);
                mut_out.back().set_par_nuc(new_par_nuc);
            }else {
                add_mut_not_parent(*src_mut_iter, node, radius_left,
                       left_mut, left_lower_bound, right_mut,
                        right_lower_bound);
            }
            if (!(mut.get_sensitive_increment() &
                            src_mut_iter->get_incremented())) {
                                next_split_par_score++;
            }
            // split LCA
            add_node_split(mut, src_mut_iter->get_incremented(), major_allele,
                           split_allele_count_change_out,
                           par_score_change_split_LCA);
            src_mut_iter++;
        } else {
            if (mut.is_valid()) {
                mut_out.emplace_back(mut, 0, mut.get_mut_one_hot());
                // get range
                if (use_bound) {
                    mut_out.back().init(node);                    
                }
                //next_split_par_score++;
            }
            add_node_split(mut,mut.get_all_major_allele(),major_allele,mut.get_mut_one_hot(), split_allele_count_change_out,
                           par_score_change_split_LCA);
        }
    }
    // Finish off allele count change
    allele_cnt_change_end(src_allele_cnt_change_iter, src_allele_cnt_change_end,
                          allele_change_out, next_src_par_score,
                          next_src_par_score_lower_bound_diff);
    // Add sentinel
    allele_change_out.emplace_back();
    while (src_mut_iter->get_position()!=INT_MAX) {
        add_mut(*src_mut_iter, node, radius_left, left_mut,
                left_lower_bound, right_mut, right_lower_bound, mut_out);
        if (!(src_mut_iter->get_par_state()&src_mut_iter->get_incremented())) {
            next_split_par_score++;            
        }
        add_node_split(*src_mut_iter, split_allele_count_change_out,
                       par_score_change_split_LCA);
        src_mut_iter++;
    }
    mut_out.emplace_back();
    left_mut.emplace_back();
    right_mut.emplace_back();
    src_side.LCA = node;
    src_side.src_par_score_lower_bound =
        next_src_par_score + next_src_par_score_lower_bound_diff;
    left_lower_bound+=src_side.src_par_score_lower_bound;
    right_lower_bound+=src_side.src_par_score_lower_bound;
    next_split_par_score+=src_side.src_par_score_lower_bound;
    if (use_bound) {
    search_subtree_first_level(node,old_LCA, src_side, radius_left,
                               left_mut, left_lower_bound, right_mut,
                               right_lower_bound,next_split_par_score);

    }else {
    search_subtree_first_level(node, old_LCA, src_side, radius_left,
                               src_mut_in);

    }
    
    src_side.par_score_change_from_src_remove = next_src_par_score;
    src_side.allele_count_change_from_src = std::move(allele_change_out);
    src_side.node_stack_from_src.push_back(node);

    output_LCA(split_allele_count_change_out, par_score_change_split_LCA+next_src_par_score,
               src_side, radius_left);
    return true;
}
void find_moves_bounded(MAT::Node* src,output_t& out,int search_radius){
    Bounded_Mut_Change_Collection src_mut;
    Bounded_Mut_Change_Collection src_mut_next;
    src_mut.reserve(src->mutations.size());
    src_side_info src_side{out,0,0,src,src};
    for (const auto& mut : src->mutations) {
        src_mut.emplace_back(mut,0,mut.get_all_major_allele());
        if (use_bound) {
            src_mut.back().init(src);    
        }
    }
    //sentinel
    src_mut.emplace_back();
    upward_integrated(src_side, src_mut, search_radius, src_mut_next, src->mutations);
    for (int radius_left=search_radius-1; radius_left>0; radius_left--) {
        src_mut.swap(src_mut_next);
        src_mut_next.clear();
        upward_integrated(src_side, src_mut, radius_left, src_mut_next, src_side.allele_count_change_from_src);
    }
    
}