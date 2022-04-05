#include "Bounded_Search.hpp"
#include "process_individual_mutation.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <signal.h>
typedef Bounded_Mut_Change_Collection::iterator Bounded_Mut_Iter;

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
    check_parsimony_score_change_above_LCA(
        actual_LCA, par_score_change, allele_count_change_from_splitting_LCA,
        parent_of_parent_added);
    std::vector<MAT::Node *> ignored;
#ifdef CHECK_PAR_MAIN
    output_t temp;
    auto ref_out=individual_move(src_side.src, actual_LCA, actual_LCA, temp);
    if (ref_out>=0) {
        if (par_score_change<0) {
            fprintf(stderr, "%s\n",src_side.src->identifier.c_str());
#ifdef STOP_ON_ERROR
            raise(SIGTRAP);
#endif
        }
        //assert(par_score_change>=0);
    } else {
        if (ref_out!=par_score_change) {
            fprintf(stderr, "%s\n",src_side.src->identifier.c_str());
#ifdef STOP_ON_ERROR
            raise(SIGTRAP);
#endif
        }
        //assert(ref_out==par_score_change);
    }
#endif
    output_result(src_side.src, actual_LCA, actual_LCA, par_score_change,
                  src_side.out, radius_left);
}
template<typename T>
static void search_subtree_first_level(MAT::Node *node, MAT::Node *to_exclude,
                                       const src_side_info &src_side, int radius_left,
                                       Bounded_Mut_Change_Collection &either, T& ignored_range,Reachable reachable) {
    for (auto child : node->children) {
        if (child == to_exclude) {
            continue;
        }
        auto reachable_child=reachable;
        if (!reachable.always_search) {
            if (!child->get_descendent_changed()) {
                continue;
            }
            if (reachable.reachable_change&&child->get_self_changed()) {
                reachable_child.always_search=true;
            }
        }
        search_subtree_bounded(child, src_side, radius_left - 1, either,
                               src_side.src_par_score_lower_bound,ignored_range,reachable_child
#ifdef CHECK_BOUND
                               ,true
#endif
                              );
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
    int &next_src_par_score,int&) {
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
        src_allele_cnt_change_iter++;
    }
    return major_allele;
}
static void allele_cnt_change_end(
    Mutation_Count_Change_Collection::const_iterator
    &src_allele_cnt_change_iter,
    Mutation_Count_Change_Collection::const_iterator &src_allele_cnt_change_end,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &) {
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
    int &next_src_par_score, int &src_side_lower_bound) {
    nuc_one_hot major_allele = mut.get_all_major_allele();
    while (src_allele_cnt_change_iter!=src_allele_cnt_change_end&&src_allele_cnt_change_iter->get_position() < mut.get_position()) {
        if (src_allele_cnt_change_iter->is_valid()) {
            next_src_par_score--;
            src_side_lower_bound--;
        }
        src_allele_cnt_change_iter++;
    }
    if (src_allele_cnt_change_iter!=src_allele_cnt_change_end&&src_allele_cnt_change_iter->get_position() == mut.get_position()) {
        assert(src_allele_cnt_change_iter->get_all_major_allele()!=0xf);
        int score_change = -1;
        major_allele = decrement_mutation_count(
                           allele_change_out, mut,
                           src_allele_cnt_change_iter->get_all_major_allele(), score_change);
        // this node incrementing some allele that can reduce parsimony
        // score among ancestors of  parent node
        next_src_par_score += score_change;
        if (!(src_allele_cnt_change_iter->get_all_major_allele() &
                mut.get_sensitive_decrement())) {
            src_side_lower_bound--;
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
            src_side_lower_bound--;
        }
    }
    return major_allele;
}
static void allele_cnt_change_end(
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_iter,
    MAT::Mutations_Collection::const_iterator &src_allele_cnt_change_end,
    Mutation_Count_Change_Collection &allele_change_out,
    int &next_src_par_score, int &src_side_lower_bound) {
    while (src_allele_cnt_change_iter != src_allele_cnt_change_end) {
        if (src_allele_cnt_change_iter->is_valid()) {
            next_src_par_score--;
            src_side_lower_bound--;
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
MAT::Mutation* find_mut(MAT::Mutations_Collection& muts, int pos) {
    auto iter=muts.find(pos);
    if (iter==muts.end()) {
        return nullptr;
    } else {
        return &(*iter);
    }
}
Mutation_Count_Change* find_mut(Mutation_Count_Change_Collection& muts, int pos) {
    for (auto &m:muts ) {
        if (m.get_position()==pos) {
            return &m;
        }
    }
    return nullptr;
}
Mutation_Count_Change_W_Lower_Bound_Downward* find_mut(Bounded_Mut_Change_Collection& muts, int pos) {
    for (auto &m:muts ) {
        if (m.get_position()==pos) {
            return &m;
        }
    }
    return nullptr;
}
static void src_mut_no_match(
    Mutation_Count_Change_W_Lower_Bound_to_ancestor &src_mut, const MAT::Node *node,
    std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> &mut_out, int &par_score_change_split_LCA,
    Mutation_Count_Change_Collection &split_allele_count_change_out,Bounded_Mut_Change_Collection& sibling_out) {
    mut_out.emplace_back(src_mut,node,sibling_out);
    add_node_split(src_mut, split_allele_count_change_out,
                   par_score_change_split_LCA);
}
template <typename T,typename S>
static bool
upward_integrated(src_side_info &src_side,
                  std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> &src_mut_in, int radius_left,
                  std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> &mut_out, const T &src_branch, S ignore_iter,Reachable reachable) {
    auto old_LCA=src_side.LCA;
    MAT::Node *node = old_LCA->parent;
    if (!node) {
        return false;
    }
    // IN: alelle cnt change from this node (that will change the major allele
    // at parent node)
    auto src_allele_cnt_change_iter = src_branch.begin();
    auto src_allele_cnt_change_end = src_branch.end();

    // IN: Mutations on src node if placed as children of this node to keep its
    // and its descendant unchanged
    auto src_mut_iter = src_mut_in.begin();

    // OUT: Mutations on src node if placed as children of parent node, have two
    // copies, with same mutations, but different range Left side: nodes whose
    // DFS index is less than this node
    auto max_mut_size = node->mutations.size() + src_mut_in.size();
    mut_out.reserve(max_mut_size);
    Bounded_Mut_Change_Collection sibling_muts;
    sibling_muts.reserve(src_mut_in.size());
    // allele cnt change from parent node
    Mutation_Count_Change_Collection allele_change_out;
    allele_count_reserve(allele_change_out, node->mutations.size(), src_branch);
    // par socre mchange from src node removal on src->parent node (LCA) branch
    int next_src_par_score = src_side.par_score_change_from_src_remove;
    // Less increments from old_LCA that can improve parsimony at ancestors of
    // parent node as lower bound

    // Only include increment form adding src between parent node and its
    // parent, par score change on the src branch node is next_src_par_score
    // calculated in the same loop
    int par_score_change_split_LCA = 0;
    Mutation_Count_Change_Collection split_allele_count_change_out;
    for (const auto &mut : node->mutations) {
        if (ignore_iter(mut.get_position())) {
            //in ignored range
            continue;
        }
        // Calculate allele count change for parent (and new major allele at
        // parent node for splitting between parent node and its parent )
        nuc_one_hot major_allele = allele_cnt_change_middle(
                                       src_allele_cnt_change_iter, src_allele_cnt_change_end, mut,
                                       allele_change_out, next_src_par_score,
                                       src_side.src_par_score_lower_bound);
        // Get mutations if placed as children of parent node, and whether it is
        // profitable to place between parent node and parent of parent node
        while (src_mut_iter->get_position() < mut.get_position()) {
            src_mut_no_match(*src_mut_iter, node, mut_out, par_score_change_split_LCA, split_allele_count_change_out,sibling_muts);
            src_mut_iter++;
        }
        if (src_mut_iter->get_position() == mut.get_position()) {
            // accumulate mutation for placement
            if ( mut.get_par_one_hot() != src_mut_iter->get_incremented()) {
                mut_out.emplace_back(*src_mut_iter,node,mut,sibling_muts);
            } else {
                src_mut_iter->to_ancestor_adjust_range(node, sibling_muts,mut);
            }
            // split LCA
            add_node_split(mut, src_mut_iter->get_incremented(), major_allele,
                           split_allele_count_change_out,
                           par_score_change_split_LCA);
            src_mut_iter++;
        } else {
            if (mut.is_valid()) {
                mut_out.emplace_back(mut, node);
            }
            add_node_split(mut,mut.get_all_major_allele(),major_allele,mut.get_mut_one_hot(), split_allele_count_change_out,
                           par_score_change_split_LCA);
        }
    }
    // Finish off allele count change
    allele_cnt_change_end(src_allele_cnt_change_iter, src_allele_cnt_change_end,
                          allele_change_out, next_src_par_score,
                          src_side.src_par_score_lower_bound);
    // Add sentinel
    allele_change_out.emplace_back();
    while (src_mut_iter->get_position()!=INT_MAX) {
        src_mut_no_match(*src_mut_iter, node, mut_out, par_score_change_split_LCA, split_allele_count_change_out,sibling_muts);
        src_mut_iter++;
    }
    mut_out.emplace_back();
    src_side.LCA = node;
    sibling_muts.emplace_back();
    search_subtree_first_level(node, old_LCA, src_side, radius_left,
                               sibling_muts,ignore_iter,reachable);



    src_side.par_score_change_from_src_remove = next_src_par_score;
    src_side.allele_count_change_from_src = std::move(allele_change_out);

    output_LCA(split_allele_count_change_out, par_score_change_split_LCA+next_src_par_score,
               src_side, radius_left);
    return true;
}
template<typename T>
void __find_moves_bounded(MAT::Node *&src, int &search_radius,
                          std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> &src_mut,
                          std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> &src_mut_next,
                          src_side_info &src_side, T ignored_pos,Reachable reachable) {
    upward_integrated(src_side, src_mut, search_radius, src_mut_next,
                      src->mutations,ignored_pos,reachable);
    for (int radius_left = search_radius - 1; radius_left > 0; radius_left--) {
        src_mut.swap(src_mut_next);
        src_mut_next.clear();
        if (!src_side.LCA) {
            break;
        }
        if (!reachable.always_search) {
            if (src_side.LCA->get_self_changed()&&reachable.reachable_change) {
                reachable.always_search=true;
            }
            if(!src_side.LCA->get_ancestor_changed()) {
                break;
            }
        }
        upward_integrated(src_side, src_mut, radius_left, src_mut_next,
                          src_side.allele_count_change_from_src,ignored_pos,reachable);
    }
}
void find_moves_bounded(MAT::Node *src, output_t &out, int search_radius,bool do_drift,Reachable reachable
#ifdef CHECK_BOUND
                        ,
                        counters &count
#endif
                       ) {
    int par_score_change_base=do_drift?-1:0;
    std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> src_mut;
    std::vector<Mutation_Count_Change_W_Lower_Bound_to_ancestor> src_mut_next;
    src_mut.reserve(src->mutations.size());
#ifdef CHECK_BOUND
    src_side_info src_side {out,count,0,0,src,src};
#else
    src_side_info src_side {out,par_score_change_base,par_score_change_base,src,src};
#endif
    if (reachable.reachable_change&&src->get_self_moved()) {
        reachable.always_search=true;
    }
    for (const auto& mut : src->mutations) {
        if (mut.get_all_major_allele()==0xf||mut.get_par_one_hot()==mut.get_all_major_allele()) {
            continue;
        }
        src_mut.emplace_back(mut,src,Mutation_Count_Change_W_Lower_Bound_to_ancestor::source_node());
    }
    //sentinel
    src_mut.emplace_back();
    if (src->ignore.empty()) {
        __find_moves_bounded(src, search_radius, src_mut, src_mut_next, src_side,ignore_ranger_nop{},reachable);
    } else {
        __find_moves_bounded(src, search_radius, src_mut, src_mut_next, src_side,ignore_ranger{src->ignore},reachable);
    }
}