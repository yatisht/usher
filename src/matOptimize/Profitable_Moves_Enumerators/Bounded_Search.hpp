#include "Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <array>
#include <iterator>
#include <utility>
#include <vector>
typedef std::vector<node_info>::const_iterator iter_t;
typedef std::pair<iter_t, iter_t> iter_range_t;
typedef std::array<iter_range_t, 4> range_per_allele_t;
struct src_side_info{
    int par_score_change_from_src_remove;
    int src_par_score_lower_bound;
    MAT::Node *LCA;
    MAT::Node *src;
    Mutation_Count_Change_Collection allele_count_change_from_src;
    std::vector<MAT::Node *> node_stack_from_src;
    output_t& out;
};
struct Mutation_Count_Change_W_Lower_Bound : public Mutation_Count_Change {
    range_per_allele_t ranges;
    bool have_content;
    bool offsetable;
    using Mutation_Count_Change::Mutation_Count_Change;
    void init(const MAT::Node *node) {
        have_content = true;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            if (get_incremented() & (1 << nu_idx)) {
                auto start_iter = std::lower_bound(all.begin(), all.end(),
                                                   node_info{node->dfs_index});
                auto end_iter = std::lower_bound(
                    start_iter, all.end(), node_info{node->dfs_end_index});
                if (end_iter != all.end() && end_iter > start_iter &&
                    end_iter->dfs_idx > node->dfs_end_index) {
                    end_iter--;
                }
                ranges[nu_idx] = std::make_pair(start_iter, end_iter);
            } else {
                ranges[nu_idx] = std::make_pair(all.end(), all.end());
            }
        }
    }
    void set_iter(iter_t start, iter_t end, const MAT::Node *node,
                  int max_level, int nu_idx) {
        auto start_iter =
            std::lower_bound(start, end, node_info{node->dfs_index});
        while (start_iter != end && start_iter->dfs_idx < node->dfs_end_index) {
            if (start_iter->level <= max_level) {
                ranges[nu_idx].first = start_iter;
                start_iter++;
                have_content = true;
                break;
            }
            start_iter++;
        }
        while (start_iter != end && start_iter->dfs_idx < node->dfs_end_index) {
            if (start_iter->level <= max_level) {
                have_content = true;
                ranges[nu_idx].second = start_iter;
            }
            start_iter++;
        }
    }
    bool init(const MAT::Node *node, int radius_left) {
        int max_level = node->level + radius_left;
        have_content = false;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            if (get_incremented() & (1 << nu_idx)) {
                set_iter(all.begin(), all.end(), node, max_level, nu_idx);
            } else {
                ranges[nu_idx] = std::make_pair(all.end(), all.end());
            }
        }
        return have_content;
    }

    bool to_decendent(const MAT::Node *node, int radius_left) {
        if (!have_content) {
            return false;
        }
        int max_level = node->level + radius_left;
        have_content = false;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            if (ranges[nu_idx].second != ranges[nu_idx].first) {
                set_iter(ranges[nu_idx].first, ranges[nu_idx].second, node,
                         max_level, nu_idx);
            }
        }
        return have_content;
    }
    std::tuple<Mutation_Count_Change_W_Lower_Bound,
               Mutation_Count_Change_W_Lower_Bound, bool, bool>
    to_parent(const MAT::Node *node, int radius_left) {
        int max_level = node->level + radius_left;
        std::tuple<Mutation_Count_Change_W_Lower_Bound,
                   Mutation_Count_Change_W_Lower_Bound, bool, bool>
            to_return{*this, *this, false, false};
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            // first half
            auto start_iter = ranges[nu_idx].first;
            auto start_in_range_iter = ranges[nu_idx].first;
            while (start_iter >= all.begin() &&
                   start_iter->dfs_idx >= node->dfs_index) {
                if (start_iter->level <= max_level) {
                    std::get<2>(to_return) = true;
                    start_in_range_iter = start_iter;
                }
                start_iter--;
            }
            std::get<0>(to_return).ranges[nu_idx] =
                std::make_pair(start_in_range_iter, ranges[nu_idx].first);
            // second half
            auto end_iter = ranges[nu_idx].second;
            auto end_in_range_iter = ranges[nu_idx].second;
            while (end_iter < all.end() &&
                   end_iter->dfs_idx <= node->dfs_end_index) {
                if (end_iter->level <= max_level) {
                    std::get<3>(to_return) = true;
                    end_in_range_iter = end_iter;
                }
                end_iter++;
            }
            std::get<1>(to_return).ranges[nu_idx] =
                std::make_pair(end_in_range_iter, ranges[nu_idx].second);
            ranges[nu_idx] = std::make_pair(start_iter, end_iter);
        }
        return to_return;
    }
};
typedef std::vector<Mutation_Count_Change_W_Lower_Bound>
    Bounded_Mut_Change_Collection;
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                    int radius_left,
                    const Bounded_Mut_Change_Collection &par_muts,
                    int lower_bound);