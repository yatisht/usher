#include "Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <utility>
#include <vector>
typedef const node_info* iter_t;
typedef std::pair<iter_t, iter_t> iter_range_t;
typedef std::array<iter_range_t, 4> range_per_allele_t;
struct src_side_info{
    output_t& out;
    int par_score_change_from_src_remove;
    int src_par_score_lower_bound;
    MAT::Node *LCA;
    MAT::Node *src;
    Mutation_Count_Change_Collection allele_count_change_from_src;
    std::vector<MAT::Node *> node_stack_from_src;
};
struct Mutation_Count_Change_W_Lower_Bound : public Mutation_Count_Change {
    range_per_allele_t ranges;
    bool have_content;
    uint8_t par_sensitive_increment;
    uint8_t par_sensitive_decrement;
    using Mutation_Count_Change::Mutation_Count_Change;
    void init(const MAT::Node *node) {
        assert(use_bound);
        have_content = true;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            auto start=all.data();
            auto end=all.data()+all.size();
            if (get_incremented() & (1 << nu_idx)) {
                auto start_iter = std::lower_bound(start,end,
                                                   node_info{node->dfs_index});
                auto end_iter = std::upper_bound(
                    start_iter, end, node_info{node->dfs_end_index});
                ranges[nu_idx] = std::make_pair(start_iter, end_iter);
            } else {
                ranges[nu_idx] = std::make_pair(end, end);
            }
            assert(ranges[nu_idx].first<=ranges[nu_idx].second);
        }
    }
    void set_iter(iter_t start, iter_t end, const MAT::Node *node,
                  size_t max_level, int nu_idx) {
        assert(use_bound);
        assert(start<=end);
        auto start_iter =
            std::lower_bound(start, end, node_info{node->dfs_index});
        ranges[nu_idx].first = end;
        ranges[nu_idx].second = end;
        while (start_iter != end && start_iter->dfs_idx <= node->dfs_end_index) {
            if (start_iter->level <= max_level) {
                ranges[nu_idx].first = start_iter;
                start_iter++;
                have_content = true;
                break;
            }
            start_iter++;
        }
        while (start_iter != end && start_iter->dfs_idx <= node->dfs_end_index) {
            if (start_iter->level <= max_level) {
                have_content = true;
                ranges[nu_idx].second = start_iter+1;
            }
            start_iter++;
        }
        assert(ranges[nu_idx].first<=ranges[nu_idx].second);
        assert(ranges[nu_idx].second<=end);
    }
    bool init(const MAT::Node *node, int radius_left) {
        assert(use_bound);
        int max_level = node->level + radius_left;
        have_content = false;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            auto end=all.data()+all.size();
            if (get_incremented() & (1 << nu_idx)) {
                set_iter(all.data(), end, node, max_level, nu_idx);
            } else {
                ranges[nu_idx] = std::make_pair(end, end);
            }
            assert(ranges[nu_idx].first<=ranges[nu_idx].second);
            
        }
        return have_content||(get_incremented()|get_par_state());
    }

    bool to_decendent(const MAT::Node *node, int radius_left) {
        assert(use_bound);
        if (!have_content) {
            return get_incremented()|get_par_state();
        }
        int max_level = node->level + radius_left;
        have_content = false;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            if (ranges[nu_idx].second - ranges[nu_idx].first) {
                set_iter(ranges[nu_idx].first, ranges[nu_idx].second, node,
                         max_level, nu_idx);
            }
            assert(ranges[nu_idx].first<=ranges[nu_idx].second);
        }
        return have_content||(get_incremented()|get_par_state());
    }
    std::tuple<Mutation_Count_Change_W_Lower_Bound,
               Mutation_Count_Change_W_Lower_Bound>
    to_parent(const MAT::Node *node, int radius_left) {
        assert(use_bound);
        size_t max_level = node->level + radius_left;
        std::tuple<Mutation_Count_Change_W_Lower_Bound,
                   Mutation_Count_Change_W_Lower_Bound>
            to_return{*this, *this};
        bool left_have_content=false;
        bool right_have_content=false;
        for (int nu_idx = 0; nu_idx < 4; nu_idx++) {
            const auto &all = addable_idxes[get_position()][nu_idx];
            auto end=all.data()+all.size();
            if (all.empty()) {
                continue;
            }
            // first half
            auto start_iter = ranges[nu_idx].first;
            auto start_in_range_iter = start_iter;
            auto start_to_set = start_iter;
            start_iter--;
            while (start_iter >= all.data() &&
                   start_iter->dfs_idx >= node->dfs_index) {
                if (start_iter->level <= max_level) {
                    left_have_content = true;
                    start_in_range_iter = start_iter;
                }
                start_to_set=start_iter;
                start_iter--;
            }
            std::get<0>(to_return).ranges[nu_idx] =
                std::make_pair(start_in_range_iter, ranges[nu_idx].first);
                assert(start_in_range_iter<=ranges[nu_idx].first);
            // second half
            auto end_iter = ranges[nu_idx].second;
            auto end_in_range_iter = end_iter;
            while (end_iter < end &&
                   end_iter->dfs_idx <= node->dfs_end_index) {
                if (end_iter->level <= max_level) {
                    right_have_content = true;
                    end_in_range_iter = end_iter+1;
                }
                end_iter++;
            }
            assert(end_in_range_iter<=end);
            std::get<1>(to_return).ranges[nu_idx] =
                std::make_pair(ranges[nu_idx].second,end_in_range_iter);
                assert(ranges[nu_idx].second<=end_in_range_iter);
            ranges[nu_idx] = std::make_pair(start_to_set, end_iter);
            assert(ranges[nu_idx].first<=ranges[nu_idx].second);
            std::get<0>(to_return).have_content=left_have_content;
            std::get<1>(to_return).have_content=right_have_content;
        }
        return to_return;
    }
};
typedef std::vector<Mutation_Count_Change_W_Lower_Bound>
    Bounded_Mut_Change_Collection;
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                    int radius_left,
                    const Bounded_Mut_Change_Collection &par_muts,
                    int lower_bound
#ifdef CHECK_BOUND
                                      ,bool do_continue
#endif
                    );