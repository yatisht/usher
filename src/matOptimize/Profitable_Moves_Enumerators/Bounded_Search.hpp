#include "Profitable_Moves_Enumerators.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <utility>
#include <vector>
struct src_side_info {
    output_t &out;
#ifdef CHECK_BOUND
    counters &savings;
#endif
    int par_score_change_from_src_remove;
    int src_par_score_lower_bound;
    MAT::Node *LCA;
    MAT::Node *src;
    Mutation_Count_Change_Collection allele_count_change_from_src;
    std::vector<MAT::Node *> node_stack_from_src;
};
class Mutation_Count_Change_W_Lower_Bound : public Mutation_Count_Change {
    static const uint8_t LEVEL_END = -1;
    uint8_t par_sensitive_increment;
    uint8_t next_level;
    uint16_t idx;
    void init(const MAT::Node *node) {
        assert(use_bound);
        const auto& this_aux=addable_idxes[get_position()];
        idx =this_aux.find_idx(node);
        if (idx != EMPTY_POS) {
            next_level =this_aux.nodes[idx].level;
        } else {
            next_level = LEVEL_END;
        }
    }
    void to_descendant_adjust_range(Mutation_Count_Change_W_Lower_Bound &in,
                                    const MAT::Node *node, int level_left) {
        if (in.next_level > node->level) {
            return;
        }
        const auto &addable_idxes_this_pos = addable_idxes[get_position()];
        auto end_idx = node->dfs_end_index;
        // skip through irrelevent vertexes (shouldn't happen)
        for (;
             addable_idxes_this_pos.nodes[in.idx].dfs_end_idx < node->dfs_index;
             in.idx++) {
            fputc('a', stderr);
        }
        if (addable_idxes_this_pos.nodes[in.idx].dfs_start_idx > end_idx) {
            in.next_level = LEVEL_END;
            idx = in.idx;
            level_left=LEVEL_END;
            return;
        }
        next_level = level_left;
        // accumulate closest level
        for (; addable_idxes_this_pos.nodes[in.idx].dfs_start_idx < end_idx;
             in.idx++) {
            const auto &addable = addable_idxes_this_pos.nodes[in.idx];
            assert(addable.dfs_end_idx < end_idx);
            if (test_level(level_left, addable)) {
                if (addable.level < next_level) {
                    next_level = addable.level;
                }
            }
        }
    }
    // going to descendant
    bool test_level(int radius_left, const range_tree_node &node) {
        for (int idx = 0; idx < 4; idx++) {
            if (get_incremented() & (1 << idx)) {
                if (node.min_level[idx] < radius_left) {
                    return true;
                }
            }
        }
        return false;
    }

  public:
    struct to_ancestor {};
    struct to_descendent {};
    struct src_node {};
    uint8_t get_senesitive_increment()const{
        return par_sensitive_increment;
    }
    void set_sensitive_increment(uint8_t in){
        par_sensitive_increment=in;
    }

    void to_ancestor_adjust_range(const MAT::Node *node) {
        const auto &addable_idxes_this_pos = addable_idxes[get_position()];
        auto par_idx = addable_idxes_this_pos.nodes[idx].parent_idx;
        uint32_t idx = 0;
        while (par_idx != EMPTY_POS &&
               addable_idxes_this_pos.nodes[par_idx].dfs_start_idx >
                   node->dfs_index) {
            par_idx = addable_idxes_this_pos.nodes[par_idx].parent_idx;
        }
        if (par_idx != EMPTY_POS) {
            idx = addable_idxes_this_pos.nodes[par_idx].children_start_idx;
        }
        for (; addable_idxes_this_pos.start_idxes[idx] < node->dfs_index;
             idx++) {
        }
        next_level = node->level;
    }
    // Going down no coincide
    Mutation_Count_Change_W_Lower_Bound(){}
    Mutation_Count_Change_W_Lower_Bound(Mutation_Count_Change_W_Lower_Bound &in,
                                        const MAT::Node *node, int level_left,
                                        to_descendent)
        : Mutation_Count_Change_W_Lower_Bound(in) {
        to_descendant_adjust_range(in, node, level_left);
        set_sensitive_increment(in.get_par_state());
    }
    // Going Down coincide
    Mutation_Count_Change_W_Lower_Bound(Mutation_Count_Change_W_Lower_Bound &in,
                                        const MAT::Node *node, int level_left,
                                        const MAT::Mutation &coincided_mut,
                                        to_descendent)
        : Mutation_Count_Change_W_Lower_Bound(in) {
        to_descendant_adjust_range(in, node, level_left);
        set_sensitive_increment( coincided_mut.get_sensitive_increment() &
                                  (coincided_mut.get_all_major_allele() |
                                   coincided_mut.get_boundary1_one_hot()));
        assert(get_senesitive_increment() != coincided_mut.get_mut_one_hot());
    }
    /*
    // Going Up not Coincide
    Mutation_Count_Change_W_Lower_Bound(const Mutation_Count_Change_W_Lower_Bound &in,
                                        const MAT::Node *node, to_ancestor)
        : Mutation_Count_Change_W_Lower_Bound(in) {
        to_ancestor_adjust_range(in, node);
        par_sensitive_increment = in.get_par_state();
    }
    // Going up coincide
    Mutation_Count_Change_W_Lower_Bound(const Mutation_Count_Change_W_Lower_Bound &in,
                                        const MAT::Node *node,
                                        const MAT::Mutation &coincided_mut,
                                        to_ancestor)
        : Mutation_Count_Change_W_Lower_Bound(in) {
        to_ancestor_adjust_range(in, node);
        par_sensitive_increment = coincided_mut.get_sensitive_increment() &
                                  (coincided_mut.get_all_major_allele() |
                                   coincided_mut.get_boundary1_one_hot());
        assert(par_sensitive_increment != coincided_mut.get_mut_one_hot());
    }
    */
    // Going down new
    Mutation_Count_Change_W_Lower_Bound(const MAT::Mutation &in,
                                        const MAT::Node *node, int level_left,
                                        to_descendent)
        : Mutation_Count_Change(in, 0, in.get_par_one_hot()) {
        assert(in.is_valid());
        set_par_nuc(in.get_mut_one_hot());
        set_sensitive_increment(in.get_sensitive_increment() &
                                  (in.get_all_major_allele() |
                                   in.get_boundary1_one_hot()));
        init(node);
        if (next_level != LEVEL_END) {
            for (; addable_idxes[get_position()].nodes[idx].dfs_start_idx <
                   node->dfs_end_index;
                 idx++) {
                const auto &addable = addable_idxes[get_position()].nodes[idx];
                if (test_level(level_left, addable)) {
                    if (addable.level < next_level) {
                        next_level = addable.level;
                    }
                }
            }
        }
    }
    // Going up new
    Mutation_Count_Change_W_Lower_Bound(const MAT::Mutation &in,
                                        const MAT::Node *node,
                                        to_ancestor)
        : Mutation_Count_Change(in, 0, in.get_mut_one_hot()) {
        assert(in.is_valid());
        set_par_nuc(in.get_par_one_hot());
        set_sensitive_increment(in.get_sensitive_increment() &
                                  (in.get_all_major_allele() |
                                   in.get_boundary1_one_hot()));
        init(node);
    }
    // Going up src_node, no coincide
    Mutation_Count_Change_W_Lower_Bound(const MAT::Mutation &in,
                                        const MAT::Node *node, src_node)
        : Mutation_Count_Change(in, 0, in.get_all_major_allele()) {
        set_par_nuc(in.get_par_one_hot());
        set_sensitive_increment(in.get_par_one_hot());
        init(node);
    }
    /*
    // Going up src_node, no coincide
    Mutation_Count_Change_W_Lower_Bound(const MAT::Mutation &src_mut,
                                        const MAT::Node *node,
                                        const MAT::Mutation &coincided_mut,
                                        src_node)
        : Mutation_Count_Change(src_mut, 0, src_mut.get_all_major_allele()) {
        assert(coincided_mut.get_par_one_hot() !=
               src_mut.get_all_major_allele());
        set_par_nuc(coincided_mut.get_par_one_hot());
        par_sensitive_increment = coincided_mut.get_sensitive_increment() &
                                  (coincided_mut.get_all_major_allele() |
                                   coincided_mut.get_boundary1_one_hot());
        init(node);
    }*/
    bool valid_on_subtree()const{
        return next_level!=LEVEL_END||(get_incremented()&get_par_state());
    }
};
typedef std::vector<Mutation_Count_Change_W_Lower_Bound>
    Bounded_Mut_Change_Collection;
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                    int radius_left,
                    Bounded_Mut_Change_Collection &par_muts,
                    int lower_bound
#ifdef CHECK_BOUND
                                      ,bool do_continue,bool first_level
#endif
                    );