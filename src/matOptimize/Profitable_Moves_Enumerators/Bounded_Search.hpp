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
#define LEVEL_END UINT8_MAX
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
};
struct ignore_ranger_nop {
    ignore_ranger_nop(const Mutation_Annotated_Tree::ignored_t& in) {}
    ignore_ranger_nop() {}
    bool operator()(int) {
        return false;
    }
};
struct ignore_ranger {
    Mutation_Annotated_Tree::ignored_t::const_iterator iter;
    ignore_ranger(const Mutation_Annotated_Tree::ignored_t& in) {
        iter=in.begin();
    }
    bool operator()(int pos) {
        while (iter->second<pos) {
            iter++;
        }
        return iter->first<=pos&&iter->second>=pos;
    }
};
class Mutation_Count_Change_W_Lower_Bound_Downward : public Mutation_Count_Change {
    uint8_t par_sensitive_increment;
    uint8_t next_level;
    uint32_t idx;
    bool descend(const range_tree & addable_idxes_this_pos,uint8_t level,std::vector<int>& start_useful_idxes,std::vector<int>& end_useful_idxes,uint8_t max_level,uint32_t end_dfs_idx) {
        if (addable_idxes_this_pos.nodes[idx].level==level) {
            auto end_idx=idx;
            bool continue_find_end_idx=true;
            while (idx != EMPTY_POS&&addable_idxes_this_pos.nodes[idx].level==level) {
                for(; addable_idxes_this_pos.nodes[idx].dfs_end_idx<=end_dfs_idx; idx++) {
                    const auto& this_idx_tree_node=addable_idxes_this_pos.nodes[idx];
                    if (this_idx_tree_node.children_start_idx!=EMPTY_POS&&test_level(max_level,this_idx_tree_node)) {
                        goto FOUND;
                    }
                }
                next_level=LEVEL_END;
                return true;
FOUND:
                if (continue_find_end_idx) {
                    end_idx=std::max(end_idx,idx);
                    bool found_at_least_one=false;
                    for (uint32_t prode_idx=end_idx; addable_idxes_this_pos.nodes[prode_idx].dfs_end_idx<=end_dfs_idx; prode_idx++) {
                        const auto& this_idx_tree_node=addable_idxes_this_pos.nodes[prode_idx];
                        if (this_idx_tree_node.children_start_idx!=EMPTY_POS&&test_level(max_level,this_idx_tree_node)) {
                            end_idx=std::max(end_idx,prode_idx);
                            found_at_least_one=true;
                        }
                    }
                    if (found_at_least_one) {
                        end_idx=addable_idxes_this_pos.nodes[end_idx].children_start_idx;
                    } else {
                        continue_find_end_idx=false;
                    }
                }
                idx=addable_idxes_this_pos.nodes[idx].children_start_idx;

            }
            start_useful_idxes.push_back(addable_idxes_this_pos.nodes[idx].dfs_start_idx);
            end_useful_idxes.push_back(addable_idxes_this_pos.nodes[idx].dfs_end_idx);
            next_level=0;
            return true;
        }
        return false;
    }
    bool forward_useless_idx(const MAT::Node* node, int level_left) {
        for (; addable_idxes[get_position()].nodes[idx].dfs_start_idx <=
                node->dfs_end_index;
                idx++) {
            const auto &addable = addable_idxes[get_position()].nodes[idx];
            if (test_level(level_left+node->level, addable)) {
                return false;;
            }
        }
        next_level=LEVEL_END;
        return true;
    }

    void to_descendant_adjust_range(Mutation_Count_Change_W_Lower_Bound_Downward &in,
                                    const MAT::Node *node, int level_left,std::vector<int>& start_useful_idx,std::vector<int>& end_useful_idx) {
        if (!use_bound) {
            return;
        }
        /*if (get_position()==11782&&node->dfs_index==14445) {
            fputc('a', stderr);
        }
        if (get_position()==21724&&node->dfs_index==28242) {
            fputc('a', stderr);
        }*/
        if (in.next_level==LEVEL_END) {
            return;
        }
        next_level=LEVEL_END;
        const auto &addable_idxes_this_pos = addable_idxes[get_position()];
        auto end_idx = node->dfs_end_index;
        for (;
                addable_idxes_this_pos.end_idxes[in.idx] < node->dfs_index;
                in.idx++) {
        }
        if (addable_idxes_this_pos.nodes[in.idx].dfs_start_idx > end_idx) {
            idx = EMPTY_POS;
        }
        if (idx==EMPTY_POS) {
#ifdef CHECK_IDX
            auto test_idx=addable_idxes[get_position()].find_idx(node);
            if(test_idx!=EMPTY_POS&&addable_idxes[get_position()].nodes[test_idx].level>=node->level) {
                for(; addable_idxes[get_position()].start_idxes[test_idx]<=node->dfs_end_index; test_idx++) {
                    if (test_level(node->level+level_left, addable_idxes[get_position()].nodes[test_idx])) {
                        assert(addable_idxes[get_position()].nodes[test_idx].level>=node->level);
                    }
                }
            }
#endif
            return;
        }
        //assert(in.idx==addable_idxes_this_pos.find_idx(node));
        assert(in.idx==EMPTY_POS||in.idx<addable_idxes_this_pos.nodes.size());
        idx=in.idx;
        if(descend(addable_idxes_this_pos, node->level,start_useful_idx,end_useful_idx,level_left+node->level,node->dfs_end_index)) {
            return;
        }
        if (forward_useless_idx(node, level_left)) {
            return;
        }
        assert(addable_idxes_this_pos.nodes[in.idx].level>node->level);
        set_next_level(node, level_left,start_useful_idx,end_useful_idx,in.idx);
        assert(in.idx==EMPTY_POS||in.idx<addable_idxes_this_pos.nodes.size());
    }
    // going to descendant
    bool test_level(int radius_left, const range_tree_node &node) {
        for (int idx = 0; idx < 4; idx++) {
            if (get_incremented() & (1 << idx)) {
                if (node.min_level[idx] <= radius_left) {
                    return true;
                }
            }
        }
        return false;
    }
    void set_next_level(const MAT::Node* node, int level_left,std::vector<int>& start_useful_idxes,std::vector<int>& end_useful_idxes,uint32_t& end_idx) {
        next_level = LEVEL_END;
        bool found=false;
        uint32_t end_useful_idx=0;
        const auto& this_useful_idx=addable_idxes[get_position()];
        for (end_idx=idx; this_useful_idx.nodes[end_idx].dfs_start_idx <=
                node->dfs_end_index;
                end_idx++) {
            const auto &addable = this_useful_idx.nodes[end_idx];
            next_level=std::min((uint8_t)addable.level,next_level);
            if (test_level(level_left+node->level, addable)) {
                found=true;
                end_useful_idx=std::max(end_useful_idx,end_idx);
            }
        }
        if (!found) {
            next_level=LEVEL_END;
        } else {
            start_useful_idxes.push_back(this_useful_idx.nodes[idx].dfs_start_idx);
            end_useful_idxes.push_back(this_useful_idx.nodes[end_useful_idx].dfs_start_idx);
        }
    }
  public:
    uint8_t get_senesitive_increment()const {
        return par_sensitive_increment;
    }
    void set_sensitive_increment(uint8_t in) {
        par_sensitive_increment=in;
    }
    Mutation_Count_Change_W_Lower_Bound_Downward(
        Mutation_Count_Change& base,
        uint8_t par_sensitive_increment,uint8_t next_level,uint16_t idx)
        :Mutation_Count_Change(base),par_sensitive_increment(par_sensitive_increment)
        ,next_level(next_level),idx(idx) {}
    // Going down no coincide
    Mutation_Count_Change_W_Lower_Bound_Downward() {}
    Mutation_Count_Change_W_Lower_Bound_Downward(Mutation_Count_Change_W_Lower_Bound_Downward &in,
            const MAT::Node *node, int level_left,std::vector<int>& start_useful_idx,std::vector<int>& end_useful_idx)
        : Mutation_Count_Change_W_Lower_Bound_Downward(in) {
        to_descendant_adjust_range(in, node, level_left,start_useful_idx,end_useful_idx);
        set_sensitive_increment(in.get_par_state());
    }
    // Going Down coincide
    Mutation_Count_Change_W_Lower_Bound_Downward(Mutation_Count_Change_W_Lower_Bound_Downward &in,
            const MAT::Node *node, int level_left,
            const MAT::Mutation &coincided_mut,std::vector<int>& start_useful_idx,std::vector<int>& end_useful_idx)
        : Mutation_Count_Change_W_Lower_Bound_Downward(in) {
        to_descendant_adjust_range(in, node, level_left,start_useful_idx,end_useful_idx);
        set_par_nuc(coincided_mut.get_mut_one_hot());
        set_sensitive_increment( coincided_mut.get_sensitive_increment() &
                                 (coincided_mut.get_all_major_allele() |
                                  coincided_mut.get_boundary1_one_hot()));
        //assert(get_senesitive_increment() != coincided_mut.get_mut_one_hot());
    }

    // Going down new
    Mutation_Count_Change_W_Lower_Bound_Downward(const MAT::Mutation &in,
            const MAT::Node *node, int level_left,std::vector<int>& start_useful_idx,std::vector<int>& end_useful_idx)
        : Mutation_Count_Change(in, 0, in.get_par_one_hot()) {
        /*if (get_position()==12809&&node->dfs_index==32047) {
            fputc('a', stderr);
        }*/
        assert(in.is_valid());
        set_par_nuc(in.get_mut_one_hot());
        set_sensitive_increment(in.get_sensitive_increment() &
                                (in.get_all_major_allele() |
                                 in.get_boundary1_one_hot()));
        if (!use_bound) {
            return;
        }
        const auto& this_possition_addable_idx=addable_idxes[get_position()];
        auto start_idx=this_possition_addable_idx.find_idx(node);
        if (start_idx==EMPTY_POS) {
            next_level=LEVEL_END;
            return;
        }
        idx=start_idx;
        if(descend(addable_idxes[get_position()], node->level,start_useful_idx,end_useful_idx,node->level+level_left,node->dfs_end_index)) {
            return;
        }
        assert(this_possition_addable_idx.nodes[start_idx].level>node->level);
        if (forward_useless_idx(node, level_left)) {
            return;
        }
        uint32_t ignored;
        set_next_level(node, level_left,start_useful_idx,end_useful_idx,ignored);
    }
    bool valid_on_subtree()const {
        return next_level!=LEVEL_END||(get_incremented()&get_senesitive_increment());
    }
};
class Mutation_Count_Change_W_Lower_Bound_to_ancestor : public Mutation_Count_Change {
    uint32_t idx;
    void init(const MAT::Node *node) {
        if (!use_bound) {
            return;
        }
        const auto& this_aux=addable_idxes[get_position()];
        uint32_t next_idx=EMPTY_POS;
        /*if (get_position()==2061) {
            fputc('a', stderr);
        }*/
        idx =this_aux.find_idx(node,next_idx);
    }
    void output_absent_sibling(std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>& sibling_out) {
        sibling_out.emplace_back(*this,get_par_state(),LEVEL_END,EMPTY_POS);
    }
  public:
    void to_ancestor_adjust_range(const MAT::Node *node,std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>& sibling_out) {
        uint32_t start_idx=EMPTY_POS;
        if (!use_bound) {
            sibling_out.emplace_back(*this,get_par_state(),0,0);
            return;
        }
        /*if (get_position()==241&&node->dfs_index==73037) {
            fputc('a', stderr);
        }*/
        assert(get_incremented()!=get_par_state());
        const auto &addable_idxes_this_pos = addable_idxes[get_position()];
        if (idx==EMPTY_POS) {
            idx =addable_idxes_this_pos.find_idx(node);
            assert(idx!=EMPTY_POS||(get_incremented()&&get_par_state()));
            start_idx=idx;
        } else {
            while (idx!=EMPTY_POS&&addable_idxes_this_pos.nodes[idx].level >=
                    node->level) {
                idx = addable_idxes_this_pos.nodes[idx].parent_idx;
            }
            if (idx!=EMPTY_POS) {
                idx = addable_idxes_this_pos.nodes[idx].children_start_idx;
            }
            if (idx==EMPTY_POS) {
                idx=0;
            }
            for (; addable_idxes_this_pos.end_idxes[idx]<node->dfs_index; idx++ ) {}
            if (addable_idxes_this_pos.nodes[idx].dfs_start_idx>node->dfs_end_index) {
                idx=EMPTY_POS;
                output_absent_sibling(sibling_out);
                return;
            }
            start_idx=idx;
        }
        while (start_idx!=EMPTY_POS&&addable_idxes_this_pos.nodes[start_idx].level<=node->level) {
            /*for (; addable_idxes_this_pos.start_idxes[start_idx]<node->dfs_index;idx++ ) {}
                    if (addable_idxes_this_pos.start_idxes[idx]>node->dfs_end_index) {
            if (start_idx>0) {
                start_idx--;
            }
            if (addable_idxes_this_pos.nodes[start_idx].dfs_end_idx<node->dfs_index||
            addable_idxes_this_pos.start_idxes[start_idx]>node->dfs_end_index) {
                output_absent_sibling(sibling_out);
                return;
            }
            }*/
            start_idx=addable_idxes_this_pos.nodes[start_idx].children_start_idx;
        }
        /*uint32_t next_idx=start_idx;
        next_idx=addable_idxes_this_pos.nodes[start_idx].children_start_idx;
        while (next_idx!=EMPTY_POS&&addable_idxes_this_pos.nodes[next_idx].level==node->level+1) {
            start_idx=next_idx;
            next_idx=addable_idxes_this_pos.nodes[start_idx].children_start_idx;
        }*/
        assert(start_idx==EMPTY_POS||addable_idxes_this_pos.nodes[start_idx].dfs_start_idx>=node->dfs_index);
        if (start_idx == EMPTY_POS) {
            output_absent_sibling(sibling_out);
            return;
        } else {
            sibling_out.emplace_back(*this,get_par_state(),0,start_idx);
        }
    }
    void to_ancestor_adjust_range(const MAT::Node *node,std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>& sibling_out,const MAT::Mutation &coincided_mut) {
        to_ancestor_adjust_range(node,sibling_out);
        auto par_sensitive_increment = coincided_mut.get_sensitive_increment() &
                                       (coincided_mut.get_all_major_allele() |
                                        coincided_mut.get_boundary1_one_hot());
        sibling_out.back().set_sensitive_increment(par_sensitive_increment);
    }
    Mutation_Count_Change_W_Lower_Bound_to_ancestor() {}
    // Going Up not Coincide
    Mutation_Count_Change_W_Lower_Bound_to_ancestor(const Mutation_Count_Change_W_Lower_Bound_to_ancestor &in,
            const MAT::Node *node,std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>& sibling_out)
        : Mutation_Count_Change_W_Lower_Bound_to_ancestor(in) {
        assert(get_incremented()!=get_par_state());
        to_ancestor_adjust_range(node,sibling_out);
    }
    // Going up coincide
    Mutation_Count_Change_W_Lower_Bound_to_ancestor(const Mutation_Count_Change_W_Lower_Bound_to_ancestor &in,
            const MAT::Node *node,
            const MAT::Mutation &coincided_mut,std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>& sibling_out
                                                   )
        : Mutation_Count_Change_W_Lower_Bound_to_ancestor(in) {
        assert(get_incremented()!=get_par_state());
        to_ancestor_adjust_range(node,sibling_out);
        set_par_nuc(coincided_mut.get_par_one_hot());
        auto par_sensitive_increment = coincided_mut.get_sensitive_increment() &
                                       (coincided_mut.get_all_major_allele() |
                                        coincided_mut.get_boundary1_one_hot());
        sibling_out.back().set_sensitive_increment(par_sensitive_increment);
        //assert(par_sensitive_increment != coincided_mut.get_mut_one_hot());
    }
    // Going up new
    Mutation_Count_Change_W_Lower_Bound_to_ancestor(const MAT::Mutation &in,
            const MAT::Node *node)
        : Mutation_Count_Change(in, 0, in.get_mut_one_hot()) {
        assert(get_incremented()!=get_par_state());
        set_par_nuc(in.get_par_one_hot());
        init(node);
    }
    struct source_node {};
    Mutation_Count_Change_W_Lower_Bound_to_ancestor(const MAT::Mutation &in,
            const MAT::Node *node,source_node)
        : Mutation_Count_Change(in, 0, in.get_all_major_allele()) {
        assert(get_incremented()!=get_par_state());
        set_par_nuc(in.get_par_one_hot());
        init(node);
    }
};
typedef std::vector<Mutation_Count_Change_W_Lower_Bound_Downward>
Bounded_Mut_Change_Collection;
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                            int radius_left,
                            Bounded_Mut_Change_Collection &par_muts,
                            int lower_bound,ignore_ranger_nop,Reachable,Move_Found_Callback& callback
#ifdef CHECK_BOUND
                            ,bool do_continue
#endif
                           );
void search_subtree_bounded(MAT::Node *node, const src_side_info &src_side,
                            int radius_left,
                            Bounded_Mut_Change_Collection &par_muts,
                            int lower_bound,ignore_ranger,Reachable,Move_Found_Callback& callback
#ifdef CHECK_BOUND
                            ,bool do_continue
#endif
                           );