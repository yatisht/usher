#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <iterator>
#include <linux/limits.h>
#include <stack>
#include <tbb/concurrent_vector.h>
#include <utility>
#include <iterator>
#include <vector>
#define MAX_DIST 16

struct make_tree_temp{
    uint32_t out_start_idx;
    uint32_t out_end_idx;
    uint32_t root_dfs_idx;
};
struct range_tree_temp{
    uint32_t dfs_start_idx;
    uint32_t dfs_end_idx;
    std::array<uint32_t, 4> min_level;
    std::vector<range_tree_temp*> children;
    uint32_t level;
    range_tree_temp(){}
    void init(MAT::Node* node){
        level=node->level;
        for (int i=0; i<4; i++) {
            min_level[i]=INT_MAX;
        }
    }

    range_tree_temp(MAT::Node* node){
        dfs_start_idx=node->dfs_index;
        dfs_end_idx=node->dfs_end_index;
        init(node);
    }

    range_tree_temp(MAT::Node* node,
        std::vector<range_tree_temp*>::const_iterator start,
        std::vector<range_tree_temp*>::const_iterator end): children(std::make_reverse_iterator(end),std::make_reverse_iterator(start)){
            init(node);
            assert(end-start>1);
            dfs_start_idx=(*end)->dfs_start_idx;
            dfs_end_idx=(*start)->dfs_end_idx;
            for (; start<end; start++) {
                update_level(*start);
            }
    }

    void update_level(const range_tree_temp* other){
        for (int i=0; i<4; i++) {
            min_level[i]=std::min(min_level[i],other->min_level[i]);
        }
    }
};
static bool to_swap(range_tree_temp* possible_descendant,size_t range_start, size_t range_end){
    auto ele_idx=possible_descendant->dfs_start_idx;
    assert(ele_idx!=range_start);
    if (ele_idx>range_start&&ele_idx<=range_end) {
        assert(possible_descendant->dfs_end_idx<=range_end);
        return true;
    }else {
        assert(possible_descendant->dfs_end_idx<range_start||ele_idx>range_end);
        return false;
    }
}

//Make a node correspond to at most a node in the actual tree
void cover_all_remaining(std::vector<range_tree_temp *> &to_merge,
               MAT::Node *&covering_node_parent) {
    auto temp = new range_tree_temp(covering_node_parent, to_merge.begin(),
                                    to_merge.end());
    to_merge.clear();
    to_merge.push_back(temp);
}
static MAT::Node *
merge_children(const std::vector<MAT::Node *> &dfs_ordered_nodes,
               std::vector<range_tree_temp *> &to_merge, size_t min_idx) {
    auto out_iter=to_merge.begin();
    auto scanning_iter=to_merge.begin();
    MAT::Node* covering_node_parent;
    while(scanning_iter<to_merge.end()){
        auto iter=scanning_iter;
        auto covering_node=dfs_ordered_nodes[(*iter)->dfs_start_idx];
        assert(covering_node->parent);
        covering_node_parent=covering_node->parent;
        auto par_range_start=covering_node_parent->dfs_index;
        auto par_range_end=covering_node_parent->dfs_end_index;
        while (par_range_start>min_idx||par_range_end<min_idx) {
            while (iter<to_merge.end()&&(*iter)->dfs_start_idx>par_range_start) {
                assert((*iter)->dfs_end_idx<par_range_end);
                iter++;
            }
            covering_node=covering_node->parent;
            covering_node_parent=covering_node->parent;
            par_range_start=covering_node_parent->dfs_index;
            par_range_end=covering_node_parent->dfs_end_index;
        }
        assert(!(min_idx>=covering_node->dfs_index&&min_idx<=covering_node->dfs_end_index));
        assert(out_iter<scanning_iter);
        if (iter==scanning_iter) {
            *out_iter=*iter;
        }else {
            *out_iter=new range_tree_temp(covering_node,scanning_iter,iter);
        }
        out_iter++;
        scanning_iter=iter;
        assert(out_iter<to_merge.end());
    }
    if (out_iter==to_merge.end()) {
        cover_all_remaining(to_merge, covering_node_parent);
        return covering_node_parent;
    }else {
        to_merge.erase(out_iter,to_merge.end());
        return nullptr;
    }
}
static MAT::Node* condense_children(std::vector<range_tree_temp*>& to_condense,const node_info& range,const std::vector<MAT::Node*>& dfs_ordered_nodes){
    auto range_start=range.dfs_idx;
    auto range_end=dfs_ordered_nodes[range_start]->dfs_end_index;
    range_tree_temp* new_ele=new range_tree_temp(dfs_ordered_nodes[range_start]);
    new_ele->min_level[range.base]=range.level;
    auto ele_end=to_condense.end()-1;
    auto ele_start=to_condense.begin();
    while (ele_start<ele_end) {
        if (to_swap(*ele_start, range_start, range_end)) {
            new_ele->update_level(*ele_start);
            std::swap(*ele_start,*ele_end);
            ele_end--;
        }else {
            ele_start++;
        }
    }
    if (!to_swap(*ele_end, range_start, range_end)) {
        ele_end++;
    }else {
        new_ele->update_level(*ele_end);
    }
    new_ele->children=std::vector<range_tree_temp*>(ele_end,to_condense.end());
    to_condense.erase(ele_end,to_condense.end());
    MAT::Node* out=nullptr;
    if (to_condense.size()==MAX_DIST) {
        out=merge_children(dfs_ordered_nodes, to_condense, range.dfs_idx);
    }
    to_condense.push_back(new_ele);
    return out;
}
static void flatten(std::vector<range_tree_temp*>& in, std::vector<range_tree_node>& out){}
void make_range_tree(const std::vector<MAT::Node*>& dfs_ordered_nodes,const tbb::concurrent_vector<node_info>& in,std::vector<range_tree_node>& out){
    std::vector<range_tree_temp*> curr_top;
    std::stack<MAT::Node*> merging_points({dfs_ordered_nodes[0]});
    auto stack_top_ramge_start=0;
    auto stack_top_ramge_end=merging_points.top()->dfs_end_index;
    curr_top.reserve(MAX_DIST);
    for (auto iter=in.rbegin(); iter<in.rend(); iter++) {
        if (iter->dfs_idx<stack_top_ramge_start||iter->dfs_idx>stack_top_ramge_end) {
            assert(curr_top.size()>1);
            cover_all_remaining(curr_top, merging_points.top());
            merging_points.pop();
            stack_top_ramge_start=merging_points.top()->dfs_index;
            stack_top_ramge_end=merging_points.top()->dfs_end_index;
        }
        if (auto node=condense_children(curr_top, *iter, dfs_ordered_nodes)) {
            assert(node->dfs_index>merging_points.top()->dfs_index);
            merging_points.push(node);
            stack_top_ramge_start=merging_points.top()->dfs_index;
            stack_top_ramge_end=merging_points.top()->dfs_end_index;
        }
    }
    
}


uint32_t range_tree::find_idx(const MAT::Node* node,uint32_t probe_start_idx_in) const{
    auto probe_start_idx=probe_start_idx_in;
    for (; start_idxes[probe_start_idx]<node->dfs_end_index; probe_start_idx++) {}
    if (start_idxes[probe_start_idx]>node->dfs_end_index) {
        probe_start_idx--;
    }
    if (probe_start_idx<probe_start_idx_in||nodes[probe_start_idx].dfs_end_idx<node->dfs_index) {
        return INT_MAX;
    }
    if (nodes[probe_start_idx].dfs_start_idx<node->dfs_index) {
        find_idx(node,nodes[probe_start_idx].children_start_idx);
    }
    return probe_start_idx;
}
