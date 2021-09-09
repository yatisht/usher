#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include "Profitable_Moves_Enumerators.hpp"
#include <iterator>
#include <cstdio>
#include <linux/limits.h>
#include <queue>
#include <stack>
#include <tbb/concurrent_vector.h>
#include <utility>
#include <iterator>
#include <vector>
#define MAX_DIST 16
std::vector<range_tree> addable_idxes;
struct make_tree_temp{
    uint32_t out_start_idx;
    uint32_t out_end_idx;
    uint32_t root_dfs_idx;
};
struct range_tree_temp{
    uint32_t dfs_start_idx;
    uint32_t dfs_end_idx;
    std::array<uint16_t, 4> min_level;
    std::vector<range_tree_temp*> children;
    uint16_t level;
    MAT::Node* covering_node;
    range_tree_temp(){
        assert(false);
    }
    void init(MAT::Node* node){
        level=node->level;
        for (int i=0; i<4; i++) {
            min_level[i]=UINT16_MAX;
        }
        covering_node=node;
    }

    range_tree_temp(MAT::Node* node){
        dfs_start_idx=node->dfs_index;
        dfs_end_idx=node->dfs_end_index;
        init(node);
                if (dfs_start_idx==16030) {
            fputc('a',stderr);
        }
    }

    range_tree_temp(MAT::Node* node,
        std::vector<range_tree_temp*>::const_iterator start,
        std::vector<range_tree_temp*>::const_iterator end): children(std::make_reverse_iterator(end),std::make_reverse_iterator(start)){
            init(node);
            assert(end-start>1);
            #ifndef NDEBUG
            size_t last_end_idx=node->dfs_index;
            #endif
            for (const auto& child : children) {
                assert(child->dfs_start_idx>=last_end_idx);
                update_level(child);
                #ifndef NDEBUG
                    last_end_idx=child->dfs_end_idx;
                #endif
            }
            assert(last_end_idx<=node->dfs_end_index);
            dfs_start_idx=children.front()->dfs_start_idx;
            dfs_end_idx=children.back()->dfs_end_idx;
        if (dfs_start_idx==16030) {
            fputc('a',stderr);
        }
    }

    void update_level(const range_tree_temp* other){
        for (int i=0; i<4; i++) {
            min_level[i]=std::min(min_level[i],other->min_level[i]);
        }
    }

    range_tree_node out(uint16_t par_idx){
        return range_tree_node{dfs_start_idx,dfs_end_idx,min_level,par_idx,UINT16_MAX,level};
    }
};

//Make a node correspond to at most a node in the actual tree
static void cover_all_remaining(std::vector<range_tree_temp *> &to_merge,
               MAT::Node *&covering_node_parent,size_t& node_count) {
    auto temp = new range_tree_temp(covering_node_parent, to_merge.begin(),
                                    to_merge.end());
    to_merge.clear();
    to_merge.push_back(temp);
    node_count++;
}
static MAT::Node *
merge_children(const std::vector<MAT::Node *> &dfs_ordered_nodes,
               std::vector<range_tree_temp *> &to_merge, size_t min_idx,size_t& node_count) {
    auto out_iter=to_merge.begin();
    auto scanning_iter=to_merge.begin();
    MAT::Node* covering_node_parent;
    while(scanning_iter<to_merge.end()){
        auto iter=scanning_iter;
        if ((*iter)->dfs_start_idx==24368) {
            fputc('a',stderr);
        }
        auto covering_node=(*iter)->covering_node;
        assert(covering_node->parent);
        covering_node_parent=covering_node->parent;
        auto par_range_start=covering_node_parent->dfs_index;
        auto par_range_end=covering_node_parent->dfs_end_index;
        while (par_range_start>min_idx||par_range_end<min_idx) {
            while (iter<to_merge.end()&&(*iter)->dfs_start_idx>=par_range_start) {
                assert((*iter)->dfs_end_idx<=par_range_end);
                iter++;
            }
            covering_node=covering_node_parent;
            covering_node_parent=covering_node->parent;
            par_range_start=covering_node_parent->dfs_index;
            par_range_end=covering_node_parent->dfs_end_index;
        }
        assert(!(min_idx>=covering_node->dfs_index&&min_idx<=covering_node->dfs_end_index));
        assert(out_iter<=scanning_iter);
        assert(out_iter<to_merge.end());
        if (iter==scanning_iter) {
            iter++;
        }
        if (iter==scanning_iter+1) {
            *out_iter=*scanning_iter;
        }else {
            *out_iter=new range_tree_temp(covering_node,scanning_iter,iter);
            node_count++;
        }
        out_iter++;
        scanning_iter=iter;
    }
    if (out_iter==to_merge.end()) {
        cover_all_remaining(to_merge, covering_node_parent,node_count);
        return covering_node_parent;
    }else {
        to_merge.erase(out_iter,to_merge.end());
        return nullptr;
    }
}
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
static MAT::Node* condense_children(std::vector<range_tree_temp*>& to_condense,tbb::concurrent_vector<node_info>::const_reverse_iterator& range_iter,const std::vector<MAT::Node*>& dfs_ordered_nodes,size_t& node_count){
    const auto range=*range_iter;
    auto range_start=range.dfs_idx;
    auto range_end=dfs_ordered_nodes[range_start]->dfs_end_index;
    range_tree_temp* new_ele=new range_tree_temp(dfs_ordered_nodes[range_start]);
    while (range_iter->dfs_idx==range.dfs_idx) {
        new_ele->min_level[range_iter->base]=range.level;
        range_iter++;
    }
    if (range.dfs_idx==16030) {
        fputc('a', stderr);
    }
    MAT::Node* out=nullptr;
    if (!to_condense.empty()) {
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
    #ifndef NDEBUG
    size_t last_end_idx=new_ele->dfs_start_idx;
    for ( const auto child:new_ele->children ) {
        assert(child->dfs_start_idx>=last_end_idx);
        last_end_idx=child->dfs_end_idx;  
    }
    assert(new_ele->dfs_end_idx>=last_end_idx);
    #endif
    to_condense.erase(ele_end,to_condense.end());
    if (to_condense.size()==MAX_DIST) {
        out=merge_children(dfs_ordered_nodes, to_condense, range.dfs_idx,node_count);
    }
    }
    to_condense.push_back(new_ele);
    return out;
}
struct bfs_queue_content{
    range_tree_temp* child_to_push;
    size_t par_idx;
};
struct priority_bfs_queue_comp{
    bool operator()(const bfs_queue_content& a, const bfs_queue_content& b) const{
        if (a.child_to_push->level<b.child_to_push->level) {
            return true;
        }else if (a.child_to_push->level==b.child_to_push->level&&a.child_to_push->dfs_start_idx<b.child_to_push->dfs_start_idx) {
            return true;
        }
        return false;
    }
};

typedef std::priority_queue<bfs_queue_content,std::vector<bfs_queue_content>,priority_bfs_queue_comp> priority_bfs_queue;
static void flatten(std::vector<range_tree_temp*>& in, std::vector<range_tree_node>& out,size_t node_count){
    out.reserve(node_count*2);
    std::vector<bfs_queue_content> container;
    container.reserve(node_count);
    for(auto iter=in.rbegin();iter<in.rend();iter++){
        container.push_back(bfs_queue_content{*iter,out.size()});
        out.push_back((*iter)->out(UINT16_MAX));
        assert(out.size()==1||out[out.size()-1].dfs_start_idx>out[out.size()-2].dfs_end_idx);
    }
    priority_bfs_queue queue(priority_bfs_queue_comp(),std::move(container));
    while (!queue.empty()) {
        auto& top=queue.top();
        out[top.par_idx].children_start_idx=in.size();
        for(auto& top_content:top.child_to_push->children){
            queue.push(bfs_queue_content{top_content,in.size()});
            if (top_content->dfs_start_idx<out.back().dfs_end_idx) {
                assert(top_content->level>out.back().level);
                out.push_back(range_tree_node{UINT16_MAX});
            }
            out.push_back((top_content)->out(top.par_idx));
        }
        delete top.child_to_push;
        queue.pop();
    }

}
void make_range_tree(const std::vector<MAT::Node*>& dfs_ordered_nodes,tbb::concurrent_vector<node_info>& in,range_tree& out){
    if (in.empty()) {
        return;
    }
    in.push_back(node_info{0});
    std::sort(in.begin(),in.end());
    size_t node_count=0;
    struct stact_content_t{
        MAT::Node* node;
        std::vector<range_tree_temp*> overflow;
        stact_content_t(MAT::Node* node):node(node){}
        stact_content_t(MAT::Node* node,std::vector<range_tree_temp*>&in):node(node){
            overflow.reserve(MAX_DIST);
            push_no_check(in);
        }
        void output(std::vector<range_tree_temp*>&out,size_t& node_count){
            if(out.size()==16){
                cover_all_remaining(overflow, node,node_count);
            };
            if (overflow.size()+out.size()<=16) {
                out.insert(out.end(),overflow.begin(),overflow.end());
            }else {
                cover_all_remaining(overflow,node,node_count);
                out.push_back(overflow[0]);
            }
        }
        void push(std::vector<range_tree_temp*>& in,size_t& node_count){
            if (overflow.size()==16) {
                cover_all_remaining(overflow,node,node_count);
            }
            push_no_check(in);
        }
        private:
        void push_no_check(std::vector<range_tree_temp*>&in){
            assert(in.size()==1);
            overflow.push_back(in[0]);
            in.clear();
        }
    };
    std::vector<range_tree_temp*> curr_top;
    std::stack<stact_content_t> merging_points({dfs_ordered_nodes[0]});
    auto stack_top_ramge_start=0;
    auto stack_top_ramge_end=merging_points.top().node->dfs_end_index;
    curr_top.reserve(MAX_DIST);
    tbb::concurrent_vector<node_info>::const_reverse_iterator iter = in.rbegin();
    while( iter < in.rend()-1) {
        if (iter->dfs_idx < stack_top_ramge_start ||
            iter->dfs_idx > stack_top_ramge_end) {
            assert(curr_top.size() > 1);
            merging_points.top().output(curr_top,node_count);
            cover_all_remaining(curr_top, merging_points.top().node,node_count);
            merging_points.pop();
            stack_top_ramge_start = merging_points.top().node->dfs_index;
            stack_top_ramge_end = merging_points.top().node->dfs_end_index;
        }
        if (auto node = condense_children(curr_top, iter, dfs_ordered_nodes,node_count)) {
            if (node == merging_points.top().node) {
                merging_points.top().push(curr_top,node_count);
            } else {
                assert(node->dfs_index > merging_points.top().node->dfs_index);
                merging_points.emplace(node, curr_top);
                stack_top_ramge_start = merging_points.top().node->dfs_index;
                stack_top_ramge_end = merging_points.top().node->dfs_end_index;
            }
        }
    }
    assert(node_count<INT16_MAX);
    flatten(curr_top, out.nodes,node_count);
    out.start_idxes.reserve(out.nodes.size());
    for(const auto& node: out.nodes){
        assert(out.start_idxes.empty()||out.start_idxes.back()==UINT16_MAX||out.start_idxes.back()<node.dfs_start_idx);
        out.start_idxes.push_back(node.dfs_start_idx);
    }
}

uint16_t range_tree::find_idx(const MAT::Node* node,uint32_t probe_start_idx_in) const{
    auto probe_start_idx=probe_start_idx_in;
    for (; start_idxes[probe_start_idx]<node->dfs_end_index; probe_start_idx++) {}
    if (start_idxes[probe_start_idx]>node->dfs_end_index) {
        probe_start_idx--;
    }
    if (probe_start_idx<probe_start_idx_in||nodes[probe_start_idx].dfs_end_idx<node->dfs_index) {
        return -1;
    }
    if (nodes[probe_start_idx].dfs_start_idx<node->dfs_index) {
        find_idx(node,nodes[probe_start_idx].children_start_idx);
    }
    return probe_start_idx;
}
