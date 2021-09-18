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
#include <memory>
#define MAX_DIST 16
std::vector<range_tree> addable_idxes;
struct range_tree_temp{
    uint32_t dfs_start_idx;
    uint32_t dfs_end_idx;
    std::array<uint16_t, 4> min_level;
    std::vector<std::shared_ptr<range_tree_temp>> children;
    uint16_t level;
    uint16_t secondary_level;
    bool is_ori_node;
    const MAT::Node* covering_node;
    range_tree_temp(){
        assert(false);
    }
    void init(const MAT::Node* node){
        level=node->level;
        secondary_level=0;
        for (int i=0; i<4; i++) {
            min_level[i]=UINT16_MAX;
        }
        covering_node=node;
    }
    range_tree_temp(const MAT::Node* node){
        assert(node->parent);
        dfs_start_idx=node->dfs_index;
        dfs_end_idx=node->dfs_end_index;
        init(node);
        is_ori_node=false;
    }
    range_tree_temp(const node_info& in,const std::vector<MAT::Node*>dfs_ordered_nodes):range_tree_temp(dfs_ordered_nodes[in.dfs_idx]){
        min_level[in.base]=in.level;
        is_ori_node=true;
    }
    void process_child(const MAT::Node* node){
        assert(children.size()<=MAX_DIST);
        #ifndef NDEBUG
        int last_end=node->dfs_index;
        for(const auto& child:children){
            assert(child->dfs_start_idx>=last_end);
            last_end=child->dfs_end_idx;
        }
        assert(last_end<=node->dfs_end_index);
        #endif
        dfs_start_idx=children.front()->dfs_start_idx;
        dfs_end_idx=children.back()->dfs_end_idx;
    }
    template<typename iter_t>
    range_tree_temp(iter_t iter,iter_t end,const MAT::Node* node){
        init(node);
        children.reserve(end-iter);
        for(;iter<end;iter++){
            update_level(*iter);
            secondary_level=std::max(secondary_level,(*iter)->secondary_level);
            children.push_back(*iter);
        }
        secondary_level++;
        process_child(node);
        is_ori_node=false;
    }
    range_tree_temp(std::vector<std::shared_ptr<range_tree_temp>>& new_children,int node_size,const MAT::Node* node){
            struct idx_cmp{
                bool operator()(const std::shared_ptr<range_tree_temp>& a, const std::shared_ptr<range_tree_temp>& b) const{
                    return a->dfs_start_idx<b->dfs_start_idx;
                }
            };
        init(node);
        if (node_size<=MAX_DIST) {
        std::sort(new_children.begin(),new_children.end(),idx_cmp());
        children.reserve(node_size);
        for(auto iter=new_children.begin();iter<new_children.end();iter++){
            update_level(*iter);
            if ((*iter)->is_ori_node) {
                children.push_back(*iter);
            }else{
            children.insert(children.end(),(*iter)->children.begin(),(*iter)->children.end());
            }
        }
        assert(children.size()<=node_size);
        }else if (new_children.size()<=MAX_DIST) {
            children.reserve(MAX_DIST);
            struct size_cmp{
                bool operator()(const std::shared_ptr<range_tree_temp>a, const std::shared_ptr<range_tree_temp> b) const{
                    return a->children.size()>b->children.size();
                }
            };

            std::sort(new_children.begin(),new_children.end(),size_cmp());
            for (const auto& child : new_children) {
                update_level(child);
                if (node_size>MAX_DIST) {
                    children.push_back(child);
                    node_size-=(child->children.size()-1);
                }else {
                    if (child->is_ori_node) {
                        children.push_back(child);
                    }else {
                        children.insert(children.end(),child->children.begin(),child->children.end());                    
                    }
                }
            }
            std::sort(children.begin(),children.end(),idx_cmp());
        assert(children.size()<=node_size);
        }else {
            std::sort(new_children.begin(),new_children.end(),idx_cmp());
            while (true) {
                auto iter=new_children.begin();
                auto end=new_children.end();
                while (end-iter>MAX_DIST) {
                    children.emplace_back(new range_tree_temp(iter,iter+MAX_DIST,node));
                    iter+=MAX_DIST;
                }
                if (iter!=end) {
                    children.emplace_back(new range_tree_temp(iter,end,node));
                }
                if (children.size()>MAX_DIST) {
                    children.swap(new_children);
                    children.clear();
                }else {
                    break;
                }
            }
            for (auto child : children) {
                update_level(child);
            }
        }
        process_child(node);
        is_ori_node=false;
    }

    void update_level(const std::shared_ptr<range_tree_temp>& other){
        for (int i=0; i<4; i++) {
            min_level[i]=std::min(min_level[i],other->min_level[i]);
        }
    }
    void cover(const std::shared_ptr<range_tree_temp>& other){
        update_level(other);
        if (children.empty()) {
            if (other->is_ori_node) {
                children.push_back(other);
            }else {
                assert(!other->children.empty());
                children.swap(other->children);
            }
            return;
        }else if (other->is_ori_node) {
            children.push_back(other);
            return;
        }
        //Merge
        std::vector<std::shared_ptr<range_tree_temp>> new_children;
        new_children.reserve(children.size()+other->children.size());
        auto iter=children.begin();
        auto other_start_idx=other->dfs_start_idx;
        while (iter<children.end()) {
            if ((*iter)->dfs_end_idx<other_start_idx) {
                new_children.push_back(*iter);
            }else {
                assert(other_start_idx!=(*iter)->dfs_end_idx);
                break;
            }
            iter++;
        }
        new_children.insert(new_children.end(),other->children.begin(),other->children.end());
            assert(iter==children.end()||other->dfs_end_idx<(*iter)->dfs_start_idx);        
        new_children.insert(new_children.end(),iter,children.end());
        children.swap(new_children);
        
    }
    bool can_cover(const std::shared_ptr<range_tree_temp>& other) const{
        return dfs_start_idx<other->dfs_start_idx&&dfs_end_idx>=other->dfs_end_idx;
    }

    range_tree_node out(uint16_t par_idx){
        return range_tree_node{dfs_start_idx,dfs_end_idx,min_level,UINT16_MAX,par_idx,level};
    }
};

struct bfs_queue_content{
    std::shared_ptr<range_tree_temp> child_to_push;
    size_t par_idx;
};
struct priority_bfs_queue_comp{
    bool operator()(const bfs_queue_content& a, const bfs_queue_content& b) const{
        if (a.child_to_push->level>b.child_to_push->level) {
            return true;
        }else if (a.child_to_push->level==b.child_to_push->level&&a.child_to_push->secondary_level<b.child_to_push->secondary_level) {
            return true;
        }else if (a.child_to_push->level==b.child_to_push->level&&a.child_to_push->secondary_level==b.child_to_push->secondary_level&&a.child_to_push->dfs_start_idx>b.child_to_push->dfs_start_idx) {
            return true;
        }
        return false;
    }
};
/*static void set_secondary_level(std::shared_ptr<range_tree_temp>& node,uint8_t par_level,uint8_t par_secondary_level){
    if (node->level==par_level) {
        node->secondary_level=par_secondary_level+1;
    }else{
        node->secondary_level=0;
    }
    for (auto& child : node->children) {
        set_secondary_level(child, node->level, node->secondary_level);
    }
}*/
typedef std::priority_queue<bfs_queue_content,std::vector<bfs_queue_content>,priority_bfs_queue_comp> priority_bfs_queue;
static void flatten(std::shared_ptr<range_tree_temp> top, std::vector<range_tree_node>& out,size_t node_count){
    out.reserve(node_count*2);
    //set_secondary_level(top, -1, -1);
    std::vector<bfs_queue_content> container;
    container.reserve(node_count);
    out.push_back(top->out(UINT16_MAX));
    if (top->children.size()==0) {
        return;
    }
    container.push_back(bfs_queue_content{top,0});
    priority_bfs_queue queue(priority_bfs_queue_comp(),std::move(container));
#ifndef NDEBUG
    int last_level=-1;
#endif
    while (!queue.empty()) {
        auto top=queue.top();
        queue.pop();
        if (top.child_to_push->children[0]->dfs_start_idx<out.back().dfs_end_idx) {
            //assert(top.child_to_push->children[0]->level>last_level);
            out.push_back(range_tree_node{UINT32_MAX,UINT32_MAX});
        }
        out[top.par_idx].children_start_idx=out.size();
        for(auto& top_content:top.child_to_push->children){
            assert(top_content->dfs_start_idx>=out.back().dfs_end_idx||out.back().dfs_end_idx==UINT32_MAX);
            if (!top_content->children.empty()) {
                queue.push(bfs_queue_content{top_content,out.size()});                
            }
            out.push_back((top_content)->out(top.par_idx));
        }
        #ifndef NDEBUG
            last_level=top.child_to_push->level;
        #endif
    }
}
struct temp_tree_build_comp{
    bool operator()(const std::shared_ptr<range_tree_temp>& a,const std::shared_ptr<range_tree_temp>& b)const{
        auto a_idx=a->covering_node->dfs_index;
        auto b_idx=b->covering_node->dfs_index;
        if (a_idx<b_idx) {
            return true;
        }else if (a_idx==b_idx&&a->dfs_start_idx>b->dfs_start_idx) {
            return true;;
        }
        return false;
    }
};
typedef std::priority_queue<std::shared_ptr<range_tree_temp>,std::vector<std::shared_ptr<range_tree_temp>>,temp_tree_build_comp> tree_build_heap_t;
void make_range_tree(const std::vector<MAT::Node*>& dfs_ordered_nodes,tbb::concurrent_vector<node_info>& in,range_tree& out,size_t idx){
    if (in.empty()) {
        return;
    }
    size_t node_count=0;
    std::vector<std::shared_ptr<range_tree_temp>> content;
    content.reserve(in.size());
    for (const auto& node : in) {
        content.emplace_back(new range_tree_temp(node,dfs_ordered_nodes));
    }

    tree_build_heap_t heap(temp_tree_build_comp(),std::move(content));
    std::shared_ptr<range_tree_temp> top;
    while (true) {
        top=heap.top();
        std::vector<std::shared_ptr<range_tree_temp>> child_nodes({top});
        heap.pop();
        int node_count=std::max(top->children.size(),(size_t)1);
        if (top->covering_node->dfs_index==73050&&idx==241) {
            fputc('a', stderr);
        }
        while (!heap.empty()&&heap.top()->covering_node==top->covering_node) {
            if (heap.top()->dfs_start_idx==top->dfs_start_idx) {
                assert(heap.top()->dfs_end_idx==top->dfs_end_idx);
                top->update_level(heap.top());
            }else if (top->can_cover(heap.top())) {
                node_count+=std::max(heap.top()->children.size(),(size_t)1);
                if (top->children.empty()) {
                    node_count--;                    
                }
                top->cover(heap.top());
            }else if(heap.top()->can_cover(top)){
                assert(false);
            }else {
                child_nodes.push_back(heap.top());
                node_count+=std::max(heap.top()->children.size(),(size_t)1);
            }
            heap.pop();
        }

        if (child_nodes.size()!=1)  {
            top=std::shared_ptr<range_tree_temp>(new range_tree_temp(child_nodes,node_count,top->covering_node));
        }

        if(!top->covering_node->parent){
            assert(heap.empty());
            break;
        }
        top->covering_node=top->covering_node->parent;
        heap.push(top);
    }
    assert(node_count<INT16_MAX);
    if (idx==66) {
        fputc('a', stderr);
    }
    flatten(top, out.nodes,node_count);
    out.nodes.push_back(range_tree_node{UINT32_MAX,UINT32_MAX});
    out.end_idxes.reserve(out.nodes.size());
    for(const auto& node: out.nodes){
        assert(out.end_idxes.empty()||out.end_idxes.back()==UINT32_MAX||out.end_idxes.back()<=node.dfs_end_idx);
        out.end_idxes.push_back(node.dfs_end_idx);
    }
}

uint16_t range_tree::find_idx(const MAT::Node* node,uint16_t& last_idx,uint16_t probe_start_idx_in) const{
    auto probe_start_idx=probe_start_idx_in;
    if (end_idxes.empty()) {
        return UINT16_MAX;
    }
    for (; end_idxes[probe_start_idx]<node->dfs_index; probe_start_idx++) {}
    if (nodes[probe_start_idx].dfs_start_idx>node->dfs_end_index) {
        return UINT16_MAX;
    }
    if (nodes[probe_start_idx].level<node->level) {
        if (nodes[probe_start_idx].children_start_idx==UINT16_MAX) {
            return UINT16_MAX;
        }
        last_idx=probe_start_idx;
        return find_idx(node,last_idx,nodes[probe_start_idx].children_start_idx);
    }
    return probe_start_idx;
}
