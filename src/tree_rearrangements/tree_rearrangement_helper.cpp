#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "Twice_Bloom_Filter.hpp"
#include <cstddef>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output){
    Twice_Bloom_Filter filter;
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            filter.insert(m.position);
        }
    }
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            if(filter.query(m.position)){
                output.push_back(n);
                break;
            }
        }
    }
}
#ifndef NDEBUG
void BFS(MAT::Node* src,MAT::Node* ori_src, int radius,std::vector<MAT::Node*>& out);
#endif
/*
static bool check_in_radius_or_parent(MAT::Node* later,MAT::Node* earlier, int radius){
    bool parent_possible=true;
    assert(later->index>earlier->index);
    while (radius>0||parent_possible) {
        if (later->index>earlier->index) {
            later=later->parent;
            if (!later) {
                return false;
            }
        }else if(later->index<earlier->index){
            earlier=earlier->parent;
            if (!earlier) {
                return false;
            }
            parent_possible=false;
        }else {
            assert(later->index==earlier->index);
            return true;
        }
        radius--;
    }
    return false;
}
*/
typedef std::vector<unsigned short> Reachable_Map_t;
/*
#define OUTSIDE ((char)0)
#define IN2RADIUS ((char)1)
#define INRADIUS ((char)2)
#define NODE_SELF ((char)3)
static bool check_in_radius(MAT::Node* this_node, int radius,Reachable_Map_t& reachable_map){
    size_t this_node_idx=this_node->index;
    if (reachable_map[this_node_idx]) {
        return true;
    }
    MAT::Node* ancestor=this_node->parent;
    std::vector<size_t> in_radius_idx;
    for (char n_hop=0; n_hop<radius; n_hop++) {
        if (!ancestor) {
            break;
        }
        in_radius_idx.push_back(ancestor->index);
        if (reachable_map[ancestor->index]) {
            return true;
        }
        ancestor=ancestor->parent;
    }
    std::vector<size_t> in_2radius_idx;
    for (char n_hop=0; n_hop<radius; n_hop++) {
        if (!ancestor) {
            break;
        }
        if (reachable_map[ancestor->index]>IN2RADIUS) {
            return true;
        }
        in_2radius_idx.push_back(ancestor->index);
        ancestor=ancestor->parent;
    }

    for(auto this_idx:in_radius_idx){
        reachable_map[this_idx]=std::max(reachable_map[this_idx],INRADIUS);
    }

    for(auto this_idx:in_2radius_idx){
        reachable_map[this_idx]=std::max(reachable_map[this_idx],IN2RADIUS);
    }

    reachable_map[this_node_idx]=NODE_SELF;
    return false;
}
*/
static bool check_in_radius(MAT::Node* this_node, int radius,Reachable_Map_t& reachable_map){
    std::vector<size_t> parent_idx;
    parent_idx.reserve(radius);
    for (size_t dist=0; dist<radius; dist++) {
        if (!this_node) {
            break;
        }
        size_t this_node_idx=this_node->index;
        if(dist+reachable_map[this_node_idx]<radius){
            return true;
        }
        parent_idx.push_back(this_node_idx);
        this_node=this_node->parent;
    }
    for (unsigned short idx=0; idx<parent_idx.size(); idx++) {
        reachable_map[parent_idx[idx]]=std::min(reachable_map[parent_idx[idx]],idx);
    }
    return false;
}
void feed_nodes(int radius,std::vector<MAT::Node *> &to_feed,
                         std::deque<MAT::Node*>& out_queue,
                         std::mutex& out_mutex,
                         std::condition_variable& out_pushed,
                         const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;
    #ifdef DETAIL_DEBUG_NODE_INPUT
    std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> pushed_nodes;
    #endif

#define __output_node(this_node) \
{std::lock_guard<std::mutex> lk(out_mutex);\
out_queue.push_back(this_node);\
out_pushed.notify_one();}\
this_node=this_node->parent;\
if(this_node->parent){\
    next_round_idx.push_back(this_node->index);\
}

#ifdef DETAIL_DEBUG_NODE_INPUT
#define output_node(this_node)\
assert(pushed_nodes.insert(this_node).second);\
__output_node(this_node)
#else
#define output_node(this_node) __output_node(this_node)
#endif

    {
        this_round_idx.reserve(to_feed.size());
        for (auto node : to_feed) {
            auto parent = node->parent;
            if (parent)
                this_round_idx.push_back(node->index);
        }
    }
    while (!this_round_idx.empty()) {
        int pushed_node=0;
        //Reachable_Map_t reachable_map(dfs_ordered_nodes.size(),OUTSIDE);
        Reachable_Map_t reachable_map(dfs_ordered_nodes.size(),2*radius+2);
    #ifdef DETAIL_DEBUG_NODE_INPUT
        std::unordered_map<MAT::Node*,MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> pushed_nodes_none_neighbor;
        #endif
    std::sort(this_round_idx.begin(), this_round_idx.end());
    for (auto iter = this_round_idx.begin(); iter < this_round_idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (check_grand_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter])) {
            next_round_idx.push_back(*iter);
        } else {
            MAT::Node *this_node = dfs_ordered_nodes[*iter];
    #ifdef DETAIL_DEBUG_NODE_INPUT
                std::vector<MAT::Node *> other_dsts;
                BFS(this_node, this_node, radius, other_dsts);
                std::vector<char> state_before(other_dsts.size());
                for (size_t node_idx=0; node_idx<other_dsts.size(); node_idx++) {
                    state_before[node_idx]=reachable_map[other_dsts[node_idx]->index];
                }
            #endif
            if (check_in_radius(dfs_ordered_nodes[*iter], 2*radius+1,reachable_map)) {
                next_round_idx.push_back(*iter);
            } else {
                output_node(this_node);
#ifdef DETAIL_DEBUG_NODE_INPUT
                for(size_t reachable_node_idx=0;reachable_node_idx<other_dsts.size();reachable_node_idx++){
                    const auto reachable_node=other_dsts[reachable_node_idx];
                    auto insert_result=pushed_nodes_none_neighbor.emplace(reachable_node,this_node);
                    if(!insert_result.second){
                        fprintf(stderr, "At Node index %zu: Neighbor node index %zu already reachable from %zu, state before insert %d, state after insert %d.\n",this_node->index,reachable_node->index,insert_result.first->second->index,(int)state_before[reachable_node_idx],(int)reachable_map[reachable_node->index]);
                        assert(false);
                    }
                }
#endif
                pushed_node++;
            }
        }
    }
    MAT::Node* last_node=dfs_ordered_nodes[this_round_idx.back()];
    output_node(last_node);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    #ifdef DETAIL_DEBUG_NODE_INPUT
    fprintf(stderr,"pushed %d nodes, %zu nodes left\n",pushed_node,this_round_idx.size());
    #endif
    }

    {
        std::lock_guard<std::mutex> lk(out_mutex);
        out_queue.push_back(nullptr);
        out_pushed.notify_one();
    }
}