#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "Twice_Bloom_Filter.hpp"
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
#define OUTSIDE 0
#define INRADIUS 1
#define IN2RADIUS 2
#define NODE_SELF 3
typedef std::unordered_map<size_t, char> Reachable_Map_t;
static bool check_in_radius(MAT::Node* this_node, int radius,Reachable_Map_t& reachable_map){
    Reachable_Map_t::iterator this_node_iter=reachable_map.emplace(this_node->index,OUTSIDE).first;
    if (this_node_iter->second) {
        return true;
    }
    MAT::Node* ancestor=this_node->parent;
    std::vector<Reachable_Map_t::iterator> in_radius_iters;
    for (char n_hop=0; n_hop<radius; n_hop++) {
        if (!ancestor) {
            break;
        }
        Reachable_Map_t::iterator ancestor_iter=reachable_map.emplace(ancestor->index,OUTSIDE).first;
        in_radius_iters.push_back(ancestor_iter);
        if (ancestor_iter->second) {
            return true;
        }
        ancestor=ancestor->parent;
    }
    std::vector<Reachable_Map_t::iterator> in_2radius_iters;
    for (char n_hop=0; n_hop<radius; n_hop++) {
        if (!ancestor) {
            break;
        }
        Reachable_Map_t::iterator ancestor_iter=reachable_map.emplace(ancestor->index,IN2RADIUS).first;
        if (ancestor_iter->second==NODE_SELF) {
            return true;
        }else {
            ancestor_iter->second=IN2RADIUS;
        }
        ancestor=ancestor->parent;
    }
    for(auto this_iter:in_radius_iters){
        this_iter->second=INRADIUS;
    }
    for(auto this_iter:in_2radius_iters){
        this_iter->second=IN2RADIUS;
    }
    this_node_iter->second=NODE_SELF;
    return false;
}
void feed_nodes(int radius,std::vector<MAT::Node *> &to_feed,
                         std::deque<MAT::Node*>& out_queue,
                         std::mutex& out_mutex,
                         std::condition_variable& out_pushed,
                         const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;
    #ifndef NDEBUG
    std::unordered_set<MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> pushed_nodes;
    #endif

#define output_node(this_node) \
assert(pushed_nodes.insert(this_node).second);\
{std::lock_guard<std::mutex> lk(out_mutex);\
out_queue.push_back(this_node);\
out_pushed.notify_one();}\
this_node=this_node->parent;\
if(this_node->parent){\
    next_round_idx.push_back(this_node->index);\
}

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
        Reachable_Map_t reachable_map;
    std::sort(this_round_idx.begin(), this_round_idx.end());
    for (auto iter = this_round_idx.begin(); iter < this_round_idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (check_grand_parent(dfs_ordered_nodes[next_ele],
                                       dfs_ordered_nodes[*iter])||check_in_radius(
                                       dfs_ordered_nodes[*iter], radius,
                                       reachable_map)) {
            next_round_idx.push_back(*iter);
        } else {
            MAT::Node *this_node = dfs_ordered_nodes[*iter];
            output_node(this_node);
            pushed_node++;
        }
    }
    MAT::Node* last_node=dfs_ordered_nodes[this_round_idx.back()];
    output_node(last_node);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    //fprintf(stderr,"pushed %d nodes\n",pushed_node);
    }

    {
        std::lock_guard<std::mutex> lk(out_mutex);
        out_queue.push_back(nullptr);
        out_pushed.notify_one();
    }
}