#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "Twice_Bloom_Filter.hpp"
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
void feed_nodes(int radius,std::vector<MAT::Node *> &to_feed,
                         std::deque<MAT::Node*>& out_queue,
                         std::mutex& out_mutex,
                         std::condition_variable& out_pushed,
                         const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;

#define output_node(this_node) \
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
    std::sort(this_round_idx.begin(), this_round_idx.end());
    for (auto iter = this_round_idx.begin(); iter < this_round_idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (!check_in_radius_or_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter],radius)) {
            assert(!check_grand_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter]));
            MAT::Node* this_node=dfs_ordered_nodes[*iter];
            output_node(this_node);
        } else {
            next_round_idx.push_back(*iter);
        }
    }
    MAT::Node* last_node=dfs_ordered_nodes[this_round_idx.back()];
    output_node(last_node);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    }

    {
        std::lock_guard<std::mutex> lk(out_mutex);
        out_queue.push_back(nullptr);
        out_pushed.notify_one();
    }
}