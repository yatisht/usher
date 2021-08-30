#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <iostream>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <stack>
#include <queue>
#include <tbb/task.h>
#include <utility>
#include <vector>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8
using Mutation_Annotated_Tree::Node;
tbb::concurrent_unordered_map<std::string, uint8_t>  Mutation_Annotated_Tree::Mutation::chromosome_map;
std::vector<std::string>  Mutation_Annotated_Tree::Mutation::chromosomes;
std::mutex Mutation_Annotated_Tree::Mutation::ref_lock;
std::vector<nuc_one_hot> Mutation_Annotated_Tree::Mutation::refs;



/* === Tree === */
std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;

    if (nid == "") {
        if (root == NULL) {
            return traversal;
        }
        nid = root->identifier;
    }

    Node* node = all_nodes[nid];
    size_t idx=0;
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        curr_node->bfs_index=idx++;
        traversal.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }

    return traversal;
}

static void depth_first_expansion_helper(Mutation_Annotated_Tree::Node* node, std::vector<Mutation_Annotated_Tree::Node*>& vec, size_t& index,size_t level) {
#ifdef DETAIL_DEBUG_NO_LOOP
    assert(std::find(vec.begin(),vec.end(),node)==vec.end());
#endif
    vec.push_back(node);
    node->level=level;
    //assert(vec.size()-1==index);
    node->dfs_index=index;
    index++;
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec,index,level+1);
    }
    node->dfs_end_index=index-1;
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::depth_first_expansion(Mutation_Annotated_Tree::Node* node) const {
    TIMEIT();
    std::vector<Node*> traversal;
    if (node == NULL) {
        node = root;
    }
    size_t index=0;
    if (node == NULL) {
        return traversal;
    }
    depth_first_expansion_helper(node, traversal,index,0);
    return traversal;
}

size_t Mutation_Annotated_Tree::Tree::get_parsimony_score() {
    size_t score = 0;
    auto dfs = depth_first_expansion();
    for (auto n: dfs) {
        for (const auto& mut:n->mutations) {
            score+=mut.is_valid();
        }
    }
    return score;
}

void Mutation_Annotated_Tree::Tree::uncondense_leaves() {
    for (auto cn = condensed_nodes.begin(); cn != condensed_nodes.end(); cn++) {

        auto n = get_node(cn->first);
        auto par = (n->parent != NULL) ? n->parent : n;

        size_t num_samples = cn->second.size();

        if (num_samples > 0) {
            rename_node(n->identifier, cn->second[0]);
        }

        for (size_t s = 1; s < num_samples; s++) {
            create_node(cn->second[s], par, n->branch_length);
        }
    }
    condensed_nodes.clear();
}

void Node::delete_this() {
    for(Node* n:children) {
        n->delete_this();
    }
    delete this;
}
void Mutation_Annotated_Tree::Tree::delete_nodes() {
    root->delete_this();
}
static void get_leaves_helper(Node* root, std::vector<Node*>& out){
    for(auto child:root->children){
        if (child->is_leaf()) {
            out.push_back(child);
        }else {
            get_leaves_helper(child, out);
        }
    }
}
std::vector<Node*> Mutation_Annotated_Tree::Tree::get_leaves() const{
    std::vector<Node*> out;
    get_leaves_helper(root, out);
    return out;
}

void Mutation_Annotated_Tree::Tree::populate_ignored_range(){
    auto leaves=get_leaves();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,leaves.size()),[&leaves](const tbb::blocked_range<size_t>& range){
        for (auto idx=range.begin(); idx<range.end(); idx++) {
            auto leaf=leaves[idx];
            leaf->populate_ignored_range();
        }
    });
}