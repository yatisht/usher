#include "tree_rearrangement_internal.hpp"
#include "check_samples.hpp"
#include "Fitch_Sankoff.hpp"
#include "src/matOptimize/Twice_Bloom_Filter.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <unordered_set>
#include <random>
#include <utility>
namespace MAT=Mutation_Annotated_Tree;
bool search_all_dir=true;
//add a root above current root, so nodes can move above the current node
void add_root(MAT::Tree *tree) {
    MAT::Node *old_root = tree->root;
    MAT::Node *new_root = new MAT::Node();
    new_root->identifier = std::to_string(++tree->curr_internal_node);
    std::string &node_name = new_root->identifier;
    new_root->children.push_back(old_root);
    tree->all_nodes.emplace(node_name, new_root);
    old_root->parent=new_root;
    tree->root=new_root;
}
//Usher expect parent of condensed node have no mutation, so before outputing usher compatible protobuf,
//add intermediate nodes to carry mutations of condensed nodes
void fix_condensed_nodes(MAT::Tree *tree) {
    std::vector<MAT::Node *> nodes_to_fix;
    for (auto iter : tree->all_nodes) {
        if (tree->condensed_nodes.count(iter.first) &&
                (!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for (auto node : nodes_to_fix) {
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier,
                          std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier, node);
    }
}

static void mark_changed_neighbor_descendent(int radius_left,MAT::Node* root){
    for (auto child : root->children) {
        child->set_ancestor_changed();
        if (radius_left>0&&(!child->get_self_changed())) {
            mark_changed_neighbor_descendent(radius_left-1, child);        
        }
    }
}

//see whether within radius of root have nodes that have same mutation as root and have changed in previous iteration (or this is the first iteration)
static void mark_changed_neighbor_self(int radius, MAT::Node* root) {
    root->set_self_changed();
    for (int radius_left=radius;radius_left>0;radius_left--){
        root=root->parent;
        if (!root||root->get_self_changed()) {
            break;
        }
        root->set_descendent_changed();
        mark_changed_neighbor_descendent(radius_left, root);
    }
}

//find src node to search
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                        std::vector<MAT::Node *> &output,
                        bool search_all_node, int radius_in,MAT::Tree &tree) {
    auto start=std::chrono::steady_clock::now();
    unsigned int radius=abs(radius_in);
    output.clear();
    if (!search_all_dir) {
        for(auto node:bfs_ordered_nodes) {
            node->clear_changed();
        }
        //mark surrounding of changed nodes
        tbb::parallel_sort(changed_nodes.begin(),changed_nodes.end());
        auto end_iter=std::unique(changed_nodes.begin(),changed_nodes.end());
        changed_nodes.erase(end_iter,changed_nodes.end());
        fprintf(stderr, "%zu changed nodes \n",changed_nodes.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0,changed_nodes.size()),[radius,&tree](const tbb::blocked_range<size_t>& r) {
            for (size_t idx=r.begin(); idx<r.end(); idx++) {
                auto node=tree.get_node(changed_nodes[idx]);
                if (!node) {
                    continue;
                }
                mark_changed_neighbor_self(radius, node);
            }
        });
    }
    if (search_all_node) {
        output=bfs_ordered_nodes;
        fprintf(stderr, "Search all nodes\n");
    } else {
        fprintf(stderr, "find max_level\n");
        auto max_level=tree.get_max_level();
        fprintf(stderr, "Max level %zu\n",max_level);
        if (radius>2*max_level) {
            output=bfs_ordered_nodes;
            fprintf(stderr, "Search all nodes\n");
        }
        output.reserve(bfs_ordered_nodes.size());
        for(auto node:bfs_ordered_nodes) {
            if (node->have_change_in_neighbor()) {
                output.push_back(node);
            }
        }
    }
    changed_nodes.clear();
    fprintf(stderr, "Will search %f of nodes\n",(double)output.size()/(double)bfs_ordered_nodes.size());
    std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now()-start;
    fprintf(stderr, "Took %f s to find nodes to move\n",elapsed_seconds.count());
}
void save_final_tree(MAT::Tree &t,
                     const std::string &output_path) {
    std::vector<MAT::Node *> dfs = t.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
    [&dfs](tbb::blocked_range<size_t> r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            dfs[i]->mutations.remove_invalid();
        }
    });
    fix_condensed_nodes(&t);
    fprintf(stderr, "%zu condensed_nodes\n",t.condensed_nodes.size());
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, output_path);
}
