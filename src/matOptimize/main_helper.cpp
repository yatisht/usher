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
#include <vector>
namespace MAT=Mutation_Annotated_Tree;
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

static void mark_changed_neighbor_descendent(int radius_left,MAT::Node* root,bool all_node_reachable) {
    for (auto child : root->children) {
        child->set_ancestor_changed();
        if (radius_left>0
                &&(!child->get_self_changed())
                /*&&(!(all_node_reachable&&child->get_ancestor_changed()))*/
           ) {
            mark_changed_neighbor_descendent(radius_left-1, child,all_node_reachable);
        }
    }
}

//see whether within radius of root have nodes that have same mutation as root and have changed in previous iteration (or this is the first iteration)
static void mark_changed_neighbor_self(int radius, MAT::Node* root,bool all_node_reachable) {
    root->set_self_changed();
    root->set_ancestor_changed();
    root->set_descendent_changed();
    mark_changed_neighbor_descendent(radius, root,all_node_reachable);
    for (int radius_left=radius; radius_left>0; radius_left--) {
        root=root->parent;
        if (!root) {
            return;
        }
        if (root->get_self_changed()) {
            return;
        }
        /*if (all_node_reachable&&root->get_descendent_changed()) {
            return;
        }*/
        root->set_descendent_changed();
        mark_changed_neighbor_descendent(radius_left, root,all_node_reachable);
    }
}

//find src node to search
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                        std::vector<MAT::Node *> &output,
                        bool search_all_node, bool search_all_dir,int radius_in,MAT::Tree &tree) {
    auto start=std::chrono::steady_clock::now();
    unsigned int radius=abs(radius_in);
    output.clear();
    fprintf(stderr, "find max_level\n");
    auto max_level=tree.get_max_level();
    fprintf(stderr, "Max level %zu\n",max_level);
    bool all_node_reachable=(radius>=(2*max_level));
    if (!search_all_dir) {
        //mark surrounding of changed nodes
        FILE* fh=fopen("changed_nodes","w");
        std::vector<MAT::Node*> changed_nodes_ptr;
        changed_nodes_ptr.reserve(bfs_ordered_nodes.size());
        for (const auto node : bfs_ordered_nodes) {
            if (node->get_self_changed()) {
                changed_nodes_ptr.push_back(node);
                fprintf(fh, "%s\n", node->identifier.c_str());
            }
        }
        fclose(fh);
        fprintf(stderr, "%zu changed nodes \n",changed_nodes_ptr.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0,changed_nodes_ptr.size()),[radius,&changed_nodes_ptr,all_node_reachable](const tbb::blocked_range<size_t>& r) {
            for (size_t idx=r.begin(); idx<r.end(); idx++) {
                mark_changed_neighbor_self(radius, changed_nodes_ptr[idx],all_node_reachable);
            }
        });
    }
    int both_count=0;
    int upward_only_count=0;
    int downward_only_count=0;
    for(auto node:bfs_ordered_nodes) {
        auto upward=node->get_ancestor_changed();
        auto downward=node->get_descendent_changed();
        if (upward&&downward) {
            both_count++;
        } else if (upward) {
            upward_only_count++;
        } else if (downward) {
            downward_only_count++;
        }
    }
    fprintf(stderr, "Upward: %d,downward %d, both %d \n",upward_only_count,downward_only_count,both_count);
    if (search_all_node) {
        output=bfs_ordered_nodes;
        fprintf(stderr, "Search all nodes\n");
    } else {
        if (all_node_reachable) {
            output=bfs_ordered_nodes;
            fprintf(stderr, "Search all nodes\n");
        } else {
            output.reserve(bfs_ordered_nodes.size());
            for(auto node:bfs_ordered_nodes) {
                if (node->have_change_in_neighbor()) {
                    output.push_back(node);
                }
            }
        }
    }
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
