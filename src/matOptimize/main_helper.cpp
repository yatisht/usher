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
#include <unordered_set>
#include <random>
#include <utility>
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

//see whether within radius of root have nodes that have same mutation as root and have changed in previous iteration (or this is the first iteration)
static void check_changed_neighbor(int radius, MAT::Node* root, MAT::Node* exclude,bool& found,bool is_first,const std::unordered_set<int>& pos) {
    if (root->changed) {
        found=true;
        return;
    }
    if (!radius||found) {
        return;
    }
    if (root->parent&&root->parent!=exclude) {
        check_changed_neighbor(radius-1, root->parent, root, found,is_first,pos);
    }
    for(auto node:root->children) {
        if (found) {
            return;
        }
        if(node!=exclude) {
            check_changed_neighbor(radius-1, node, root, found,is_first,pos);
        }
    }
}

//find src node to search
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                        tbb::concurrent_vector<MAT::Node *> &output,
                        bool is_first, int radius) {
    std::vector<MAT::Node*> first_pass_nodes(bfs_ordered_nodes);
    auto start=std::chrono::steady_clock::now();
    std::shuffle(first_pass_nodes.begin(), first_pass_nodes.end(), rng);
    //first_pass_nodes.reserve(nodes_with_recurrent_mutations.size()*2);
    //add_ancestors(nodes_with_recurrent_mutations, pushed_nodes, first_pass_nodes);
    fprintf(stderr, "First pass nodes: %zu \n",first_pass_nodes.size());
    if (is_first) {
        output=tbb::concurrent_vector<MAT::Node*>(first_pass_nodes.begin(),first_pass_nodes.end());
    } else {
        output.reserve(first_pass_nodes.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0,first_pass_nodes.size()),[&first_pass_nodes,is_first,radius,&output](const tbb::blocked_range<size_t>& r) {
            for (size_t idx=r.begin(); idx<r.end(); idx++) {
                auto node=first_pass_nodes[idx];
                std::unordered_set<int> pos;
                pos.reserve(node->mutations.size()*4);
                for(const auto& mut:node->mutations) {
                    pos.insert(mut.get_position());
                }
                bool found=false;
                check_changed_neighbor(radius, node, nullptr, found, is_first, pos);
                if (found) {
                    output.push_back(node);
                }
            }
        });
    }
    fprintf(stderr, "Will search %f of nodes\n",(double)output.size()/(double)bfs_ordered_nodes.size());
    std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now()-start;
    fprintf(stderr, "Took %f s to find nodes to move\n",elapsed_seconds.count());
}
