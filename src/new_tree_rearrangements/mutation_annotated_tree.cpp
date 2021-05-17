#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <iostream>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <stack>
#include <queue>
#include <vector>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8
using Mutation_Annotated_Tree::Node;
tbb::concurrent_unordered_map<std::string, uint8_t>  Mutation_Annotated_Tree::Mutation::chromosome_map;
std::vector<std::string>  Mutation_Annotated_Tree::Mutation::chromosomes;
std::mutex Mutation_Annotated_Tree::Mutation::ref_lock;
std::vector<nuc_one_hot> Mutation_Annotated_Tree::Mutation::refs;



/* === Tree === */
        
std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::get_leaves(std::string nid) {
    std::vector<Node*> leaves;
    if (nid == "") {
        if (root == NULL) {
            return leaves;
        }
        nid = root->identifier;
    }
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves;
}

std::vector<std::string> Mutation_Annotated_Tree::Tree::get_leaves_ids(std::string nid) {
    std::vector<std::string> leaves_ids;
    if (nid == "") {
        if (root == NULL) {
            return leaves_ids;
        }
        nid = root->identifier;
    }
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves_ids.push_back(curr_node->identifier);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves_ids;
}

size_t Mutation_Annotated_Tree::Tree::get_num_leaves(Node* node) {
    if (node == NULL) {
        node = root;
    }

    if (node->is_leaf()) {
        return 1;
    }
    size_t num_leaves = 0;
    for (auto c: node->children) {
        num_leaves += get_num_leaves(c);
    }
    return num_leaves;
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::rsearch (const std::string& nid, bool include_self) const {
    std::vector<Node*> ancestors;
    Node* node = get_node(nid);
    if (node==NULL) {
        return ancestors;
    }    
    if (include_self) {
        ancestors.push_back(node);
    }
    while (node->parent != NULL) {
        ancestors.push_back(node->parent);
        node = node->parent;
    }
    return ancestors;
}


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

static void depth_first_expansion_helper(Mutation_Annotated_Tree::Node* node, std::vector<Mutation_Annotated_Tree::Node*>& vec, size_t& index) {
    #ifdef DETAIL_DEBUG_NO_LOOP
    assert(std::find(vec.begin(),vec.end(),node)==vec.end());
    #endif
    vec.push_back(node);
    assert(vec.size()-1==index);
    node->dfs_index=index;
    index++;
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec,index);
    }
    node->dfs_end_index=index;
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
    depth_first_expansion_helper(node, traversal,index);
    return traversal;
}

size_t Mutation_Annotated_Tree::Tree::get_parsimony_score() {
    size_t score = 0;
    auto dfs = depth_first_expansion();
    for (auto n: dfs) {
        for (const auto& mut:n->mutations){
            score+=mut.is_valid();
        }
    }
    return score;
}

void Mutation_Annotated_Tree::Tree::condense_leaves(std::vector<std::string> missing_samples) {
    if (condensed_nodes.size() > 0) {
        fprintf(stderr, "WARNING: tree contains condensed nodes. It may be condensed already!\n");
    }

    auto tree_leaves = get_leaves_ids();
    for (auto l1_id: tree_leaves) {
        std::vector<Node*> polytomy_nodes;

        auto l1 = get_node(l1_id);
        if (l1 == NULL) {
            continue;
        }
        if (std::find(missing_samples.begin(), missing_samples.end(), l1->identifier) != missing_samples.end()) {
            continue;
        }
        if (l1->mutations.size() > 0) {
            continue;
        }

        for (auto l2: l1->parent->children) {
                if (std::find(missing_samples.begin(), missing_samples.end(), l2->identifier) != missing_samples.end()) {
                    continue;
                }
            if (l2->is_leaf() && (get_node(l2->identifier) != NULL) && (l2->mutations.size() == 0)) {
                polytomy_nodes.push_back(l2);
            }
        }

        if (polytomy_nodes.size() > 1) {
            std::string new_node_name = "node_" + std::to_string(1+condensed_nodes.size()) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
            
            auto curr_node = get_node(l1->identifier);
            auto new_node = create_node(new_node_name, curr_node->parent, l1->branch_length);

            new_node->clear_mutations();
            
            condensed_nodes[new_node_name] = std::vector<std::string>(polytomy_nodes.size());

            for (size_t it = 0; it < polytomy_nodes.size(); it++) {
                condensed_nodes[new_node_name][it] = polytomy_nodes[it]->identifier;
                remove_node(polytomy_nodes[it]->identifier, false);
            }
        }
    }
}

void Mutation_Annotated_Tree::Tree::uncondense_leaves() {
    for (size_t it = 0; it < condensed_nodes.size(); it++) {
        auto cn = condensed_nodes.begin();
        std::advance(cn, it);

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
    condensed_leaves.clear();
}
// Merge nodes that have no mutations comparing to parent into parent node 
void Mutation_Annotated_Tree::Tree::collapse_tree() {
    auto bfs = breadth_first_expansion();

    for (size_t idx = 1; idx < bfs.size(); idx++) {
        auto node = bfs[idx];
        auto mutations = node->mutations;        
        if (mutations.size() == 0) {
            auto parent = node->parent;
            auto children = node->children;
            for (auto child: children) {
                move_node(child->identifier, parent->identifier);
            }
        }
        //If internal node has one child, the child can be moved up one level
        else if (node->children.size() == 1) {
            auto child = node->children.front();
            auto parent = node->parent;
            for (auto m: mutations) {
                child->add_mutation(m);
            }
            move_node(child->identifier, parent->identifier);
        }
    }
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_tree_copy(const Mutation_Annotated_Tree::Tree& tree, const std::string& identifier) {
    TIMEIT();
    auto root = tree.root;
    if (identifier != "") {
        root = tree.get_node(identifier);
    }
    
    Tree copy = create_tree_from_newick_string (get_newick_string(tree, root, true, true));

    std::vector<Node*> dfs1;
    std::vector<Node*> dfs2;

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, 2),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
              if (k==0) {
                dfs1 = tree.depth_first_expansion(root);
              }
              else {
                dfs2 = copy.depth_first_expansion();
              }
            }
            }, ap);


    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs1.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
              auto n1 = dfs1[k];
              auto n2 = dfs2[k];
              n2->clade_annotations.resize(n1->clade_annotations.size());
              for (size_t i=0; i<n1->clade_annotations.size(); i++) {
                 n2->clade_annotations[i] = n1->clade_annotations[i];
              }
              for (auto m: n1->mutations) {
                Mutation m2 = m.copy();
                n2->add_mutation(m2);
                }
              }
            }, ap);

    size_t num_condensed_nodes = static_cast<size_t>(tree.condensed_nodes.size());
    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_condensed_nodes),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t idx = r.begin(); idx < r.end(); idx++) {
               auto cn = tree.condensed_nodes.begin(); 
               std::advance(cn, idx);
               copy.condensed_nodes.insert(std::pair<std::string, std::vector<std::string>>(cn->first, std::vector<std::string>(cn->second.size())));
               for (size_t k = 0; k < cn->second.size(); k++) {
                  copy.condensed_nodes[cn->first][k] = cn->second[k];
                  copy.condensed_leaves.insert(cn->second[k]);
               }
            }
    }, ap);

    return copy;
}

// Extract the subtree consisting of the specified set of samples. This routine
// maintains the internal node names of the input tree. Mutations are copied
// from the tree such that the path of mutations from root to the sample is
// same as the original tree.
Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_subtree (const Mutation_Annotated_Tree::Tree& tree, const std::vector<std::string>& samples) {
    TIMEIT();
    Tree subtree;

    // Set of leaf and internal nodes corresponding to the subtree
    tbb::concurrent_unordered_set<Node*> subtree_nodes;
    // Maintain a set of all ancestors of a sample for each sample
    std::vector<tbb::concurrent_unordered_set<Node*>> all_ancestors(samples.size());

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
               subtree_nodes.insert(tree.get_node(samples[k]));
               for (auto anc: tree.rsearch(samples[k], true)) {
                   all_ancestors[k].insert(anc);
               }
            }
    }, ap);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
               for (size_t j=i+1; j<samples.size(); ++j){
                   for (auto anc: tree.rsearch(samples[i], true)) {
                      if (all_ancestors[j].find(anc) != all_ancestors[j].end()) {
                         subtree_nodes.insert(anc);
                         break;
                      }
                   }
               }
            }
    }, ap);
    
    auto dfs = tree.depth_first_expansion();
    size_t num_annotations = tree.get_num_annotations();

    std::stack<Node*> last_subtree_node;
    for (auto n: dfs) {
        // If the node is in subtree_nodes, it should be added to the subtree
        if (subtree_nodes.find(n) != subtree_nodes.end()) {
            Node* subtree_parent = NULL;
            if (last_subtree_node.size() > 0) {
                while (!tree.is_ancestor(last_subtree_node.top()->identifier, n->identifier)) {
                    last_subtree_node.pop();
                }
                subtree_parent = last_subtree_node.top();
            }
            // Add as root of the subtree
            if (subtree_parent == NULL) {
                // for root node, need to size the annotations vector
                Node* new_node = subtree.create_node(n->identifier, -1.0, num_annotations);
                // need to assign any clade annotations which would belong to that root as well
                for (size_t k = 0; k < num_annotations; k++) {
                    if (n->clade_annotations[k] != "") {
                        new_node->clade_annotations[k] = n->clade_annotations[k];
                    }
                }
                std::vector<Node*> root_to_node = tree.rsearch(n->identifier, true); 
                std::reverse(root_to_node.begin(), root_to_node.end());
                root_to_node.emplace_back(n);

                for (auto curr: root_to_node) {
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            // Add to the parent identified
            else {
                Node* new_node = subtree.create_node(n->identifier, subtree_parent->identifier);

                auto par_to_node = tree.rsearch(n->identifier, true);
                std::reverse(par_to_node.begin(), par_to_node.end());
                par_to_node.erase(par_to_node.begin(), std::find(par_to_node.begin(), par_to_node.end(), subtree_parent)+1);


                for (auto curr: par_to_node) {
                    for (size_t k = 0; k < num_annotations; k++) {
                        if (curr->clade_annotations[k] != "") {
                            new_node->clade_annotations[k] = curr->clade_annotations[k];
                        }
                    }
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            last_subtree_node.push(n);
        }
    }

    subtree.curr_internal_node = tree.curr_internal_node;

    return subtree;
}

void Node::delete_this(){
    for(Node* n:children){
        n->delete_this();
    }
    delete this;
}
void Mutation_Annotated_Tree::Tree::delete_nodes(){
    root->delete_this();
}