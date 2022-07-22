#include "init_pipeline.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <iostream>
#include <random>

// Same preorder traversal as Chronumental performs to map
// from ripples node ids -> Chronumental node ids
void preorder_traversal(MAT::Tree &T) {
    MAT::Node *root = T.root;
    std::stack<MAT::Node *> s;
    std::vector<MAT::Node *> preorder;
    if (root == NULL) {
        std::cout << "ERROR: Empty tree!"
                  << "\n";
        exit(1);
    }
    s.push(root);
    while (!s.empty()) {
        auto node = s.top();
        s.pop();
        preorder.push_back(node);
        for (auto &child : node->children) {
            s.push(child);
        }
    }
    auto dfs = T.depth_first_expansion();
    if (dfs.size() != preorder.size()) {
        std::cout << "ERROR: Traversal sizes not matching."
                  << "\n";
        exit(1);
    }
    std::ofstream node_file("ripples_to_chron_ids.txt");
    node_file << "MAT_node_id"
              << "\t"
              << "chronumental_node_id"
              << "\n";
    for (size_t i = 0; i < dfs.size(); ++i) {
        node_file << dfs[i]->identifier << "\t";
        node_file << preorder[i]->identifier << "\n";
    }

    node_file.close();
}

// Temporary borrowed from ripples code to find number of long branches in the tree
// that will be considered for searching.
int find_long_branches(MAT::Tree &T, uint32_t branch_len, uint32_t num_descendants) {
    auto bfs = T.breadth_first_expansion();
    static tbb::affinity_partitioner ap;
    tbb::concurrent_unordered_set<std::string> nodes_to_consider;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs.size()),
    [&](const tbb::blocked_range<size_t> r) {
        for (size_t i = r.begin(); i < r.end(); ++i) {
            auto n = bfs[i];
            if (n == T.root) {
                continue;
            }
            if (n->mutations.size() >= branch_len) {
                if (T.get_num_leaves(n) >= num_descendants) {
                    nodes_to_consider.insert(n->identifier);
                }
            }
        }
    },
    ap);
    std::vector<std::string> nodes_to_consider_vec;
    for (const auto &elem : nodes_to_consider) {
        nodes_to_consider_vec.emplace_back(elem);
    }
    std::sort(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end());
    std::shuffle(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end(),
                 std::default_random_engine(0));

    return nodes_to_consider.size();
}
