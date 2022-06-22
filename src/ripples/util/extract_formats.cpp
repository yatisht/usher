#include "extract_formats.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

void get_parents(Mutation_Annotated_Tree::Tree *T,
                 std::unordered_set<std::string> &need_parents,
                 std::unordered_set<std::string> &all_nodes) {

    // Workflow run from "usher/scripts/recombination"
    FILE *node_to_parent_fp = fopen("filtering/data/nodeToParent.txt", "w");
    FILE *node_to_parent_no_underscore_fp =
        fopen("filtering/data/nodeToParent_no_underscore.txt", "w");

    if (node_to_parent_fp == NULL) {
        fprintf(stderr, "Error. ripplesUtils <MAT_name.pb>\n");
    }
    fprintf(node_to_parent_fp, "node\tparent\n");

    for (const auto &id : need_parents) {
        auto node = T->get_node(id);
        if (node == NULL) {
            continue;
        }
        auto parent = node->parent;
        if (parent == NULL) {
            continue;
        }
        std::string node_id = node->identifier;
        std::string parent_id = parent->identifier;

        // Insert parent node_id into collection of all relevant nodes
        all_nodes.insert(parent_id);

        // Format nodeToParent w/out "node_", temporary fix for simulations
        bool format_flag = true;
        if (format_flag == true) {
            if (node_id[0] == 'n') {
                node_id = node_id.substr(5, node_id.size());
            }
            if (parent_id[0] == 'n') {
                parent_id = parent_id.substr(5, parent_id.size());
            }
        }

        // Write to nodeToParent.txt file
        fprintf(node_to_parent_fp, "%s\t%s\n", node->identifier.c_str(),
                parent_id.c_str());

        // Write to nodeToParent_no_underscore.txt file
        fprintf(node_to_parent_no_underscore_fp, "%s\t%s\n", node_id.c_str(),
                parent_id.c_str());
    }
    fclose(node_to_parent_fp);
    fclose(node_to_parent_no_underscore_fp);
}
// This function is borrowed from src/matUtils/describe.cpp, except it simply
// outputs node ids in previously supported format eg.(1) instead of node_1
std::vector<std::string>
mutation_paths_no_label(MAT::Tree *T, std::vector<std::string> samples) {
    std::vector<std::string> mpaths;
    mpaths.push_back("sample_id\tpath_from_root");
    for (auto sample : samples) {
        std::string cpath = sample + "\t";
        auto root_to_sample = T->rsearch(sample, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        for (auto n : root_to_sample) {
            auto node_id = n->identifier;
            for (size_t i = 0; i < n->mutations.size(); i++) {
                cpath += n->mutations[i].get_string();
                if (i + 1 < n->mutations.size()) {
                    cpath += ",";
                }
            }
            if (n != root_to_sample.back()) {
                // note, for uncondensed samples with parsimony score 0,
                // this will leave a > at the end. That's basically fine.
                cpath += " (" + node_id.substr(5, node_id.size()) + ") > ";
            }
        }
        mpaths.push_back(cpath);
    }
    return mpaths;
}

// For now, write all samples paths (no samples explicitly named)
// Borrowed from src/matUtils/extract.cpp
void generate_sample_paths(MAT::Tree &T) {
    std::string sample_paths = "filtering/data/sample_paths.txt";
    std::vector<std::string> samples;
    samples = T.get_leaves_ids();

    // no samples given, subtree = tree
    MAT::Tree *subtree;
    subtree = &T;
    std::ofstream outfile(sample_paths);
    auto mpaths = mutation_paths_no_label(subtree, samples);
    for (auto mstr : mpaths) {
        outfile << mstr << "\n";
    }
    outfile.close();

    fprintf(stdout, "Finished generating sample paths file. Written to "
                    "filtering/data/sample_paths.txt.\n");
}

void leaves_per_node(MAT::Tree &T) {
    std::ofstream leaves_file("filtering/data/leaves.txt");

    auto dfs = T.depth_first_expansion();

    for (size_t i = 0; i < dfs.size(); ++i) {
        auto node_id = dfs[i]->identifier;
        if (node_id.find("node_") != std::string::npos) {
            leaves_file << node_id.substr(5, node_id.size()) << "\t";
        } else {
            leaves_file << node_id << "\t";
        }
        size_t num_leaves = T.get_num_leaves(dfs[i]);
        leaves_file << num_leaves << "\n";
    }
    leaves_file.close();
}
