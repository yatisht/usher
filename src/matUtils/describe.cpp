#include "describe.hpp"

std::vector<std::string> mutation_paths(MAT::Tree* T, std::vector<std::string> samples) {
    std::vector<std::string> mpaths;
    mpaths.push_back("sample_id\tpath_from_root");
    for (auto sample: samples) {
        std::string cpath = sample + "\t";
        auto root_to_sample = T->rsearch(sample, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        for (auto n: root_to_sample) {
            for (size_t i=0; i<n->mutations.size(); i++) {
                cpath += n->mutations[i].get_string();
                if (i+1 < n->mutations.size()) {
                    cpath += ",";
                }
            }
            if (n != root_to_sample.back()) {
                //note, for uncondensed samples with parsimony score 0,
                //this will leave a > at the end. That's basically fine.
                cpath += " (" + n->identifier + ") > ";
            }
        }
        mpaths.push_back(cpath);
    }
    return mpaths;
}

std::vector<std::string> clade_paths(MAT::Tree* T) {
    //get the set of clade path strings for printing
    //similar to the above, construct a series of strings to be printed or redirected later on
    std::vector<std::string> clpaths;
    clpaths.push_back("clade\troot_id\tfrom_tree_root\n");
    //do a breadth-first search
    //clades are annotated only at the root, so when we see an annotation, add it to the list.

    auto bfs = T->breadth_first_expansion();
    for (auto n: bfs) {
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                std::string curpath;
                //record the name of the clade
                curpath += ann + "\t";
                curpath += n->identifier + "\t";
                //get the ancestral mutations back to the root
                std::string root = "";
                auto root_to_sample = T->rsearch(n->identifier, true);
                std::reverse(root_to_sample.begin(), root_to_sample.end());
                for (auto an: root_to_sample) {
                  for (size_t i=0; i<an->mutations.size(); i++) {
                      root += an->mutations[i].get_string();
                      if (i+1 < an->mutations.size()) {
                          root += ",";
                      }
                  }
                  if (an != root_to_sample.back()) {
                      root += " > ";
                  }
                }
                //save values to the string and save the string
                curpath += root;
                clpaths.push_back(curpath + "\n");
            }
        }
    }
    return clpaths;
}

std::vector<std::string> all_nodes_paths(MAT::Tree* T) {
    std::vector<std::string> dfs_strings;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        std::string node_path = n->identifier + ": ";
        for (size_t i=0; i<n->mutations.size(); i++) {
            node_path += n->mutations[i].get_string();
            if (i+1 < n->mutations.size()) {
                node_path += ",";
            }
        }
        dfs_strings.push_back(node_path);
    }
    return dfs_strings;
}
