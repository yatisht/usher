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

void sort_if_necessary(std::vector<MAT::Mutation>& mutations) {
    // If mutations are not already sorted by position, sort them.
    bool is_sorted = true;
    int prev_pos = 0;
    for (auto mut: mutations) {
        if (mut.position < prev_pos) {
            is_sorted = false;
            break;
        }
    }
    if (! is_sorted) {
        std::sort(mutations.begin(), mutations.end());
    }
}

std::vector<MAT::Mutation> add_mutations(const std::vector<MAT::Mutation>& parent_muts,
                                         const std::vector<MAT::Mutation>& node_muts) {
    // Return a new vector that includes both parent_muts and node_muts, but collapsing multiple
    // mutations at the same position.  Inputs must be sorted by position.  Output is sorted.
    std::vector<MAT::Mutation> combined_muts;
    if (parent_muts.size() == 0) {
        combined_muts = node_muts;
    } else if (node_muts.size() == 0) {
        combined_muts = parent_muts;
    } else {
        size_t px = 0;
        for (auto n: node_muts) {
            while (parent_muts[px].position < n.position && px < parent_muts.size()) {
                combined_muts.emplace_back(parent_muts[px]);
                px++;
            }
            if (px < parent_muts.size()) {
                if (parent_muts[px].position == n.position) {
                    if (n.mut_nuc == parent_muts[px].par_nuc) {
                        // They cancel each other out; don't add either to combined_muts.
                    } else {
                        MAT::Mutation mut;
                        mut.par_nuc = parent_muts[px].par_nuc;
                        mut.mut_nuc = n.mut_nuc;
                        combined_muts.emplace_back(mut);
                    }
                    px++;
                } else {
                    combined_muts.emplace_back(n);
                }
            } else {
                combined_muts.emplace_back(n);
            }
        }
        while (px < parent_muts.size()) {
            combined_muts.emplace_back(parent_muts[px]);
            px++;
        }
    }
    
    return combined_muts;
}

size_t count_reversions(const std::vector<MAT::Mutation>& clade_muts,
                        const std::vector<MAT::Mutation>& node_muts) {
    // Return the number of reversions to reference from clade_muts in node_muts.
    // Inputs must be sorted by position.
    size_t rev_count = 0;
    if (clade_muts.size() > 0 && node_muts.size() > 0) {
        size_t cx = 0;
        for (auto n: node_muts) {
            while (clade_muts[cx].position < n.position && cx < clade_muts.size()) {
                cx++;
            }
              if (cx < clade_muts.size() &&
                  clade_muts[cx].position == n.position &&
                  n.mut_nuc == clade_muts[cx].par_nuc) {
                  rev_count++;
              }
        }
    }
    return rev_count;
}

void print_node_stats_r(Mutation_Annotated_Tree::Node* node, size_t& leaf_count, size_t& mut_count,
                        const std::vector<MAT::Mutation>& parent_clade_muts,
                        const std::vector<MAT::Mutation>& parent_muts, size_t parent_rev_count,
                        std::ofstream& outfile) {
    // Recursively descend from node, counting number of reversions since last annotated clade on
    // the way down and tallying up leaf counts and total mut counts of descendants on the way up.
    // Print stats for each leaf and internal node.
    sort_if_necessary(node->mutations);
    const std::vector<MAT::Mutation> my_muts = add_mutations(parent_muts, node->mutations);
    bool is_clade_root = std::any_of(node->clade_annotations.begin(), node->clade_annotations.end(),
                                     [](std::string& clade){ return clade != ""; });
    size_t rev_count = is_clade_root ? 0 :
      (parent_rev_count + count_reversions(parent_clade_muts, node->mutations));
    if (node->children.size() > 0) {
        size_t leaf_count_total = 0, mut_count_total = node->mutations.size();
        const std::vector<MAT::Mutation>& clade_muts = is_clade_root ? my_muts : parent_clade_muts;
        for (auto child: node->children) {
            size_t leaf_count_child, mut_count_child;
            print_node_stats_r(child, leaf_count_child, mut_count_child,
                               clade_muts, my_muts, rev_count, outfile);
            leaf_count_total += leaf_count_child;
            mut_count_total += mut_count_child;
        }
        outfile << node->identifier << "\t" << leaf_count_total << "\t" << mut_count_total <<
          "\t" << mut_count_total / (double)leaf_count_total << "\t" << rev_count << "\n";
        leaf_count = leaf_count_total;
        mut_count = mut_count_total;
    } else {
        leaf_count = 1;
        mut_count = node->mutations.size();
        outfile << node->identifier << "\t" << leaf_count << "\t" << mut_count <<
          "\t" << mut_count << "\t" << rev_count << "\n";
    }
}

void print_node_stats(Mutation_Annotated_Tree::Node* node, size_t& leaf_count, size_t& mut_count,
                      std::ofstream& outfile) {
    // Print out columns that can be used to determine the "aspect ratio" of each internal node's
    // branch displayed as a rectangular tree: "tall and narrow" branches have many leaves with
    // relatively few mutations, while "short and wide" branches have fewer leaves with relatively
    // many mutations.  This can be used as a sort of branch-level quality filter.
    // Also print out the number of reversions relative to the annotated clade for both internal
    // nodes and leaves.
    outfile << "node\tleaf_count\tmut_count\tmut_density\trev_from_lineage\n";
    std::vector<MAT::Mutation> no_muts;
    print_node_stats_r(node, leaf_count, mut_count, no_muts, no_muts, 0, outfile);
}
