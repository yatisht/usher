#include "describe.hpp"

std::vector<std::string> mutation_paths(const MAT::Tree& T, std::vector<std::string> samples) {
    std::vector<std::string> mpaths;
    for (auto sample: samples) {
        std::string cpath = sample + ": ";
        auto root_to_sample = T.rsearch(sample, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        //fprintf(stdout, "%s: ", sample.c_str());
        for (auto n: root_to_sample) {
            if (!n->is_root()) {
                for (size_t i=0; i<n->mutations.size(); i++) {
                    cpath += n->mutations[i].get_string();
                    //fprintf(stdout, "%s", n->mutations[i].get_string().c_str());
                    if (i+1 < n->mutations.size()) {
                        cpath += ",";
                        //fprintf(stdout, ",");
                    }
                }
                // //some samples have no internal mutations, which is really weird and concering
                // if (n->mutations.size() == 0 && n->identifier != sample) {
                //     //there are definitely internal nodes that have no mutations
                //     //most commonly its the root (node 1) that has no mutations on a given path.
                //     fprintf(stderr, "DEBUG: Node %s", n->identifier.c_str());
                //     fprintf(stderr, " on the path of %s has no mutations.", sample.c_str());
                //     fprintf(stderr, " Level %ld.", n->level);
                //     fprintf(stderr, " Branch length %f.", n->branch_length);
                //     fprintf(stderr, " Parent is %s.", n->parent->identifier.c_str());
                //     fprintf(stderr, " %ld children.\n", n->children.size());
                // }
                if (n != root_to_sample.back()) {
                    //note, for uncondensed samples with parsimony score 0,
                    //this will leave a > at the end. That's basically fine.
                    cpath += " > ";
                    //fprintf(stdout, " > ");
                }
            }
        }
        mpaths.push_back(cpath);
        //fprintf(stdout, "\n");
    }
    return mpaths;
}

std::vector<std::string> clade_paths(MAT::Tree T, std::vector<std::string> clades) {
    //get the set of clade path strings for printing
    //similar to the above, construct a series of strings to be printed or redirected later on
    std::vector<std::string> clpaths;
    //NOTE: setting aside functionality of getting unique mutations for clarity for the moment, may reimplement later.
    clpaths.push_back("clade\troot_id\tfrom_tree_root");
    //clpaths.push_back("clade\troot_id\tspecific\tfrom_tree_root");
    //do a breadth-first search
    //the first time a clade that is in clades is encountered, that's the root;
    //get the path of mutations to that root (rsearch), save the unique mutations + that path
    //unique mutations being ones that occurred in the clade root, and the path being all mutations from that root back to the tree root
    //then continue. if a clade has already been encountered in the breadth first, its
    //not clade root, and should be skipped.
    std::unordered_set<std::string> clades_seen;

    auto dfs = T.breadth_first_expansion();
    for (auto n: dfs) {
        std::string curpath;
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                //if its one of our target clades and it hasn't been seen...
                if (std::find(clades.begin(), clades.end(), ann) != clades.end() && clades_seen.find(ann) == clades_seen.end()) {
                    //record the name of the clade
                    curpath += ann + "\t";
                    curpath += n->identifier + "\t";
                    // //get its own mutations, if there are any
                    // std::string unique = "";
                    // for (size_t i=0; i<n->mutations.size(); i++) {
                    //     unique += n->mutations[i].get_string();
                    //     if (i+1 < n->mutations.size()) {
                    //         unique += ",";
                    //     }
                    // }
                    // curpath += unique + "\t";
                    //get the ancestral mutations back to the root
                    std::string root = "";
                    auto root_to_sample = T.rsearch(n->identifier, true);
                    std::reverse(root_to_sample.begin(), root_to_sample.end());
                    for (auto an: root_to_sample) {
                        //skip the root.
                        if (!an->is_root()) {
                            for (size_t i=0; i<an->mutations.size(); i++) {
                                root += an->mutations[i].get_string();
                                if (i+1 < an->mutations.size()) {
                                    root += ",";
                                }
                            }
                            // //some samples have no internal mutations, which is really weird and concering
                            // if (an->mutations.size() == 0) {
                            //     //there are definitely internal nodes that have no mutations
                            //     //most commonly its the root (node 1) that has no mutations on a given path.
                            //     fprintf(stderr, "DEBUG: Node %s", an->identifier.c_str());
                            //     fprintf(stderr, " on the path of %s has no mutations.", n->identifier.c_str());
                            //     fprintf(stderr, " Level %ld.", an->level);
                            //     fprintf(stderr, " Branch length %f.\n", an->branch_length);
                            // }
                            if (an != root_to_sample.back()) {
                                root += " > ";
                            }
                        }
                    }
                    //save values to the string, mark the clade as seen, and save the string
                    curpath += root;
                    clades_seen.insert(ann);
                    clpaths.push_back(curpath + "\n");
                }
            }
        }
    }
    return clpaths;
}