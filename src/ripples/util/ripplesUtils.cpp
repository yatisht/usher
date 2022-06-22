#include "extract_formats.hpp"
#include <iostream>
#include <string>
#include <string_view>

int main(int argc, char **argv) {

    std::string input_mat_filename;
    // Workflow run from "usher/scripts/recombination"
    const std::string &combined_pvals_filename =
        "filtering/data/combinedCatOnlyBestWithPVals.txt";

    if (argc == 2) {
        input_mat_filename = std::string{argv[1]};
    } else {
        fprintf(stderr, "Error. ripplesUtils <MAT_name.pb>\n");
        exit(1);
    }

    fprintf(stdout, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    // Load input MAT and uncondense tree
    MAT::Tree T;
    if (input_mat_filename.find(".pb\0") != std::string::npos) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
    } else {
        fprintf(
            stderr,
            "ERROR: Input file ending not recognized. Must be .json or .pb\n");
        exit(1);
    }
    fprintf(stdout, "Completed loading MAT, getting trios now.\n");
    fprintf(stdout, "Generating sample paths file.\n");
    generate_sample_paths(T);

    // Ouputs two files:  allRelevantNodeNames.txt and nodeToParent.txt
    get_trios(T, combined_pvals_filename);

		// Generate "leaves.txt" file, same as "optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt" 
		// The file contains the following information on each line for all nodes in tree:  node_id \t number of leaves
    fprintf(stdout, "Calculating number of leaves per node in the tree. Written to 'filtering/data/leaves.txt'.\n");
    leaves_per_node(T);

    return 0;
}
