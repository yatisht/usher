#include "extract_formats.hpp"
#include "text_parser.hpp"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

void get_trios(MAT::Tree T, std::string filepath) {

    // Read entire dataset into memory
    text_parser rec(filepath); // combinedCatOnlyBestWithPVals.txt

    std::unordered_set<std::string> all_nodes;
    std::unordered_set<std::string> needParent;

    // Skip over header line (starting with #)
    rec.next_line();
    std::string_view yes{"y"};

    for (; !rec.done(); rec.next_line()) {

        // Add "node_" if id only, if "node_if" already or sample name
        // don't change
        auto recomb = std::string{rec.get_value(0)};
        if (std::isdigit(recomb.at(0)) == 1) {
            recomb = "node_" + recomb;
        }
        auto donor = std::string{rec.get_value(3)};
        if (std::isdigit(donor.at(0)) == 1) {
            donor = "node_" + donor;
        }
        auto acceptor = std::string{rec.get_value(6)};
        if (std::isdigit(acceptor.at(0)) == 1) {
            acceptor = "node_" + acceptor;
        }

        // Get the recomb, donor, acceptor trios from each line
        all_nodes.insert(recomb);
        all_nodes.insert(donor);
        all_nodes.insert(acceptor);

        // Check if donor is placed as sibling
        if (rec.get_value(4) == yes) {
            needParent.insert(donor);
        }
        // Check if acceptor is placed as sibling
        if (rec.get_value(7) == yes) {
            needParent.insert(acceptor);
        }
    }
    // Call get_parents to retrieve the parents for all the donors/acceptors
    // placed as siblings
    get_parents(&T, needParent, all_nodes);

    // Create allRelevantNodeNames.txt file
    // Workflow run from "usher/scripts/recombination"
    FILE *allRelevantNodeNames_fp =
        fopen("filtering/data/allRelevantNodeNames.txt", "w");
    for (const auto &n : all_nodes) {
        fprintf(allRelevantNodeNames_fp, "%s\n", n.c_str());
    }
    fclose(allRelevantNodeNames_fp);
}


inline uint64_t str_to_uint64(std::string_view str) noexcept {
    uint64_t num = 0;
    for (auto c : str)
        (num *= 10) += c - '0';
    return num;
}

