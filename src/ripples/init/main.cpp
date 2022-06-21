#include "init_pipeline.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <random>
#include <stdio.h>
#include <utility>

namespace po = boost::program_options;

// Pre-ripples script to get general run information needed to configure parallel jobs.
int main(int argc, char **argv) {

    po::options_description desc("optimize options");
    desc.add_options()("input-mat,i", po::value<std::string>()->required(),
                       "Input mutation-annotated tree file [REQUIRED].")(
        "branch-length,l", po::value<uint32_t>()->default_value(3),
        "Minimum length of the branch to consider to recombination events")(
        "num-descendants,n", po::value<uint32_t>()->default_value(2),
        "Minimum number of leaves that node should have to be considered for "
        "recombination.");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless
        // the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    uint32_t branch_len = vm["branch-length"].as<uint32_t>();
    uint32_t num_descendants = vm["num-descendants"].as<uint32_t>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    // Output text file for mapping UShER node_ids to Chronumental node_ids
    preorder_traversal(T);

    size_t long_branches = find_long_branches(T, branch_len, num_descendants);
		std::cout << long_branches << "\n";
    return 0;
}
