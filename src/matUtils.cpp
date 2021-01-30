#define matToVCF
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
#include "mutation_annotated_tree.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

po::variables_map check_options(int argc, char** argv) {
    // Check command line options and return variable map.
    po::options_description desc{"Options"};
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file to mask [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->required(),
         "Output masked mutation-annotated tree file [REQUIRED]")
        ("restricted-samples,s", po::value<std::string>()->required(),
         "[REQUIRED]")
        ("help,h", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(p)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

int main(int argc, char** argv) {

    // Command line options
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string samples_filename = vm["restricted-samples"].as<std::string>();

    // Load input MAT and uncondense tree
    auto T = MAT::load_mutation_annotated_tree(input_mat_filename);
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }

    // Load restricted sample names from the input file and add it to the set
    std::unordered_set<std::string> restricted_samples;
    std::ifstream infile(samples_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the restricted sampels file: %s!\n", samples_filename.c_str());
        exit(1);
    }
    std::string sample;
    while (std::getline(infile, sample)) {
        if (T.get_node(sample) == NULL) {
            fprintf(stderr, "ERROR: Sample %s missing in input MAT!\n", sample.c_str());
            exit(1);
        }
        restricted_samples.insert(std::move(sample));
    }

    // Set of nodes rooted at restricted samples
    std::unordered_set<MAT::Node*> restricted_roots;
    std::unordered_map<std::string, bool> visited;
    for (auto s: restricted_samples) {
        visited[s] = false;
    }
    for (auto cn: T.breadth_first_expansion()) { 
        auto s = cn->identifier;
        if (restricted_samples.find(s) == restricted_samples.end()) {
            continue;
        }
        if (visited[s]) {
            continue;
        }
        auto curr_node = T.get_node(s);
        for (auto n: T.rsearch(s)) {
            bool found_unrestricted = false;
            for (auto l: T.get_leaves_ids(n->identifier)) {
                if (restricted_samples.find(l)  == restricted_samples.end()) {
                    found_unrestricted = true;
                    break;
                }
            }
            if (!found_unrestricted) {
                for (auto l: T.get_leaves_ids(n->identifier)) {
                    visited[l] = true;
                }
                curr_node = n;
                break;
            }
        }
        restricted_roots.insert(curr_node);
    }

    fprintf(stderr, "Restricted roots size: %zu\n\n", restricted_roots.size());

    // Map to store number of occurences of a mutation in the tree
    std::unordered_map<std::string, int> mutations_counts;
    for (auto n: T.depth_first_expansion()) {
        for (auto mut: n->mutations) {
            if (mut.is_masked()) {
                continue;
            }
            auto mut_string = mut.get_string();
            if (mutations_counts.find(mut_string) == mutations_counts.end()) {
                mutations_counts[mut_string] = 1;
            }
            else {
                mutations_counts[mut_string] += 1;
            }
        }
    }
    
    // Reduce mutation counts for mutations in subtrees rooted at 
    // restricted_roots. Mutations specific to restricted samples 
    // will now be set to 0. 
    for (auto r: restricted_roots) {
        //fprintf(stdout, "At restricted root %s\n", r->identifier.c_str());
        for (auto n: T.depth_first_expansion(r)) {
            for (auto mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                mutations_counts[mut_string] -= 1;
            }
        }
    }

    for (auto r: restricted_roots) {
        for (auto n: T.depth_first_expansion(r)) {
            for (auto& mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                if (mutations_counts[mut_string] == 0) {
                    fprintf(stderr, "Masking mutation %s at node %s\n", mut_string.c_str(), n->identifier.c_str());
                    mut.position = -1;
                    mut.ref_nuc = 0;
                    mut.par_nuc = 0;
                    mut.mut_nuc = 0;
                }
            }
        }
    }

    // Store final MAT to output file 
    MAT::save_mutation_annotated_tree(T, output_mat_filename);

    return 0;
}

