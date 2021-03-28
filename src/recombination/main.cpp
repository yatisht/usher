#include <time.h>  
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
#include "tbb/concurrent_unordered_set.h"
#include "../usher_graph.hpp"

namespace po = boost::program_options;

Timer timer;

po::variables_map check_options(int argc, char** argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::options_description desc("optimize options");
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file to optimize [REQUIRED].")
        ("branch-length,l", po::value<uint32_t>()->default_value(4), \
         "Minimum length of the branch to consider to recombination events")
        ("max-coordinate-range,R", po::value<int>()->default_value(1e7), \
         "Maximum range of the genomic coordinates of the mutations on the branch")
        ("min-coordinate-range,r", po::value<int>()->default_value(1e3), \
         "Minimum range of the genomic coordinates of the mutations on the branch")
        ("num-descendants,n", po::value<uint32_t>()->default_value(10), \
         "Minimum number of leaves that node should have to be considered for recombination.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
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
        // Return with error code 1 unless
        // the user specifies help
            if (vm.count("help"))
                exit(0);
            else
                exit(1);
    }
    return vm;
}

struct Pruned_Sample {
    std::string sample_name;
    std::vector<MAT::Mutation> sample_mutations;
    std::unordered_set<uint32_t> positions;

    // Assumes mutations are added in reverse chrono order
    void add_mutation (MAT::Mutation mut) {
        // If not reversal to reference allele 
        if ((mut.ref_nuc != mut.mut_nuc) && (positions.find(mut.position) == positions.end())) {
            auto iter = std::lower_bound(sample_mutations.begin(), sample_mutations.end(), mut);
            mut.par_nuc = mut.ref_nuc;
            sample_mutations.insert(iter, mut);
        }
        positions.insert(mut.position);
    }

    Pruned_Sample (std::string name) {
        sample_name = name;
        sample_mutations.clear();
        positions.clear();
    }
};

int main(int argc, char** argv) {
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    uint32_t branch_len = vm["branch-length"].as<uint32_t>();
    int max_range = vm["max-coordinate-range"].as<int>();
    int min_range = vm["min-coordinate-range"].as<int>();
    uint32_t num_descendants = vm["num-descendants"].as<uint32_t>();
    
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);
    srand (time(NULL));

    static tbb::affinity_partitioner ap;

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    
    fprintf(stderr, "Finding the branches with number of mutations equal to or exceeding %u.\n", branch_len);
    auto bfs = T.breadth_first_expansion();
    
    tbb::concurrent_unordered_set<std::string> nodes_to_consider;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, bfs.size()),
            [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
               auto n = bfs[i];
               if (n == T.root) {
                   continue;
               }
               if (n->mutations.size() >= branch_len) {
                   int start_range_high = n->mutations[0].position;
                   int end_range_low = n->mutations[0].position;
                   for (auto m: n->mutations) {
                       if (m.position < start_range_high) {
                           start_range_high = m.position;
                       }
                       if (m.position > end_range_low) {
                           end_range_low = m.position;
                       }
                   }
                   if ((end_range_low-start_range_high >= min_range) && (end_range_low-start_range_high <= max_range) && (T.get_num_leaves(n) >= num_descendants)) {
                       nodes_to_consider.insert(n->identifier);
                   }
               }
            }
        }, ap);
    
    fprintf(stderr, "Found %zu long branches\n\n", nodes_to_consider.size());

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    
    timer.Start();
    
    fprintf(stderr, "Running placement individually for %zu branches to identify potential recombination events.\n", nodes_to_consider.size());

    size_t num_done = 0;
    for (auto nid_to_consider: nodes_to_consider) {
        Pruned_Sample pruned_sample(nid_to_consider);

        // Find mutations on the node to prune
        auto node_to_root = T.rsearch(nid_to_consider, true); 
        for (auto curr: node_to_root) {
            for (auto m: curr->mutations) {
                pruned_sample.add_mutation(m);
            }
        }
        
        auto node_to_consider = T.get_node(nid_to_consider);
        int start_range_high = node_to_consider->mutations[0].position;
        int end_range_low = node_to_consider->mutations[0].position;
        for (auto m: node_to_consider->mutations) {
            if (m.position < start_range_high) {
                start_range_high = m.position;
            }
            if (m.position > end_range_low) {
                end_range_low = m.position;
            }
        }
        
        size_t total_nodes = bfs.size();

        std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
        std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);

        std::vector<int> node_set_difference(total_nodes);
        
        size_t best_node_num_leaves = 0;
        int best_set_difference = 1e9;

        std::vector<bool> node_has_unique(total_nodes);
        size_t best_j = 0;
        bool best_node_has_unique = false;

        size_t best_distance = 1e9;

        std::vector<size_t> best_j_vec;

        size_t num_best = 1;
        MAT::Node* best_node = T.root;
        best_j_vec.emplace_back(0);

        std::vector<size_t> node_distance(total_nodes);

        size_t num_recombinations = 0;

        tbb::mutex tbb_lock;

        tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k) {

                if ((nid_to_consider == bfs[k]->identifier) || (T.is_ancestor(nid_to_consider, bfs[k]->identifier)) || 
                        (T.is_ancestor(bfs[k]->identifier, nid_to_consider)) || (num_recombinations >= 10) ||
                        (T.get_num_leaves(bfs[k]) < num_descendants)) {
                    continue;
                }
                
                int start_range_low = 0; 
                int end_range_high = 1e9; 

                node_has_unique[k] = false;

                mapper2_input inp;
                inp.T = &T;
                inp.node = bfs[k];
                inp.missing_sample_mutations = &pruned_sample.sample_mutations;
                inp.excess_mutations = &node_excess_mutations[k];
                inp.imputed_mutations = &imputed_mutations[k];
                inp.best_node_num_leaves = &best_node_num_leaves;
                inp.best_set_difference = &best_set_difference;
                inp.best_node = &best_node;
                inp.best_j =  &best_j;
                inp.num_best = &num_best;
                inp.j = k;
                inp.has_unique = &best_node_has_unique;
                        
                inp.set_difference = &node_set_difference[k];

                inp.distance = node_distance[k];
                inp.best_distance = &best_distance;

                inp.best_j_vec = &best_j_vec;
                inp.node_has_unique = &(node_has_unique);

                mapper2_body(inp, true);

                bool is_recomb = false;

                // Is placement as sibling
                if (bfs[k]->is_leaf() || node_has_unique[k]) {
                    std::vector<MAT::Mutation> l2_mut;
                    
                    // Compute l2_mut
                    for (auto m1: node_excess_mutations[k]) {
                        bool found = false;
                        for (auto m2: bfs[k]->mutations) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            if ((m1.position >= start_range_high) && (m1.position <= end_range_low)) {
                                MAT::Mutation m = m1.copy();
                                l2_mut.emplace_back(m);
                            }
                            if ((m1.position < start_range_high) && (m1.position > start_range_low)) {
                                start_range_low = m1.position;
                            }
                            if ((m1.position > end_range_low) && (m1.position < end_range_high)) {
                                end_range_high = m1.position;
                            }
                        }
                    }

                    if (l2_mut.size() == 0) {
                        is_recomb = true;
                    }
                }
                // Else placement as child
                else {
                    std::vector<MAT::Mutation> node_mut;

                    for (auto m1: node_excess_mutations[k]) {
                        bool found = false;
                        for (auto m2: bfs[k]->mutations) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            if ((m1.position >= start_range_high) && (m1.position <= end_range_low)) {
                                MAT::Mutation m = m1.copy();
                                node_mut.emplace_back(m);
                            }
                            if ((m1.position < start_range_high) && (m1.position > start_range_low)) {
                                start_range_low = m1.position;
                            }
                            if ((m1.position > end_range_low) && (m1.position < end_range_high)) {
                                end_range_high = m1.position;
                            }
                        }
                    }
                    
                    if (node_mut.size() == 0) {
                        is_recomb = true;
                    }
                }

                if (is_recomb) {
                    tbb_lock.lock();
                    if (num_recombinations <= 10) {
                        num_recombinations++;
                        std::string end_range_high_str = (end_range_high == 1e9) ? "GENOME_SIZE" : std::to_string(end_range_high);
                        fprintf(stdout, "Node %s is potentially recombinant with node %s.\n"
                                "Breakpoints are between interval [%i, %i] to interval [%i, %s]. Mutation list on the node follows:\n",  
                                nid_to_consider.c_str(), bfs[k]->identifier.c_str(), start_range_low, start_range_high, end_range_low, 
                                end_range_high_str.c_str());

                        size_t tot_mutations = node_to_consider->mutations.size();
                        for (size_t idx = 0; idx < tot_mutations; idx++) {
                            auto m = node_to_consider->mutations[idx];
                            fprintf(stdout, "%s", m.get_string().c_str());
                            if (idx+1 <tot_mutations) {
                                fprintf(stdout, ",");
                            }
                            else {
                                fprintf(stdout, "\n");
                            }
                        }

                        fprintf(stdout, "Descendants of the recombinant node are:\n");
                        for (auto l: T.get_leaves(node_to_consider->identifier)) {
                            fprintf(stdout, "%s\n", l->identifier.c_str());
                        }
                    }
                    tbb_lock.unlock();
                }
            }       
        }, ap); 
        fprintf(stderr, "Done %zu/%zu.\n", ++num_done, nodes_to_consider.size());
    }

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

