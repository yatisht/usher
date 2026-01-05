#pragma once

#include "ripples.hpp"
#include "ripples_util.hpp"
#include "src/matOptimize/usher_graph.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <array>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <random>
#include <time.h>
#include <vector>

struct ripples_runner {
    ripples_runner(MAT::Tree &tree, const std::vector<Missing_Sample> &samples,
                   uint32_t num_threads, uint32_t branch_len, uint32_t num_desc,
                   const std::string &outdir)
        : T(tree), samples_(samples), num_threads_(num_threads),
          branch_len_(branch_len), num_desc_(num_desc), outdir_(outdir) {}

    void operator()() {
        // TODO: Unused constant parameters
        int start_idx = -1;
        int end_idx = -1;
        int max_range = 1e7;
        int min_range = 1e3;
        int parsimony_improvement = 3;

        Timer timer;

        tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads_);
        srand(time(NULL));
        static tbb::affinity_partitioner ap;

        T.uncondense_leaves();
        timer.Start();

        tbb::concurrent_unordered_set<MAT::Node *, idx_hash, idx_eq>
            nodes_to_consider;
        auto dfs = T.depth_first_expansion();

        for (const auto &sample : samples_) {
            for (auto anc : T.rsearch(T.get_node(T.node_name_to_node_idx(sample.name)), true)) {
                if (anc->is_root()) {
                    continue;
                }
                nodes_to_consider.insert(anc);
            }
        }

        std::vector<MAT::Node *> nodes_to_consider_vec;
        for (auto elem : nodes_to_consider) {
            nodes_to_consider_vec.emplace_back(elem);
        }
        std::sort(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end());
        std::shuffle(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end(),
                     std::default_random_engine(0));

        fprintf(stderr, "Found %zu long branches\n", nodes_to_consider.size());
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        fprintf(stderr, "Creating output files.\n");
        boost::filesystem::path path(outdir_);
        if (!boost::filesystem::exists(path)) {
            boost::filesystem::create_directory(path);
        }

        path = boost::filesystem::canonical(outdir_);
        outdir_ = path.generic_string();

        auto desc_filename = outdir_ + "/descendants.tsv";
        fprintf(stderr,
                "Creating file %s to write descendants of recombinant nodes\n",
                desc_filename.c_str());
        FILE *desc_file = fopen(desc_filename.c_str(), "w");
        fprintf(desc_file, "#node_number\tdescendants\n");

        auto recomb_filename = outdir_ + "/recombination.tsv";
        fprintf(stderr, "Creating file %s to write recombination events\n",
                recomb_filename.c_str());
        FILE *recomb_file = fopen(recomb_filename.c_str(), "w");
        fprintf(recomb_file,
                "#recomb_node_number\tbreakpoint-1_interval\tbreakpoint-2_"
                "interval\tdonor_"
                "node_number\tdonor_is_sibling\tdonor_parsimony\tacceptor_node_"
                "number\tacceptor_is_sibling\tacceptor_parsimony\toriginal_"
                "parsimony\tmin_"
                "starting_parsimony\trecomb_parsimony\n");
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        timer.Start();
        std::vector<MAT::Node *> nodes_to_search;
        nodes_to_search.reserve(dfs.size());
        for (auto &node : dfs) {
            if ((node->dfs_end_index - node->dfs_index) >= num_desc_) {
                nodes_to_search.push_back(node);
            }
        }
        std::vector<bool> do_parallel(dfs.size(), false);
        std::vector<Mapper_Info> traversal_track;
        unsigned short tree_height = 0;
        check_parallelizable(T.root, do_parallel,
                             nodes_to_search.size() / num_threads_, num_desc_,
                             tree_height, traversal_track, 0);
        std::vector<int> index_map;
        int node_to_search_idx = 0;
        index_map.reserve(dfs.size());
        for (int dfs_index = 0; dfs_index < (int)dfs.size(); dfs_index++) {
            if (node_to_search_idx != (int)nodes_to_search.size() &&
                (int)nodes_to_search[node_to_search_idx]->dfs_index == dfs_index) {
                index_map.push_back(node_to_search_idx);
                node_to_search_idx++;
            } else {
                index_map.push_back(-node_to_search_idx);
            }
        }

        fprintf(
            stderr,
            "%zu out of %zu nodes have enough descendant to be donor/acceptor",
            nodes_to_search.size(), dfs.size());
        size_t s = 0, e = nodes_to_consider.size();

        if ((start_idx >= 0) && (end_idx >= 0)) {
            s = start_idx;
            if (end_idx <= (int)e) {
                e = end_idx;
            }
        }

        fprintf(stderr,
                "Running placement individually for %zu branches to identify "
                "potential recombination events.\n",
                e - s);

        size_t num_done = 0;
        std::vector<MAT::Node *>::const_iterator cur_iter =
            nodes_to_consider_vec.begin() + s;
        std::vector<MAT::Node *>::const_iterator end =
            nodes_to_consider_vec.begin() + e;
        tbb::parallel_pipeline(
            4, 
            tbb::make_filter<void, MAT::Node *>(
                tbb::filter_mode::serial_in_order, next_node{cur_iter, end}) &
                tbb::make_filter<MAT::Node *, Ripple_Result_Pack *>(
                    tbb::filter_mode::parallel,
                    Ripple_Pipeline{nodes_to_search, do_parallel, index_map,
                                    branch_len_, min_range, max_range,
                                    num_threads_, parsimony_improvement, T,
                                    traversal_track, tree_height}) &
                tbb::make_filter<Ripple_Result_Pack *, void>(
                    tbb::filter_mode::serial_in_order,
                    Ripple_Finalizer{desc_file, recomb_file, num_done,
                                     nodes_to_consider_vec.size(), T}));
        fclose(desc_file);
        fclose(recomb_file);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    MAT::Tree &T;
    const std::vector<Missing_Sample> &samples_;
    uint32_t num_threads_;
    uint32_t branch_len_;
    uint32_t num_desc_;
    std::string outdir_;
};
