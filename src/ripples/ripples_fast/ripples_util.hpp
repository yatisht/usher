#pragma once
#include "ripples.hpp"
#include "src/matOptimize/usher_graph.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <array>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <random>
#include <time.h>
#include <vector>

#define CHECK_MAPPER

struct idx_hash {
    size_t operator()(const MAT::Node *in) const { return in->dfs_index; }
};

struct idx_eq {
    bool operator()(const MAT::Node *a, const MAT::Node *b) const {
        return a->dfs_index == b->dfs_index;
    }
};

MAT::Node *get_node_cstr(MAT::Tree &tree, char *name);

struct interval_sorter {
    bool operator()(Recomb_Interval &a, Recomb_Interval &b) {
        if (a.start_range_high < b.start_range_high) {
            return true;
        } else if (a.start_range_high == b.start_range_high &&
                   a.end_range_low < b.end_range_low) {
            return true;
        }
        return false;
    }
};

struct Ripple_Result_Pack {
    std::vector<Recomb_Interval> intervals;
    MAT::Node *node_to_consider;
    int orig_parsimony;
};

struct next_node {
    std::vector<MAT::Node *>::const_iterator &iter;
    std::vector<MAT::Node *>::const_iterator end;
    MAT::Node *operator()(tbb::flow_control &fc) const {
        if (iter == end) {
            fc.stop();
        }
        auto to_resturn = *iter;
        iter++;
        return to_resturn;
    }
};

struct Ripple_Pipeline {
    const std::vector<MAT::Node *> &nodes_to_seach;
    const std::vector<bool> &do_parallel;
    const std::vector<int> &index_map;
    uint32_t branch_len;
    int min_range;
    int max_range;
    uint32_t num_threads;
    int parsimony_improvement;
    MAT::Tree &T;
    const std::vector<Mapper_Info> &traversal_track;
    const unsigned short tree_height;
    Ripple_Result_Pack *operator()(MAT::Node *node_to_consider) const;
};

struct Ripple_Finalizer {
    FILE *desc_file;
    FILE *recomb_file;
    size_t &num_done;
    size_t total_size;
    MAT::Tree &T;

    void operator()(Ripple_Result_Pack *) const;
};
