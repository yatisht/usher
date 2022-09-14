#include "mutation_annotated_tree.hpp"
#include "Fitch_Sankoff.hpp"
#include <chrono>
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include "tbb/parallel_for_each.h"
#include <tbb/parallel_for.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <unordered_set>
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "apply_move/apply_move.hpp"
#include <tbb/queuing_rw_mutex.h>

namespace MAT=Mutation_Annotated_Tree;
void populate_mutated_pos(const Original_State_t& origin_state,MAT::Tree& tree) {
    tbb::parallel_for_each(origin_state.begin(),origin_state.end(),[&](const std::pair<size_t, Mutation_Set>& sample_mutations) {
        for (const MAT::Mutation &m : sample_mutations.second) {
            //reader lock to find this position if it is already inserted
            auto iter=mutated_positions.find(m);
            tbb::concurrent_unordered_map<std::string, nuc_one_hot>* samples;
            if (iter==mutated_positions.end()) {
                //not found, need to insert, with writer lock
                samples=new tbb::concurrent_unordered_map<std::string, nuc_one_hot>;
                auto emplace_result=mutated_positions.emplace(m,samples);
                if (!emplace_result.second) {
                    delete samples;
                    samples=emplace_result.first->second;
                }
            } else {
                samples=iter->second;
            }
            //add sample to mutation mapping
            samples->emplace(tree.get_node_name(sample_mutations.first),
                             m.get_all_major_allele());
        }
    });
    //clean up all mutexes
}

//load from usher compatible pb
Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states) {
    fputs("Start loading protobuf\n",stderr);
    Mutation_Annotated_Tree::Tree t;
    if(!Mutation_Annotated_Tree::load_mutation_annotated_tree(path,t)) {
        exit(EXIT_FAILURE);
    }
    fputs("Finished loading protobuf, start reassigning states\n",stderr);
    reassign_states(t, origin_states);
    fputs("Finished reassigning states\n",stderr);
    fprintf(stderr, "original parsimony score:%zu\n", t.get_parsimony_score());
    return t;
}
