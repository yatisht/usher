#ifndef USHER_GRAPH
#define USHER_GRAPH
//#include "tree.hpp"
#include "mutation_annotated_tree.hpp"
#include <set>
#include <cassert>
#include <unordered_set>
#include <sys/time.h>
#include <tbb/mutex.h>
#include "Instrumentor.h"

#if SAVE_PROFILE == 1
#  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__);
#else
#  define TIMEIT()
#endif

//extern std::mutex data_lock;

namespace MAT = Mutation_Annotated_Tree;

class Timer {
    private:
        struct timeval m_StartTime, m_EndTime;
    public:
        void Start() {
            gettimeofday(&m_StartTime, NULL);
        }
        long Stop() {
            long useconds, seconds, mseconds;
            gettimeofday(&m_EndTime, NULL);
            useconds = m_EndTime.tv_usec - m_StartTime.tv_usec;
            seconds = m_EndTime.tv_sec - m_StartTime.tv_sec;
            mseconds = ((seconds) * 1000 + useconds/1000.0 + 0.5);
            return mseconds;
        }
        
};

struct mapper_input {
    MAT::Tree* T;
    std::string chrom;
    int8_t ref_nuc;
    int variant_pos;
    std::vector<MAT::Node*>* bfs;
    std::unordered_map<std::string, size_t>* bfs_idx;
    std::vector<std::tuple<size_t, int8_t>> variants;
    std::vector<std::string>* variant_ids;
    
    std::vector<std::string>* missing_samples;
    std::vector<std::vector<MAT::Mutation>>* missing_sample_mutations;
};

struct mapper_body {
    int operator()(mapper_input input);
};

struct mapper2_input {
    std::string missing_sample;
    MAT::Tree* T;
    MAT::Node* node;
    std::vector<MAT::Mutation>* missing_sample_mutations;
    
    int* best_set_difference;
    int* set_difference;
    size_t* best_node_num_leaves;
    size_t j;
    size_t* best_j;
    size_t* num_best;
    MAT::Node** best_node;

    std::vector<bool>* node_has_unique;
    std::vector<size_t>* best_j_vec;
    
    bool* has_unique;

    std::vector<MAT::Mutation>* excess_mutations;
    std::vector<MAT::Mutation>* imputed_mutations;
};

void mapper2_body(mapper2_input& inp, bool compute_parsimony_scores);

#endif