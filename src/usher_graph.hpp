#include "tree.hpp"
#include <set>
#include <cassert>
#include <unordered_set>
#include <mutex>
#include <sys/time.h>
#include <tbb/mutex.h>
#include "Instrumentor.h"

#if SAVE_PROFILE == 1
#  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__);
#else
#  define TIMEIT()
#endif

//extern std::mutex data_lock;

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

struct mutation {
    std::string chrom;
    int position;
    int8_t ref_nuc;
    int8_t par_nuc;
    bool is_missing;
    std::vector<int8_t> mut_nuc;
};

static inline bool compare_by_position (const mutation & a, const mutation& b) { 
    return a.position < b.position; 
}

struct mapper_input {
    Tree* T;
    std::string chrom;
    int8_t ref_nuc;
    int variant_pos;
    std::vector<Node*>* bfs;
    std::unordered_map<std::string, size_t>* bfs_idx;
    std::vector<std::tuple<size_t, std::vector<int8_t>>> variants;
    std::vector<std::string>* variant_ids;
    
    std::vector<std::string>* missing_samples;
    std::unordered_map<Node*, std::vector<mutation>>* node_mutations;
    std::vector<std::vector<mutation>>* missing_sample_mutations;
};

struct mapper_body {
    int operator()(mapper_input input);
};

struct mapper2_input {
    std::string missing_sample;
    Tree* T;
    Node* node;
    std::unordered_map<Node*, std::vector<mutation>>* node_mutations;
    std::vector<mutation>* missing_sample_mutations;
    
    int* best_set_difference;
    int* set_difference;
    size_t* best_node_num_leaves;
    size_t j;
    size_t* best_j;
    size_t* num_best;
    Node** best_node;

#if DEBUG == 1
    std::vector<bool>* node_has_unique;
    std::vector<size_t>* best_j_vec;
#endif
    
    bool* has_unique;

    std::vector<mutation>* excess_mutations;
    std::vector<mutation>* imputed_mutations;
};

void mapper2_body(mapper2_input& inp, bool compute_parsimony_scores);
