#include "mutation_annotated_tree.hpp"
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tbb/concurrent_unordered_map.h>
#include <thread>
#include <unordered_set>
#include <vector>
#include <condition_variable>
#include "check_samples.hpp"
#include <random>
#pragma once
#define EMPTY_POS UINT16_MAX
#define MOVE_TAG 0
extern int this_rank;
extern int process_count;
extern bool use_bound; 
extern std::chrono::time_point<std::chrono::steady_clock> last_save_time;
extern bool no_write_intermediate;
extern size_t max_queued_moves;
extern uint32_t num_threads;
extern std::chrono::steady_clock::duration save_period;
namespace MAT = Mutation_Annotated_Tree;
extern std::vector<std::string> changed_nodes;
extern tbb::concurrent_unordered_map<MAT::Mutation, tbb::concurrent_unordered_map<std::string, nuc_one_hot>*,Mutation_Pos_Only_Hash,
       Mutation_Pos_Only_Comparator>
       mutated_positions;
extern std::condition_variable progress_bar_cv;
extern bool timed_print_progress;
struct Profitable_Moves {
    int score_change;
    MAT::Node* src;
    MAT::Node* dst;
    MAT::Node* LCA;
    int radius_left;
    std::unordered_set<MAT::Node*> involved_nodes;
    template<typename F>
    void apply_nodes(F f){
        for (const auto node : involved_nodes) {
            f(node);
        }
    }
    void populate_involved_nodes(){
        involved_nodes.reserve(2+2*LCA->level-src->level-dst->level);
        involved_nodes.insert(LCA);
        auto src_ancestor=src;
        auto dst_ancestor=dst;
        while (src_ancestor!=LCA) {
            involved_nodes.insert(src_ancestor);
            src_ancestor=src_ancestor->parent;
        }
        while (dst_ancestor!=LCA) {
            involved_nodes.insert(dst_ancestor);
            dst_ancestor=dst_ancestor->parent;
        }
    }
    MAT::Node* get_src() const{
        return src;
    }
    MAT::Node* get_dst()const{
        return dst;
    }
};
typedef std::shared_ptr<Profitable_Moves> Profitable_Moves_ptr_t;
struct output_t {
    int score_change;
    int radius_left;
    std::vector<Profitable_Moves_ptr_t>* moves;
    output_t():score_change(-1),radius_left(-1) {}
};
int individual_move(Mutation_Annotated_Tree::Node* src,Mutation_Annotated_Tree::Node* dst,Mutation_Annotated_Tree::Node* LCA,output_t& out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    ,MAT::Tree* tree
#endif
                   );
Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states);
void load_vcf_nh_directly( MAT::Tree& t,const std::string& vcf_path,Original_State_t& origin_states);
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t
#ifdef CHECK_STATE_REASSIGN
                 ,
                 const Original_State_t& original_state
#endif
                );
void fix_condensed_nodes(MAT::Tree *tree) ;
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                        std::vector<MAT::Node*> &output,bool is_first,int radius,MAT::Tree &tree) ;
void add_root(MAT::Tree *tree) ;
void VCF_input(const char * name,MAT::Tree& tree);

void optimize_tree_main_thread(std::vector<MAT::Node *> &nodes_to_search,
                                        MAT::Tree &t,int radius,FILE* log,bool allow_drift,int iteration,
                                        std::vector<MAT::Node*>& deferred_nodes_out,bool MPI_involved
#ifndef NDEBUG
                     , Original_State_t& origin_states
#endif
                                       );

void optimize_tree_worker_thread(MAT::Tree &t,int radius,std::vector<MAT::Node *> &nodes_to_search);
void save_final_tree(MAT::Tree &t, Original_State_t& origin_states,const std::string &output_path);
//For removing nodes with no valid mutations between rounds
void clean_tree(MAT::Tree& t);
void populate_mutated_pos(const Original_State_t& origin_state);
void add_ambuiguous_mutations(const char* path,Original_State_t& to_patch,Mutation_Annotated_Tree::Tree& tree);
void recondense_tree(MAT::Tree& t);
void add_ambiguous_mutation(const char *input_path,MAT::Tree& tree);
struct TlRng:public std::mt19937_64 {
    TlRng():std::mt19937_64(std::chrono::steady_clock::now().time_since_epoch().count()*std::hash<std::thread::id>()(std::this_thread::get_id())) {}
};
extern thread_local TlRng rng;

void adjust_all(MAT::Tree &tree) ;
size_t get_memory();
#ifdef CHECK_BOUND
struct counters{
    size_t saved;
    size_t total;
    counters():saved(0),total(0){}
};
#endif