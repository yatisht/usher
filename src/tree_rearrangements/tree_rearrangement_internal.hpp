#ifndef tree_rearrangement_internal
#define tree_rearrangement_internal
#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <condition_variable>
#include <sys/mman.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include "check_samples.hpp"
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>

/* Enumerater Nodes -(Node*)->  Find neighbors and positions to do fitch-sankoff -(Possible_Moves*)-> do fitch-sankoff -(Possible_Moves*)->  estimate profit of moves -(Profitable_Moves*)-> apply moves*/
//======================Message types========================
struct Fitch_Sankoff_Result{
    MAT::Mutation mutation;
    std::pair<size_t, size_t> range;
    //std::shared_ptr<Fitch_Sankoff::States_Type> original_state;
    char LCA_parent_state;
    Fitch_Sankoff::Scores_Type scores;
};
struct Possible_Moves{
   MAT::Node* src;
   std::vector<MAT::Node*> dsts;
};
struct Move_info{
    MAT::Node* dst;
    MAT::Node* LCA;
    std::vector<MAT::Node*> path;
    //Mutation_Annotated_Tree::Mutations_Collection mutations;
    std::vector<Fitch_Sankoff_Result*> FS_results;
};
struct Candidate_Moves{
    MAT::Node* src;
    std::vector<Move_info> moves;
    std::vector<Fitch_Sankoff_Result> container;
};

struct Fitch_Sankoff_Result_Final{
    MAT::Mutation mutation;
    Fitch_Sankoff::Scores_Type scores;
    //std::shared_ptr<Fitch_Sankoff::States_Type> original_state;
    //size_t original_state_offset;
    char LCA_parent_state;
};

struct Profitable_Move{
    int score_change;
    MAT::Node* src;
    MAT::Node* dst;
    MAT::Node* LCA;
    std::vector<MAT::Node*> path;
    std::pair<size_t, size_t> range;
    std::vector<Fitch_Sankoff_Result_Final> states;
};


struct ConfirmedMove{
    std::vector<MAT::Node*> removed;
    std::vector<MAT::Node*> added;
};


//===========================Functors for implementing the pipeline=====================
typedef tbb::flow::multifunction_node<Possible_Moves*,tbb::flow::tuple<Move_info*>> Parsimony_Score_Calcular_t;
struct Neighbors_Finder{
    int radius;
    Neighbors_Finder(int radius):radius(radius){}
    Possible_Moves* operator()(MAT::Node*)const;
};

struct Parsimony_Score_Calculator{
    const Original_State_t& original_states;
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    Candidate_Moves* operator()(Possible_Moves*)const;  
};
struct Profitable_Moves_Enumerator{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    tbb::concurrent_vector<Profitable_Move*>& profitable_moves;  
    void operator() (Candidate_Moves*)const;
};

struct Move_Executor{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    MAT::Tree& tree;
    std::vector<Profitable_Move*>& moves;
    Pending_Moves_t& tree_edits;
    const Original_State_t& ori;
    mutable std::unordered_map<void*,void*> new_parents_map;
    void operator()(tbb::blocked_range<size_t>&)const;
    private:
    MAT::Node* get_parent(MAT::Node*) const;
};
void resolve_conflict(tbb::concurrent_vector<Profitable_Move*>& candidate_moves, std::vector<Profitable_Move*>& non_conflicting_moves, std::vector<MAT::Node*>& deferred_nodes);
void finalize_children(MAT::Node* parent,ConfirmedMove& edits,MAT::Tree* tree,const Original_State_t& checker);
struct Profitable_Moves_Cacher{
    tbb::concurrent_vector<Profitable_Move>& to_monitor;
    tbb::queuing_rw_mutex& swap_lock;
    std::vector<std::pair<int, size_t>> file_offsets;
    char* filename;
    int fd;
    void* mapped_address;
    size_t length;
    std::condition_variable finish_cv;
    std::mutex finish_mutex;
    bool finished;
    std::thread this_thread;
    Profitable_Moves_Cacher(tbb::concurrent_vector<Profitable_Move>& to_monitor,tbb::queuing_rw_mutex& rw_mutex);
    void run();
    void operator()();
    void get_path(std::vector<MAT::Node*>&);
    void get_element(Profitable_Move&);
    void finish();
    ~Profitable_Moves_Cacher();
};

#endif
