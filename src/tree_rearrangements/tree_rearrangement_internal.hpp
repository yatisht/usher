#ifndef tree_rearrangement_internal
#define tree_rearrangement_internal
#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <atomic>
#include <cstdio>
#include <condition_variable>
#include <string>
#include <sys/mman.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include "check_samples.hpp"
#include <tbb/concurrent_vector.h>
#include <tbb/flow_graph.h>
#include <tbb/pipeline.h>
#include <tbb/queuing_rw_mutex.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#define MOVE_FOUND_MASK 1
#define NONE_CONFLICT_MOVE_MASK 2
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
    std::vector<std::shared_ptr<Fitch_Sankoff_Result>> FS_results;
};
struct Candidate_Moves{
    MAT::Node* src;
    size_t state_packed;
    std::vector<Move_info> moves;
    //std::vector<Fitch_Sankoff_Result> container;
};
struct Fitch_Sankoff_Result_Final{
    MAT::Mutation mutation;
    Fitch_Sankoff::Scores_Type scores;
    char LCA_parent_state;
    char LCA_parent_state_original_tree;
    States_To_Set other_move_LCA_parent_states_to_update;
    int optimized_score;
    #ifndef NDEBUG
    int original_topology_score;
    #endif
};

struct Profitable_Move{
    int score_change;
    MAT::Node* src;
    MAT::Node* dst;
    MAT::Node* LCA;
    std::vector<MAT::Node*> path;
    std::pair<size_t, size_t> range;
    all_moves_bare_type other_moves_in_subtree;
    std::vector<Fitch_Sankoff_Result_Final> states;
    #ifdef CHECK_LEAK
    bool destructed;
    ~Profitable_Move(){
        assert(destructed);
    }
    Profitable_Move():destructed(false){}
    #endif
};
    #ifdef CHECK_LEAK
    typedef std::shared_ptr<Profitable_Move> Profitable_Moves_ptr_t;
    #else
    typedef Profitable_Move* Profitable_Moves_ptr_t;
    #endif
struct Profitable_Moves_From_One_Source{
    MAT::Node* src;
    std::vector<Profitable_Moves_ptr_t> profitable_moves;
};
struct ConfirmedMove{
    std::vector<MAT::Node*> removed;
    std::vector<MAT::Node*> added;
    void swap(ConfirmedMove& other){
        removed.swap(other.removed);
        added.swap(other.added);
    }
    bool empty(){
        return removed.empty()&&added.empty();
    }
};


//===========================Functors for implementing the pipeline=====================

struct Neighbors_Finder{
    int radius;
    Neighbors_Finder(int radius):radius(radius){}
    Possible_Moves* operator()(MAT::Node*)const;
};
typedef tbb::flow::function_node<Possible_Moves*,Candidate_Moves*> Parsimony_Score_Calculator_t;
struct Parsimony_Score_Calculator{
    const Original_State_t& original_states;
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    std::atomic_long& states_in_flight;
    std::atomic_long& states_to_calculate;
    Candidate_Moves* operator()(Possible_Moves*)const;  
};
struct Profitable_Moves_Enumerator{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    const Original_State_t& original_states;
    std::atomic_long& states_in_flight;
    Profitable_Moves_From_One_Source* operator() (Candidate_Moves*)const;
};

struct Move_Executor{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    MAT::Tree& tree;
    std::vector<Profitable_Moves_ptr_t>& moves;
    Pending_Moves_t& tree_edits;
    const Original_State_t& ori;
    mutable std::unordered_map<MAT::Node*,MAT::Node*,Node_Idx_Hash,Node_Idx_Eq> new_parents_map;
    void operator()(tbb::blocked_range<size_t>&)const;
    private:
    MAT::Node* get_parent(MAT::Node*) const;
};


void finalize_children(MAT::Node* parent,ConfirmedMove& edits,MAT::Tree* tree,const Original_State_t& checker,std::vector<std::pair<MAT::Node*, MAT::Node*>>& deleted_map);

#endif
