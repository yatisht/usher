#ifndef tree_rearrangement_internal
#define tree_rearrangement_internal
#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include "check_samples.hpp"
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/* Enumerater Nodes -(Node*)->  Find neighbors and positions to do fitch-sankoff -(Possible_Moves*)-> do fitch-sankoff -(Possible_Moves*)->  estimate profit of moves -(Profitable_Moves*)-> apply moves*/
//======================Message types========================
struct Fitch_Sankoff_Result{
    MAT::Mutation mutation;
    std::pair<size_t, size_t> range;
    Fitch_Sankoff::States_Type original_state;
    char LCA_parent_state;
    Fitch_Sankoff::States_Type states;
    Fitch_Sankoff::Scores_Type scores;
};

struct Possible_Move{
    MAT::Node* src;
    MAT::Node* dst;
};

struct Move{
    int score_change;
    MAT::Node* src;
    MAT::Node* dst;
    MAT::Node* LCA;
    std::vector<MAT::Node*> path;
    std::vector<Fitch_Sankoff_Result*> states;
};


struct ConfirmedMove{
    std::vector<MAT::Node*> removed;
    std::vector<MAT::Node*> added;
};


//===========================Functors for implementing the pipeline=====================

typedef tbb::flow::multifunction_node<MAT::Node*, tbb::flow::tuple<Possible_Move*> > Neighbors_Finder_t;

struct Neighbors_Finder{
    int radius;
    Neighbors_Finder(int radius):radius(radius){}
    void operator()(MAT::Node*, Neighbors_Finder_t::output_ports_type&)const;
};

struct Profitable_Moves_Enumerator{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    tbb::concurrent_vector<Move*>& profitable_moves;  
    void operator() (Possible_Move*)const;
};

struct Move_Executor{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    MAT::Tree& tree;
    std::vector<Move*>& moves;
    Pending_Moves_t& tree_edits;
    const Sample_Mut_Type& ori;
    mutable std::unordered_map<void*,void*> new_parents_map;
    void operator()(tbb::blocked_range<size_t>&)const;
    private:
    MAT::Node* get_parent(MAT::Node*) const;
};
void resolve_conflict(tbb::concurrent_vector<Move*>& candidate_moves, std::vector<Move*>& non_conflicting_moves, std::vector<MAT::Node*>& deferred_nodes);
void finalize_children(MAT::Node* parent,ConfirmedMove& edits,MAT::Tree* tree,const Sample_Mut_Type& checker);
#endif
