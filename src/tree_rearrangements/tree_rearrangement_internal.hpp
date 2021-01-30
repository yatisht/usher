#ifndef tree_rearrangement_internal
#define tree_rearrangement_internal
#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include <unordered_map>
#include <vector>

/* Enumerater Nodes -(Node*)->  Find neighbors and positions to do fitch-sankoff -(Possible_Moves*)-> do fitch-sankoff -(Possible_Moves*)->  estimate profit of moves -(Profitable_Moves*)-> apply moves*/
//======================Message types========================
struct Fitch_Sankoff_Result{
    int pos;
    int original_tip_score;
    Fitch_Sankoff::Score_Type tip_score;
    std::vector<char> original_state;
    Fitch_Sankoff::States_Type states;
};

struct Dst_Mut{
    MAT::Node* dst;
    MAT::Mutations_Collection mut;
};

struct Fitch_Sankoff_Shared{
    std::pair<size_t, size_t> range;
};

struct Possible_Moves{
    MAT::Node* src;
    std::vector<Dst_Mut> dst;
    Fitch_Sankoff_Shared shared;
    union{
        MAT::Mutations_Collection* to_search;
        std::unordered_map<int, Fitch_Sankoff_Result>* src_tip_fs_result;
    };
};

struct Merge_Discriptor{
    MAT::Node* to_merge_with;
    MAT::Mutations_Collection shared_mutations;
    MAT::Mutations_Collection src_unique_mutations;
    MAT::Mutations_Collection to_merge_width_unique_mutations;
};
struct Move{
    MAT::Node* dst;
    std::vector<std::pair<int,Fitch_Sankoff::States_Type>> states;
    Mutation_Annotated_Tree::Mutations_Collection new_tip_mutations;
    //The child of its new parent that shares mutation, null if none
    std::vector<Merge_Discriptor*> merger;

};

struct Profitable_Moves{
    int score_change;
    MAT::Node* src;
    std::vector<Move> moves;
};

struct ConfirmedMove{
    MAT::Node* src;
    MAT::Node* dst;
    Merge_Discriptor* merge;
};

//======================For synchronization (postpone conflicting moves)===========================
/* This is for the following type of hard conflict (result in loops):
      Some common ancestor
     /      \
    A       B
    |       |
    new B   new A
Once the most profitable move of B is to be an (in)direct child of A, A will not try to move to positions that are (in)direct child of B.
*/
extern tbb::concurrent_hash_map<MAT::Node*, MAT::Node*> potential_crosses;

/* If a child of node has mutation that the node is considering, it will lead to inaccurate parsimony score calculation, but will only be off by one.*/
extern tbb::concurrent_hash_map<MAT::Node*, tbb::concurrent_vector<int>> repeatedly_mutating_loci;

extern tbb::concurrent_vector<MAT::Node*> postponed;


//===========================Functors for implementing the pipeline=====================

//First step: Enumerate movable nodes in leaf to root order (reverse BFS I guess...)
class Movable_Node_Enumerator{
    mutable std::vector<int> this_round;
    mutable std::vector<int>::iterator iter;
    mutable std::vector<int> next_round;
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    Movable_Node_Enumerator(const Movable_Node_Enumerator& other)=delete;
public:
    Movable_Node_Enumerator(std::vector<MAT::Node *>& to_check,std::vector<MAT::Node *>& dfs_ordered_nodes);

    MAT::Node* operator() (tbb::flow_control) const;
};

struct Neighbors_Finder{
    int radius;
    Neighbors_Finder(int radius):radius(radius){}
    Possible_Moves* operator()(MAT::Node*)const;
};

struct Parsimony_Score_Calculator{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    Parsimony_Score_Calculator(std::vector<MAT::Node *>& dfs_ordered_nodes):dfs_ordered_nodes(dfs_ordered_nodes){}
    Possible_Moves* operator()(Possible_Moves*)const;
};

struct Profitable_Moves_Enumerator{
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    Profitable_Moves_Enumerator(std::vector<MAT::Node *>& dfs_ordered_nodes):dfs_ordered_nodes(dfs_ordered_nodes){}
    Profitable_Moves* operator() (Possible_Moves*)const;
};


struct Move_Executor{
    std::vector<ConfirmedMove>& moves;
    std::vector<MAT::Node *>& dfs_ordered_nodes;
    MAT::Tree* tree;
    Move_Executor(std::vector<ConfirmedMove>& moves,std::vector<MAT::Node *>& dfs_ordered_nodes,MAT::Tree* tree):moves(moves),dfs_ordered_nodes(dfs_ordered_nodes),tree(tree){}
    void operator()(Profitable_Moves*)const;
};

#endif