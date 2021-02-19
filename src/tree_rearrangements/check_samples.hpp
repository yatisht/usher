#ifndef check_sample
#define check_sample
#include "../mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <utility>
struct ConfirmedMove;
struct Mutation_Pos_Only_Comparator {
    bool operator()(const Mutation_Annotated_Tree::Mutation &first,
                    const Mutation_Annotated_Tree::Mutation &second) const {
        return (first.position == second.position) &&
               (first.chrom == second.chrom);
    }
};
struct Mutation_Pos_Only_Hash {
    size_t operator()(const Mutation_Annotated_Tree::Mutation &in) const {
        return in.position;
    }
};
struct Node_Idx_Hash{
    size_t operator()(const Mutation_Annotated_Tree::Node* in) const{
        return in->index;
    }
};

struct Node_Idx_Eq{
    bool operator()(const Mutation_Annotated_Tree::Node* first,const Mutation_Annotated_Tree::Node* second)const {
        return first->index==second->index;
    }
};

typedef tbb::concurrent_unordered_map<Mutation_Annotated_Tree::Node*,ConfirmedMove,Node_Idx_Hash,Node_Idx_Eq> Pending_Moves_t;
typedef std::unordered_set<Mutation_Annotated_Tree::Mutation, Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator> Mutation_Set;
typedef std::unordered_map<std::string, Mutation_Set>
    Original_State_t;
void check_samples(
    Mutation_Annotated_Tree::Node *root,
    Original_State_t &samples,Mutation_Annotated_Tree::Tree* tree);
void get_mutation_set(Mutation_Annotated_Tree::Node* node, Mutation_Set& out);
void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,Mutation_Annotated_Tree::Tree* tree=nullptr);
void check_samples_worker_with_pending_moves(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,const Pending_Moves_t& pending_moves);
#endif