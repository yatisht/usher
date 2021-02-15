#ifndef check_sample
#define check_sample
#include "../mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <utility>
struct ConfirmedMove;
struct MutationComparator {
    bool operator()(const Mutation_Annotated_Tree::Mutation &first,
                    const Mutation_Annotated_Tree::Mutation &second) const {
        return (first.position == second.position) &&
               (first.chrom == second.chrom);
    }
};
struct MutationHash {
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
typedef std::unordered_set<Mutation_Annotated_Tree::Mutation, MutationHash, MutationComparator> Mutation_Set;
typedef std::unordered_map<std::string, Mutation_Set>
    Sample_Mut_Type;
void check_samples(
    Mutation_Annotated_Tree::Node *root,
    Sample_Mut_Type &samples);
void get_mutation_set(Mutation_Annotated_Tree::Node* node, Mutation_Set& out);
void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Sample_Mut_Type &samples);
void check_samples_worker_with_pending_moves(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Sample_Mut_Type &samples,const Pending_Moves_t& pending_moves);
#endif