#ifndef check_sample
#define check_sample
#include "mutation_annotated_tree.hpp"
#include <cstddef>
#include <cstdio>
#include <tbb/concurrent_unordered_map.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
struct ConfirmedMove;
//Identify 2 mutation as the same if position match, used for overriding mutation as going down to leaf from root
struct Mutation_Pos_Only_Comparator {
    bool operator()(const Mutation_Annotated_Tree::Mutation &first,
                    const Mutation_Annotated_Tree::Mutation &second) const {
        return (first.get_position() == second.get_position());
    }
};
struct Mutation_Pos_Only_Hash {
    size_t operator()(const Mutation_Annotated_Tree::Mutation &in) const {
        return in.get_position();
    }
};
struct Node_Idx_Hash {
    size_t operator()(const Mutation_Annotated_Tree::Node* in) const {
        return in->bfs_index;
    }
};

struct Node_Idx_Eq {
    bool operator()(const Mutation_Annotated_Tree::Node* first,const Mutation_Annotated_Tree::Node* second)const {
        return first->bfs_index==second->bfs_index;
    }
};

typedef std::unordered_set<Mutation_Annotated_Tree::Mutation, Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator> Mutation_Set;
typedef tbb::concurrent_unordered_map<size_t, Mutation_Set>
Original_State_t;
//populate state of all leaves when original state is empty, check if not
void check_samples(
    const Mutation_Annotated_Tree::Node *root,
    Original_State_t &samples,const Mutation_Annotated_Tree::Tree* tree,bool ignore_missed_samples=false);
//void get_mutation_set(Mutation_Annotated_Tree::Node* node, Mutation_Set& out);
/*void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,Mutation_Annotated_Tree::Tree* tree=nullptr);
*/
//void ins_mut(Mutation_Set &parent_mutations,const Mutation_Annotated_Tree::Mutation &m);
#endif