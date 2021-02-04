#include "mutation_annotated_tree.hpp"
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
typedef std::unordered_set<Mutation_Annotated_Tree::Mutation, MutationHash,
                           MutationComparator>
    Mutation_Set;
typedef std::unordered_map<std::string, Mutation_Set>
    Sample_Mut_Type;
void check_samples(
    Mutation_Annotated_Tree::Node *root,
    Sample_Mut_Type &samples);