#include "usher.hpp"
#include <cstdio>
#include <string>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/task.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
//add mutation m to parent_mutations, which represent the mutation of a node relative to root,
//or update major allele if already present
static void ins_mut(Mutation_Set &parent_mutations,Mutation_Annotated_Tree::Mutation &m) {
    auto temp = parent_mutations.insert(m);
    //other major allele of leaf node is also counted
    if (!temp.second) {
        //already present
        m.set_par_one_hot(temp.first->get_mut_one_hot());
        //mutate back to ref, so no more mutation
        if (m.get_mut_one_hot() == m.get_ref_one_hot()) {
            parent_mutations.erase(temp.first);
        } else {
            // update state
            const_cast<MAT::Mutation&>(*temp.first).set_mut_one_hot(m.get_mut_one_hot());
        }
    } else {
        if (m.get_mut_one_hot()==m.get_ref_one_hot()) {
            parent_mutations.erase(temp.first);
        }
    }
}
//functor for getting state of all leaves
struct fix_root_worker:public tbb::task {
    Mutation_Annotated_Tree::Node *root; //starting node whose subtree need to be processed
    Mutation_Set parent_mutations; //mutation of parent of "root" relative to the root of the entire tree
    fix_root_worker(Mutation_Annotated_Tree::Node *root,
                          const Mutation_Set &parent_mutations)
        : root(root), parent_mutations(parent_mutations){}
    tbb::task* execute() override {
        //add mutation of "root"
        for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
            ins_mut(parent_mutations, m);
        }
        //continuation
        tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
        //spawn a task for each children
        empty->set_ref_count(root->children.size());
        for (auto child : root->children) {
            assert(child->parent==root);
            empty->spawn(*new (empty->allocate_child())fix_root_worker(child, parent_mutations));
        }
        //bypass the scheduler to execute continuation task directly to fix ref count
        // if no child spawned (otherwise it will hang)
        return root->children.empty()?empty:nullptr;
    }
};
void fix_parent(Mutation_Annotated_Tree::Node *root) {
    Mutation_Set mutations;
        tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())
                                       fix_root_worker(root, mutations));
}
