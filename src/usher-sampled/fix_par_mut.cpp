#include "usher.hpp"
#include <cstdio>
#include <string>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<size_t>& changed_nodes_local,std::unordered_set<size_t>& node_with_inconsistent_state);
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
//function for processing nodes recursively
void fix_root_worker_recursive(
    Mutation_Annotated_Tree::Node *root,
    Mutation_Set parent_mutations,
    tbb::task_group &tg) {
    //add mutation of "root"
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        ins_mut(parent_mutations, m);
    }
    
    //spawn a task for each children
    for (auto child : root->children) {
        assert(child->parent==root);
        tg.run([child, parent_mutations, &tg]() {
            fix_root_worker_recursive(child, parent_mutations, tg);
        });
    }
}
void fix_parent(Mutation_Annotated_Tree::Tree &tree) {
    Mutation_Set mutations;
    tbb::task_group tg;
    fix_root_worker_recursive(tree.root, mutations, tg);
    tg.wait();
    auto dfs=tree.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dfs.size()),[&dfs](tbb::blocked_range<size_t> r) {
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            auto & mut=dfs[idx]->mutations.mutations;
            mut.erase(std::remove_if(mut.begin(), mut.end(), [](const MAT::Mutation& mut) {
                return mut.get_mut_one_hot()==mut.get_par_one_hot();
            }),mut.end());
        }
    });
    for (auto node : dfs) {
        node->branch_length=node->mutations.size();
#ifdef NDEBUG
        node->children.reserve(SIZE_MULT*node->children.size());
#endif
    }
    fprintf(stderr, "main dfs size %zu\n",dfs.size());
    {
        std::unordered_set<size_t> ignored1;
        std::unordered_set<size_t> ignored2;
        clean_up_internal_nodes(tree.root,tree,ignored1,ignored2);
    }
}
