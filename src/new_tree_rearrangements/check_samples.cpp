#include "check_samples.hpp"
#include "mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <cstdio>
#include <string>
#include <tbb/task.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#undef NDEBUG
#include <cassert>
//add mutation m to parent_mutations, which represent the mutation of a node relative to root, 
//or update major allele if already present
void ins_mut(Mutation_Set &parent_mutations,const Mutation_Annotated_Tree::Mutation &m,bool is_leaf) {
    auto temp = parent_mutations.insert(m);
    //other major allele of leaf node is also counted
    if(!is_leaf)
{    const_cast<MAT::Mutation&>(*temp.first).set_auxillary(temp.first->get_mut_one_hot(),0);
}    if (!temp.second) {
        //already present
        assert(temp.first->get_mut_one_hot()==m.get_par_one_hot());
        //mutate back to ref, so no more mutation
        if ((m.get_mut_one_hot() == m.get_ref_one_hot()&&(!is_leaf))||(m.get_all_major_allele()==m.get_ref_one_hot())) {
            parent_mutations.erase(temp.first);
        }else {
        // update state
        const_cast<MAT::Mutation&>(*temp.first).set_mut_one_hot(m.get_mut_one_hot());
        const_cast<MAT::Mutation&>(*temp.first).set_auxillary(is_leaf?m.get_all_major_allele():m.get_mut_one_hot(),0);
        }
    }else {
        assert(m.get_par_one_hot() == m.get_ref_one_hot());
        assert(m.get_mut_one_hot() != m.get_par_one_hot()||!m.is_valid());
    }
}
//functor for getting state of all leaves
struct insert_samples_worker:public tbb::task {
    Mutation_Annotated_Tree::Node *root; //starting node whose subtree need to be processed
    Mutation_Set parent_mutations; //mutation of parent of "root" relative to the root of the entire tree
    Original_State_t &samples; //output
    insert_samples_worker(Mutation_Annotated_Tree::Node *root,
                          const Mutation_Set &parent_mutations,
                          Original_State_t &samples)
        : root(root), parent_mutations(parent_mutations),
          samples(samples) {}
    tbb::task* execute() override{
        //add mutation of "root"
        for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()||root->is_leaf()){ins_mut(parent_mutations, m,root->is_leaf());}
    }
    //output
    if (root->is_leaf()) {
        samples.insert(std::make_pair(root->identifier, parent_mutations));
    }
    //continuation
    tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
    //spawn a task for each children
    empty->set_ref_count(root->children.size());
    for (auto child : root->children) {
        assert(child->parent==root);
            empty->spawn(*new (empty->allocate_child())insert_samples_worker(child, parent_mutations, samples));
    }
    //bypass the scheduler to execute continuation task directly to fix ref count
    // if no child spawned (otherwise it will hang)
    return root->children.empty()?empty:nullptr;
    }
};

//functor for checking state of all leaves
struct check_samples_worker:public tbb::task{
    const Mutation_Annotated_Tree::Node *root;
                                 Mutation_Set parent_mutations;
                                 const Original_State_t &samples;
check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 const Mutation_Set& parent_mutations,
                                 const Original_State_t &samples):root(root),parent_mutations(parent_mutations),samples(samples){}
tbb::task* execute() override{
    tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
    empty->set_ref_count(root->children.size());
    for (const Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()||root->is_leaf()){ins_mut(parent_mutations, m,root->is_leaf());};
    }

    if (root->is_leaf()) {
        auto iter = samples.find(root->identifier);
        if (iter == samples.end()) {
            fprintf(stderr, "[ERROR] Extra Sample %s \n",
                    root->identifier.c_str());
            assert(false);
        }else {
            Mutation_Set to_check(iter->second);
            for (auto m : parent_mutations) {
                auto m_iter = to_check.find(m);
                if (m_iter == to_check.end()) {
                    fprintf(
                        stderr,
                        "[ERROR] Extra mutation to\t%c\%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m.get_all_major_allele()), m.get_position(),
                        root->identifier.c_str(),root->bfs_index);
            assert(false);
                        
                } else {
                    if ((m.get_all_major_allele())!=m_iter->get_all_major_allele()) {
                        fprintf(stderr, "Mut Nuc Mismatch at \t %d of sample \t %s at bfs_index \t %zu: original \t %c , altered :\t %c \n",m.get_position(),root->identifier.c_str(),root->bfs_index,Mutation_Annotated_Tree::get_nuc(m_iter->get_all_major_allele()),MAT::get_nuc(m.get_all_major_allele()));
            assert(false);
                        
                    }
                    to_check.erase(m_iter);
                }
            }

            for (auto m_left : to_check) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.get_all_major_allele()),
                        m_left.get_position(), root->identifier.c_str(),root->bfs_index);
            assert(false);
                        
            }

            //samples.erase(iter);
        }
        return empty;
    }
    for (auto child : root->children) {
        assert(child->parent==root);
        empty->spawn(*new (empty->allocate_child())check_samples_worker(child, parent_mutations, samples));
    }
    return nullptr;
}
};
//top level
void check_samples(Mutation_Annotated_Tree::Node *root,
                   Original_State_t &samples,MAT::Tree* tree) {
    Mutation_Set mutations;
    if (samples.empty()) {
        tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) 
        insert_samples_worker(root, mutations, samples));
    } else {
        tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())check_samples_worker(root, mutations, samples));
        /*for (auto s : samples) {
            fprintf(stderr, "[ERROR] Missing Sample %s \n",
                    s.first.c_str());
        }*/
        fputs("checked\n", stderr);
    }
}