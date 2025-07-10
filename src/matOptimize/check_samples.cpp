#include "check_samples.hpp"
#include "mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <string>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/task_group.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
//add mutation m to parent_mutations, which represent the mutation of a node relative to root,
//or update major allele if already present
void ins_mut(Mutation_Set &parent_mutations,const Mutation_Annotated_Tree::Mutation &m,bool is_leaf,const MAT::Node* root) {
    auto temp = parent_mutations.insert(m);
    //other major allele of leaf node is also counted
    if(!is_leaf) {
        const_cast<MAT::Mutation&>(*temp.first).set_auxillary(temp.first->get_mut_one_hot(),0);
    }
    if (!temp.second) {
        //already present
        if(temp.first->get_mut_one_hot()!=m.get_par_one_hot()) {
            fprintf(stderr, "nuc mismatch, at pos %d, node id %zu, expect %d, got %d\n",m.get_position(),root->node_id,temp.first->get_mut_one_hot().get_nuc_no_check(),m.get_par_one_hot().get_nuc_no_check());
            //raise(SIGTRAP);
        }
        //mutate back to ref, so no more mutation
        if ((m.get_mut_one_hot() == m.get_ref_one_hot()&&(!is_leaf))||(m.get_all_major_allele()==m.get_ref_one_hot())) {
            parent_mutations.erase(temp.first);
        } else {
            // update state
            const_cast<MAT::Mutation&>(*temp.first).set_mut_one_hot(m.get_mut_one_hot());
            const_cast<MAT::Mutation&>(*temp.first).set_auxillary(is_leaf?m.get_all_major_allele():m.get_mut_one_hot(),0);
        }
    } else {
        if (m.get_par_one_hot()!=m.get_ref_one_hot()) {
            fprintf(stderr, "nuc mismatch, at pos %d\n",m.get_position());
            //raise(SIGTRAP);
        }
        if (m.get_all_major_allele()==m.get_ref_one_hot()) {
            parent_mutations.erase(temp.first);
        }
    }
}
//function for getting state of all leaves
void insert_samples_worker(const Mutation_Annotated_Tree::Node *root,
                           Mutation_Set parent_mutations,
                           Original_State_t &samples,
                           tbb::task_group &tg) {
    //add mutation of "root"
    for (const Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()||root->is_leaf()) {
            ins_mut(parent_mutations, m,root->is_leaf(),root);
        }
    }
    //output
    if (root->is_leaf()) {
        samples.insert(std::make_pair(root->node_id, parent_mutations));
    }
    //process children in parallel
    for (auto child : root->children) {
        assert(child->parent==root);
        tg.run([=, &samples, &tg]() {
            insert_samples_worker(child, parent_mutations, samples, tg);
        });
    }
}

//function for checking state of all leaves
void check_samples_worker(const Mutation_Annotated_Tree::Node *root,
                          Mutation_Set parent_mutations,
                          const Original_State_t &samples,
                          tbb::concurrent_unordered_set<size_t>& visited_samples,
                          const MAT::Tree& tree,
                          tbb::task_group &tg) {
    for (const Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()||root->is_leaf()) {
            ins_mut(parent_mutations, m,root->is_leaf(),root);
        };
    }

    if (root->is_leaf()) {
        auto iter = samples.find(root->node_id);
        if (iter == samples.end()) {
            fprintf(stderr, "[ERROR] Extra Sample %s \n",tree.get_node_name(root->node_id).c_str());
            raise(SIGTRAP);
        } else {
            Mutation_Set to_check(iter->second);
            for (auto m : parent_mutations) {
                auto m_iter = to_check.find(m);
                if (m_iter == to_check.end()) {
                    fprintf(
                        stderr,
                        "[ERROR] Extra mutation to\t%c\t%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m.get_all_major_allele()), m.get_position(),
                        tree.get_node_name(root->node_id).c_str(),root->bfs_index);
                    raise(SIGTRAP);

                } else {
                    if ((m.get_all_major_allele())!=m_iter->get_all_major_allele()) {
                        fprintf(stderr,
                                "Mut Nuc Mismatch at \t %d, address %lx, of sample \t %s "
                                "at bfs_index \t %zu: original \t %d , "
                                "altered :\t %d \n",
                                m.get_position(),(unsigned long)root,
                                tree.get_node_name(root->node_id).c_str(),
                                root->bfs_index,
                                (int)m_iter->get_all_major_allele(),
                                (int)m.get_all_major_allele());
                        raise(SIGTRAP);
                    }
                    to_check.erase(m_iter);
                }
            }

            for (auto m_left : to_check) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.get_all_major_allele()),
                        m_left.get_position(), tree.get_node_name(root->node_id).c_str(),root->bfs_index);
                raise(SIGTRAP);
            }

            visited_samples.insert(iter->first);
        }
        return;
    }
    for (auto child : root->children) {
        if (child->parent!=root) {
            fprintf(stderr, "%lx\n", (long)child);
            std::raise(SIGTRAP);
        }
        tg.run([=, &samples, &visited_samples, &tree, &tg]() {
            check_samples_worker(child, parent_mutations, samples, visited_samples, tree, tg);
        });
    }
}
//top level
void check_samples(const Mutation_Annotated_Tree::Node *root,
                   Original_State_t &samples,const MAT::Tree* tree,bool ignore_missed_samples) {
    Mutation_Set mutations;
    if (samples.empty()) {
        tbb::task_group tg;
        insert_samples_worker(root, mutations, samples, tg);
        tg.wait();
    } else {
        tbb::concurrent_unordered_set<size_t> visited_sample;
        tbb::task_group tg;
        check_samples_worker(root, mutations, samples, visited_sample, *tree, tg);
        tg.wait();
        if (!ignore_missed_samples) {
            bool have_missed=false;
            for (auto s : samples) {
                if (!visited_sample.count(s.first)) {
                    fprintf(stderr, "[ERROR] Missing Sample %s \n",
                            tree->get_node_name(s.first).c_str());
                    have_missed=true;
                }
            }
            if (have_missed) {
                fprintf(stderr, "have_missing samples\n");
            }
            assert(!have_missed);
        }
        fputs("checked\n", stderr);
    }
}
