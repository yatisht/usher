#include "check_samples.hpp"
#include "mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#undef NDEBUG
#include <cassert>
void ins_mut(Mutation_Set &parent_mutations,const Mutation_Annotated_Tree::Mutation &m,bool is_leaf) {
    auto temp = parent_mutations.insert(m);
    const_cast<MAT::Mutation&>(*temp.first).set_auxillary(0,0,0);
    if (!temp.second) {
        assert(temp.first->get_mut_one_hot()==m.get_par_one_hot());
        const_cast<MAT::Mutation&>(*temp.first).set_mut_one_hot(m.get_mut_one_hot());
        if (m.get_mut_one_hot() == m.get_ref_one_hot()) {
            parent_mutations.erase(temp.first);
        }else {
            const_cast<MAT::Mutation&>(*temp.first).set_auxillary(0,0,0);
        }
    }else {
        assert(m.get_mut_one_hot() != m.get_ref_one_hot());
        assert(m.get_mut_one_hot() != m.get_par_one_hot()||!m.is_valid());
    }
}

static void insert_samples_worker(Mutation_Annotated_Tree::Node *root,
                                  Mutation_Set parent_mutations,
                                  Original_State_t &samples) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()){ins_mut(parent_mutations, m,root->is_leaf());}
    }
    if (root->is_leaf()) {
        samples.insert(std::make_pair(root->identifier, parent_mutations));
    }
    for (auto child : root->children) {
        insert_samples_worker(child, parent_mutations, samples);
    }
}

void get_mutation_set(Mutation_Annotated_Tree::Node* node, Mutation_Set& out){
    while (node) {
        for(Mutation_Annotated_Tree::Mutation& m : node->mutations){
            out.insert(m);
        }
        node=node->parent;
    }
    for (auto iter=out.begin(); iter!=out.end();) {
        auto old_iter = iter;
        iter++;
        if (old_iter->get_ref_one_hot()==old_iter->get_mut_one_hot()) {
           out.erase(old_iter);
        }
    }
}

void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,MAT::Tree* tree) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        if(m.is_valid()){ins_mut(parent_mutations, m,root->is_leaf());};
    }

    if (root->is_leaf()) {
        auto iter = samples.find(root->identifier);
        if (iter == samples.end()) {
            fprintf(stderr, "[ERROR] Extra Sample %s \n",
                    root->identifier.c_str());
            assert(false);
        }

        else {
            for (auto m : parent_mutations) {
                auto m_iter = iter->second.find(m);
                if (m_iter == iter->second.end()) {
                    fprintf(
                        stderr,
                        "[ERROR] Extra mutation to\t%c\%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m.get_mut_one_hot()), m.get_position(),
                        root->identifier.c_str(),root->bfs_index);
            assert(false);
                        
                } else {
                    if ((m.get_mut_one_hot()|m.get_tie_one_hot())!=m_iter->get_mut_one_hot()) {
                        fprintf(stderr, "Mut Nuc Mismatch at \t %d of sample \t %s at bfs_index \t %zu: original \t %c , altered :\t %c \n",m.get_position(),root->identifier.c_str(),root->bfs_index,Mutation_Annotated_Tree::get_nuc(m_iter->get_mut_one_hot()),MAT::get_nuc(m.get_mut_one_hot()));
            assert(false);
                        
                    }
                    iter->second.erase(m_iter);
                }
            }

            for (auto m_left : iter->second) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at bfs_index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.get_mut_one_hot()),
                        m_left.get_position(), root->identifier.c_str(),root->bfs_index);
            assert(false);
                        
            }

            samples.erase(iter);
        }
    }

    for (auto child : root->children) {
        assert((!tree)||tree->get_node(child->identifier)==child);
        assert(child->parent=root);
        check_samples_worker(child, parent_mutations, samples);
    }
}

void check_samples(Mutation_Annotated_Tree::Node *root,
                   Original_State_t &samples,MAT::Tree* tree) {
    Mutation_Set mutations;
    if (samples.empty()) {
        insert_samples_worker(root, mutations, samples);
    } else {
        check_samples_worker(root, mutations, samples,tree);
        for (auto s : samples) {
            fprintf(stderr, "[ERROR] Missing Sample %s \n",
                    s.first.c_str());
        }
    }
}