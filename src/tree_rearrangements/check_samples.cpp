#include "check_samples.hpp"
#include "../mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <utility>

static void ins_mut(Mutation_Set &parent_mutations,
                    Mutation_Annotated_Tree::Mutation &m) {
    auto temp = parent_mutations.insert(m);
    if (!temp.second) {
        assert(temp.first->mut_nuc==m.par_nuc);
        temp.first->mut_nuc=m.mut_nuc;
        if (m.mut_nuc == m.ref_nuc) {
            parent_mutations.erase(temp.first);
        }
    }
}

static void insert_samples_worker(Mutation_Annotated_Tree::Node *root,
                                  Mutation_Set parent_mutations,
                                  Sample_Mut_Type &samples) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        ins_mut(parent_mutations, m);
    }
    if (root->is_leaf()) {
        samples.insert(std::make_pair(root->identifier, parent_mutations));
    }
    for (auto child : root->children) {
        insert_samples_worker(child, parent_mutations, samples);
    }
}

static void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Sample_Mut_Type &samples) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        ins_mut(parent_mutations, m);
    }
    if (root->is_leaf()) {
        auto iter = samples.find(root->identifier);
        if (iter == samples.end()) {
            fprintf(stderr, "[ERROR] Extra Sample %s ? \n",
                    root->identifier.c_str());
        } else {
            for (auto m : parent_mutations) {
                auto m_iter = iter->second.find(m);
                if (m_iter == iter->second.end()) {
                    fprintf(
                        stderr,
                        "[ERROR] Extra mutation to\t%c\%d\t of Sample\t%s at index %zu? \n",
                        Mutation_Annotated_Tree::get_nuc(m.mut_nuc), m.position,
                        root->identifier.c_str(),root->index);
                } else {
                    iter->second.erase(m_iter);
                }
            }
            for (auto m_left : iter->second) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at index %zu ? \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.mut_nuc),
                        m_left.position, root->identifier.c_str(),root->index);
            }
            samples.erase(iter);
        }
    }
    for (auto child : root->children) {
        assert(child->parent=root);
        check_samples_worker(child, parent_mutations, samples);
    }
}
void check_samples(Mutation_Annotated_Tree::Node *root,
                   Sample_Mut_Type &samples) {
    Mutation_Set mutations;
    if (samples.empty()) {
        insert_samples_worker(root, mutations, samples);
    } else {
        check_samples_worker(root, mutations, samples);
        for (auto s : samples) {
            fprintf(stderr, "[ERROR] Missing Sample %s ? \n",
                    s.first.c_str());
        }
    }
}