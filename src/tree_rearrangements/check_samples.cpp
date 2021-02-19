#include "check_samples.hpp"
#include "../mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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
                                  Original_State_t &samples) {
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
        if (old_iter->ref_nuc==old_iter->mut_nuc) {
           out.erase(old_iter);
        }
    }
}

void check_samples_worker_with_pending_moves(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,const Pending_Moves_t& pending_moves) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        ins_mut(parent_mutations, m);
    }
    if (root->is_leaf()) {
        auto leaf_iter=pending_moves.find(root);
        if(leaf_iter!=pending_moves.end()){
            std::string& old_identifier=root->identifier;
            root=leaf_iter->second.added.front();
            assert(old_identifier==root->identifier);
            for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
                ins_mut(parent_mutations, m);
            }
        }

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
                        "[ERROR] Extra mutation to\t%c\t%d\t of Sample\t%s at index \t%zu \n",
                        Mutation_Annotated_Tree::get_nuc(m.mut_nuc), m.position,
                        root->identifier.c_str(),root->index);
                } else {
                    iter->second.erase(m_iter);
                }
            }
            for (auto m_left : iter->second) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at index \t %zu ? \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.mut_nuc),
                        m_left.position, root->identifier.c_str(),root->index);
            }
            samples.erase(iter);
        }
    }
    auto iter=pending_moves.find(root);
    if(iter==pending_moves.end()){
     for (auto child : root->children) {
        assert(child->parent=root);
        check_samples_worker_with_pending_moves(child, parent_mutations, samples,pending_moves);
    }
    }else{
    const ConfirmedMove& pending_move_this_node=iter->second;
    const std::vector<MAT::Node*>& removed_nodes=pending_move_this_node.removed;
    for (auto child : root->children) {
        if(std::find(removed_nodes.begin(),removed_nodes.end(),child)==removed_nodes.end()){
        assert(child->parent=root);
        check_samples_worker_with_pending_moves(child, parent_mutations, samples,pending_moves);
        }
    }
    for(auto new_child:pending_move_this_node.added){
        check_samples_worker_with_pending_moves(new_child, parent_mutations, samples, pending_moves);
    }
    }
}

void check_samples_worker(Mutation_Annotated_Tree::Node *root,
                                 Mutation_Set parent_mutations,
                                 Original_State_t &samples,MAT::Tree* tree) {
    for (Mutation_Annotated_Tree::Mutation &m : root->mutations) {
        ins_mut(parent_mutations, m);
    }

    if (root->is_leaf()) {
        auto iter = samples.find(root->identifier);
        if (iter == samples.end()) {
            fprintf(stderr, "[ERROR] Extra Sample %s \n",
                    root->identifier.c_str());
        }

        else {
            for (auto m : parent_mutations) {
                auto m_iter = iter->second.find(m);
                if (m_iter == iter->second.end()) {
                    fprintf(
                        stderr,
                        "[ERROR] Extra mutation to\t%c\%d\t of Sample\t%s at index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m.mut_nuc), m.position,
                        root->identifier.c_str(),root->index);
                } else {
                    if (m.mut_nuc!=m_iter->mut_nuc) {
                        fprintf(stderr, "Mut Nuc Mismatch at \t %d of sample \t %s at index \t %zu: original \t %c , altered :\t %c \n",m.position,root->identifier.c_str(),root->index,Mutation_Annotated_Tree::get_nuc(m_iter->mut_nuc),MAT::get_nuc(m.mut_nuc));
                    }
                    iter->second.erase(m_iter);
                }
            }

            for (auto m_left : iter->second) {
                fprintf(stderr,
                        "[ERROR] Lost mutation to\t%c\t%d\t of Sample\t%s at index %zu \n",
                        Mutation_Annotated_Tree::get_nuc(m_left.mut_nuc),
                        m_left.position, root->identifier.c_str(),root->index);
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