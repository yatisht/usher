#include "mutation_annotated_tree.hpp"
#include "Fitch_Sankoff.hpp"
#include <cstdio>
#include <iostream>
#include <string>
#include <tbb/pipeline.h>
#include "tbb/parallel_for_each.h"
#include <tbb/parallel_for.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <unordered_set>
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "apply_move/apply_move.hpp"

namespace MAT=Mutation_Annotated_Tree;
static bool no_valid_mut(MAT::Node* node){
    for(const auto& mut:node->mutations){
        if (mut.is_valid()) {
            return false;
        }
    }
    return true;
}
static void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<std::string>& changed_nodes,std::unordered_set<std::string>& node_with_inconsistent_state){

    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::vector<MAT::Node *> this_node_ori_children = this_node->children;

    if (this_node->parent&&(this_node->children.size()==1||((!this_node->is_leaf())&&no_valid_mut(this_node)))) {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        changed_nodes.erase(this_node->identifier);
        node_with_inconsistent_state.erase(this_node->identifier);
        changed_nodes.insert(this_node->parent->identifier);
        tree.all_nodes.erase(this_node->identifier);
        for (MAT::Node *child : this_node_ori_children) {
            child->parent = this_node->parent;
            /*if (this_node_ori_children.size() == 1) {
                if (merge_mutation_single_child(child, this_node->mutations)) {
                    node_with_inconsistent_state.insert(child->identifier);
                }
                if (child->children.size() <= 1) {
                    auto &child_mut = child->mutations;
                    for (auto &mut : child_mut) {
                        mut.set_boundary_one_hot(0xf &
                                                 (~mut.get_all_major_allele()));
                    }
                    child_mut.mutations.erase(
                        std::remove_if(child_mut.begin(), child_mut.end(),
                                       [](const MAT::Mutation &mut) {
                                           return mut.get_all_major_allele() ==
                                                  mut.get_par_one_hot();
                                       }),
                        child_mut.end());
                }
            }*/
            parent_children.push_back(child);
        }
        tree.all_nodes.erase(this_node->identifier);
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        clean_up_internal_nodes(child, tree,changed_nodes,node_with_inconsistent_state);
    }
}
/*
char get_sample_state(char* sample,int pos){
    return (*mutated_positions[pos])[sample].get_nuc_no_check();
}*/
void clean_tree(MAT::Tree& t,std::unordered_set<std::string>& changed_nodes){
    std::unordered_set<std::string> node_with_inconsistent_states;
        clean_up_internal_nodes(t.root, t, changed_nodes,node_with_inconsistent_states);
#ifdef CHECK_STATE_REASSIGN
    MAT::Tree new_tree=reassign_state_full(t);
#endif
    std::vector<MAT::Node*> for_reassign;
    for(auto node_str:changed_nodes){
        auto node=t.get_node(node_str);
        for_reassign.push_back(node);
    }
    fprintf(stderr, "%zu nodes cleaned\n",for_reassign.size());
    t.depth_first_expansion();
    if(!for_reassign.empty()){
        std::vector<Altered_Node_t> nodes_with_changed_states_out;
        reassign_backward_pass(for_reassign, nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
        ,new_tree
#endif
        );
        for(const auto& node_id:node_with_inconsistent_states){
            clean_up_src_states(t.get_node(node_id), nodes_with_changed_states_out);
        }
        if (!nodes_with_changed_states_out.empty()) {
            forward_pass(nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
                     ,
                     new_tree
#endif
        );
    }}
#ifdef CHECK_STATE_REASSIGN
    compare_mutation_tree(t, new_tree);
#endif
}

static void reassign_states(MAT::Tree& t, Original_State_t& origin_states){
    auto bfs_ordered_nodes = t.breadth_first_expansion();

    for (MAT::Node *node : bfs_ordered_nodes) {
        for (const MAT::Mutation &m : node->mutations) {
            mutated_positions.emplace(
                m, new std::unordered_map<std::string, nuc_one_hot>);
        }
        node->tree = &t;
    }

    check_samples(t.root, origin_states, &t);
    std::unordered_set<std::string> ignored;
    std::unordered_set<std::string> ignored2;
    clean_up_internal_nodes(t.root,t,ignored,ignored2);
    if (t.root->children.size()>1) {
        add_root(&t);
    }
    bfs_ordered_nodes = t.breadth_first_expansion();
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&origin_states, &bfs_ordered_nodes, &output](
            const std::pair<MAT::Mutation,
                            std::unordered_map<std::string, nuc_one_hot> *>
                &pos) {
            std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
            for (auto &sample : origin_states) {
                auto iter = sample.second.find(pos.first);
                if (iter != sample.second.end()) {
                    mutated->emplace(sample.first, iter->get_all_major_allele());
                }
            }
            Fitch_Sankoff_Whole_Tree(bfs_ordered_nodes, pos.first, *mutated,
                                     output);
        });
    tbb::affinity_partitioner ap;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto &to_refill = output[i];
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size());
            }
        },
        ap);
}

Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states){
    Mutation_Annotated_Tree::Tree t =
        Mutation_Annotated_Tree::load_mutation_annotated_tree(path);
    reassign_states(t, origin_states);
    fprintf(stderr, "original score:%zu\n", t.get_parsimony_score());
    return t;
}
MAT::Tree load_vcf_nh_directly(const std::string& nh_path,const std::string& vcf_path,Original_State_t& origin_states){
    MAT::Tree ret=Mutation_Annotated_Tree::create_tree_from_newick(nh_path);
    VCF_input(vcf_path.c_str(),ret);
    ret.condense_leaves();
    fprintf(stderr, "%zu condensed_nodes",ret.condensed_nodes.size());
    ret.collapse_tree();
    reassign_states(ret, origin_states);
    return ret;
}