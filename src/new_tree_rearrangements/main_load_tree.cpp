#include "mutation_annotated_tree.hpp"
#include "Fitch_Sankoff.hpp"
#include <iostream>
#include <tbb/pipeline.h>
#include "tbb/parallel_for_each.h"
#include <tbb/parallel_for.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"

namespace MAT=Mutation_Annotated_Tree;
void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,tbb::concurrent_vector<MAT::Node*>& to_filter){

    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::vector<MAT::Node *> this_node_ori_children = this_node->children;

    if (this_node->children.size()==1&&this_node->parent) {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        tree.all_nodes.erase(this_node->identifier);
        for (MAT::Node *child : this_node_ori_children) {
            child->parent = this_node->parent;
            parent_children.push_back(child);
        }
        for (size_t node_idx=0; node_idx<to_filter.size(); node_idx++) {
            if (to_filter[node_idx]==this_node) {
                to_filter[node_idx]=this_node_ori_children[0];
            }
        }
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        clean_up_internal_nodes(child, tree,to_filter);
    }
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
    tbb::concurrent_vector<MAT::Node*> ignored;
    clean_up_internal_nodes(t.root,t,ignored);
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