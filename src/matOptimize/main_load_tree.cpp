#include "mutation_annotated_tree.hpp"
#include "Fitch_Sankoff.hpp"
#include <cstdio>
#include <iostream>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/pipeline.h>
#include "tbb/parallel_for_each.h"
#include <tbb/parallel_for.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <unordered_set>
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "apply_move/apply_move.hpp"
#include <tbb/queuing_rw_mutex.h>

namespace MAT=Mutation_Annotated_Tree;
static bool no_valid_mut(MAT::Node* node){
    for(const auto& mut:node->mutations){
        if (mut.is_valid()) {
            return false;
        }
    }
    return true;
}
/**
 * @brief Clean nodes with no valid mutation
 * @param this_node subtree rooted at this_node will be cleaned
 * @param tree
 * @param[out] changed_nodes nodes with their children set changed, need fitch sankoff backward pass
 * @param[out] node_with_inconsistent_state nodes with parent state change, need forward pass
 */
static void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<std::string>& changed_nodes,std::unordered_set<std::string>& node_with_inconsistent_state){

    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::vector<MAT::Node *> this_node_ori_children = this_node->children;
    if (this_node->parent&&(((!this_node->is_leaf())&&no_valid_mut(this_node)))) {
        //Remove this node
        this_node->parent->changed=true;
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        changed_nodes.erase(this_node->identifier);
        node_with_inconsistent_state.erase(this_node->identifier);
        //its parent have changed children set
        changed_nodes.insert(this_node->parent->identifier);
        tree.all_nodes.erase(this_node->identifier);
        //promote all its children, no need to change their mutation vector, as this_node assumed to have no valid mutations
        for (MAT::Node *child : this_node_ori_children) {
            child->changed=true;
            child->have_masked|=this_node->have_masked;
            child->parent = this_node->parent;
            parent_children.push_back(child);
        }
        tree.all_nodes.erase(this_node->identifier);
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        //recurse down
        clean_up_internal_nodes(child, tree,changed_nodes,node_with_inconsistent_state);
    }
}
//For removing nodes with no valid mutations between rounds
void clean_tree(MAT::Tree& t){
    std::unordered_set<std::string> changed_nodes;
    std::unordered_set<std::string> node_with_inconsistent_states;
        clean_up_internal_nodes(t.root, t, changed_nodes,node_with_inconsistent_states);
#ifdef CHECK_STATE_REASSIGN
    MAT::Tree new_tree=reassign_state_full(t);
#endif
    //almost the same as apply move
    std::vector<MAT::Node*> for_reassign;
    for(auto node_str:changed_nodes){
        auto node=t.get_node(node_str);
        for_reassign.push_back(node);
    }
    auto cleaned_count=for_reassign.size();
    fprintf(stderr, "%zu nodes cleaned\n",cleaned_count);
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
    if (cleaned_count) {
        clean_tree(t);
    }
}

void populate_mutated_pos(const Original_State_t& origin_state){
    tbb::parallel_for_each(origin_state.begin(),origin_state.end(),[](const std::pair<std::string, Mutation_Set>& sample_mutations){
        for (const MAT::Mutation &m : sample_mutations.second) {
            //reader lock to find this position if it is already inserted
            auto iter=mutated_positions.find(m);
            tbb::concurrent_unordered_map<std::string, nuc_one_hot>* samples;
            if (iter==mutated_positions.end()) {
                //not found, need to insert, with writer lock
                samples=new tbb::concurrent_unordered_map<std::string, nuc_one_hot>;
                auto emplace_result=mutated_positions.emplace(m,samples);
                if (!emplace_result.second) {
                    delete samples;
                    samples=emplace_result.first->second;
                }
            }else {
                samples=iter->second;
            }
            //add sample to mutation mapping
            samples->emplace(sample_mutations.first,
                                 m.get_all_major_allele());
        }
    });
    //clean up all mutexes
}
//Use Full fitch sankoff to reassign state from scratch
static void reassign_states(MAT::Tree& t, Original_State_t& origin_states,const char* transposed_vcf_path){
    auto bfs_ordered_nodes = t.breadth_first_expansion();

    check_samples(t.root, origin_states, &t);
    {
    std::unordered_set<std::string> ignored;
    std::unordered_set<std::string> ignored2;
    clean_up_internal_nodes(t.root,t,ignored,ignored2);
    //populate_mutated_pos(origin_states);
    }
    bfs_ordered_nodes = t.breadth_first_expansion();
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    if( transposed_vcf_path){
        add_ambuiguous_mutations(transposed_vcf_path,origin_states,t);
    }
    //get mutation vector
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    std::vector<std::pair<MAT::Mutation,tbb::concurrent_vector<std::pair<size_t,char>>>> pos_mutated(MAT::Mutation::refs.size());
    tbb::parallel_for_each(origin_states.begin(),origin_states.end(),[&pos_mutated,t](const std::pair<std::string, Mutation_Set>& sample_mutations){
        for (const MAT::Mutation &m : sample_mutations.second) {
            pos_mutated[m.get_position()].first=m;
            pos_mutated[m.get_position()].second.emplace_back(t.get_node(sample_mutations.first)->bfs_index,m.get_all_major_allele());
        }
    });
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0,pos_mutated.size()),
        [&output,&child_idx_range,&parent_idx,&pos_mutated](const tbb::blocked_range<size_t>& in) {
            for (size_t idx=in.begin(); idx<in.end(); idx++) {
                if (pos_mutated[idx].second.empty()) {
                    continue;
                }
                std::vector<std::pair<long,nuc_one_hot>> mutated_nodes_idx(pos_mutated[idx].second.begin(),pos_mutated[idx].second.end());
                std::sort(mutated_nodes_idx.begin(),mutated_nodes_idx.end(),mutated_t_comparator());
                mutated_nodes_idx.emplace_back(0,0xf);
                Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx, pos_mutated[idx].first, mutated_nodes_idx,
                                     output);
    
            }
        });
    tbb::affinity_partitioner ap;
    //sort and fill
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
    size_t total_mutation_size=0;
    for(const auto node:bfs_ordered_nodes){
	    total_mutation_size+=node->mutations.size();
    }
    fprintf(stderr,"Total mutation size %zu \n", total_mutation_size);

}
//load from usher compatible pb
Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states,const char* transposed_vcf_path){
    fputs("Start loading protobuf\n",stderr);
    Mutation_Annotated_Tree::Tree t =
        Mutation_Annotated_Tree::load_mutation_annotated_tree(path);
    fputs("Finished loading protobuf, start reassigning states\n",stderr);
    reassign_states(t, origin_states,transposed_vcf_path);
    fputs("Finished reassigning states\n",stderr);
    fprintf(stderr, "original score:%zu\n", t.get_parsimony_score());
    return t;
}
