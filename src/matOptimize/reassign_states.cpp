#include "mutation_annotated_tree.hpp"
#include "Fitch_Sankoff.hpp"
#include <chrono>
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_pipeline.h>
#include "tbb/parallel_for_each.h"
#include <tbb/parallel_for.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <unordered_set>
#include "check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include "apply_move/apply_move.hpp"
#include <tbb/queuing_rw_mutex.h>

static bool no_valid_mut(MAT::Node* node) {
    for(const auto& mut:node->mutations) {
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
void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<size_t>& changed_nodes_local,std::unordered_set<size_t>& node_with_inconsistent_state) {

    std::vector<MAT::Node *> this_node_ori_children = this_node->children;
    if (this_node->parent&&(((!this_node->is_leaf())&&no_valid_mut(this_node)))) {
        std::vector<MAT::Node *> &parent_children = this_node->parent->children;
        //Remove this node
        this_node->parent->set_self_changed();
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        changed_nodes_local.erase(this_node->node_id);
        node_with_inconsistent_state.erase(this_node->node_id);
        //its parent have changed children set
        changed_nodes_local.insert(this_node->parent->node_id);
        //tree.all_nodes.erase(this_node->identifier);
        //promote all its children, no need to change their mutation vector, as this_node assumed to have no valid mutations
        for (MAT::Node *child : this_node_ori_children) {
            child->set_self_changed();
            child->have_masked|=this_node->have_masked;
            child->parent = this_node->parent;
            parent_children.push_back(child);
        }
        //tree.all_nodes.erase(this_node->identifier);
        tree.erase_node(this_node->node_id);
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        //recurse down
        clean_up_internal_nodes(child, tree,changed_nodes_local,node_with_inconsistent_state);
    }
}
//For removing nodes with no valid mutations between rounds
void clean_tree(MAT::Tree& t) {
    std::unordered_set<size_t> changed_nodes;
    std::unordered_set<size_t> node_with_inconsistent_states;
    clean_up_internal_nodes(t.root, t, changed_nodes,node_with_inconsistent_states);
#ifdef CHECK_STATE_REASSIGN
    MAT::Tree new_tree=reassign_state_full(t);
#endif
    //almost the same as apply move
    std::vector<MAT::Node*> for_reassign;
    for(auto node_str:changed_nodes) {
        auto node=t.get_node(node_str);
        if(node){
            for_reassign.push_back(node);
        }
    }
    auto cleaned_count=for_reassign.size();
    fprintf(stderr, "%zu nodes cleaned\n",cleaned_count);
    t.depth_first_expansion();
    if(!for_reassign.empty()) {
        std::vector<Altered_Node_t> nodes_with_changed_states_out;
        reassign_backward_pass(for_reassign, nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
                               ,new_tree
#endif
                              );
        for(const auto& node_id:node_with_inconsistent_states) {
            clean_up_src_states(t.get_node(node_id), nodes_with_changed_states_out);
        }
        if (!nodes_with_changed_states_out.empty()) {
            forward_pass(nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
                         ,
                         new_tree
#endif
                        );
        }
    }
#ifdef CHECK_STATE_REASSIGN
    compare_mutation_tree(t, new_tree);
#endif
    if (cleaned_count) {
        clean_tree(t);
    }
}

//Use Full fitch sankoff to reassign state from scratch
void reassign_states(MAT::Tree& t, Original_State_t& origin_states) {
    auto bfs_ordered_nodes = t.breadth_first_expansion();
    auto start_time=std::chrono::steady_clock::now();
    check_samples(t.root, origin_states, &t);
    {
        std::unordered_set<size_t> ignored;
        std::unordered_set<size_t> ignored2;
        clean_up_internal_nodes(t.root,t,ignored,ignored2);
        //populate_mutated_pos(origin_states);
    }
    bfs_ordered_nodes = t.breadth_first_expansion();
    //get mutation vector
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    std::vector<std::pair<MAT::Mutation,tbb::concurrent_vector<std::pair<size_t,char>>>> pos_mutated(MAT::Mutation::refs.size());
    tbb::parallel_for_each(origin_states.begin(),origin_states.end(),[&pos_mutated,&t](const std::pair<size_t, Mutation_Set>& sample_mutations) {
        for (const MAT::Mutation &m : sample_mutations.second) {
            pos_mutated[m.get_position()].first=m;
            pos_mutated[m.get_position()].second.emplace_back(t.get_node(sample_mutations.first)->bfs_index,m.get_all_major_allele());
        }
    });
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    auto prep_end=std::chrono::steady_clock::now();
    auto prep_dur=std::chrono::duration_cast<std::chrono::milliseconds>(prep_end-start_time).count();
    fprintf(stderr, "Preparation took %zu\n",prep_dur);
    FS_result_per_thread_t FS_result;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0,pos_mutated.size()),
    [&FS_result,&child_idx_range,&parent_idx,&pos_mutated](const tbb::blocked_range<size_t>& in) {
        auto& this_result=FS_result.local();
        this_result.init(child_idx_range.size());
        for (size_t idx=in.begin(); idx<in.end(); idx++) {
            if (pos_mutated[idx].second.empty()) {
                continue;
            }
            mutated_t mutated_nodes_idx(pos_mutated[idx].second.begin(),pos_mutated[idx].second.end());
            std::sort(mutated_nodes_idx.begin(),mutated_nodes_idx.end(),mutated_t_comparator());
            mutated_nodes_idx.emplace_back(0,0xf);
            Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx, pos_mutated[idx].first, mutated_nodes_idx,
                                     this_result);

        }
    });
    auto FS_end=std::chrono::steady_clock::now();
    auto FS_dur=std::chrono::duration_cast<std::chrono::milliseconds>(FS_end-prep_end).count();
    fprintf(stderr, "FS took %zu, ratio %f\n",FS_dur,(float)FS_dur/(float)prep_dur);
    fill_muts(FS_result, bfs_ordered_nodes);
    size_t total_mutation_size=0;
    for(const auto node:bfs_ordered_nodes) {
        total_mutation_size+=node->mutations.size();
    }
    fprintf(stderr,"Total mutation size %zu \n", total_mutation_size);

}
