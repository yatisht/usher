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

void populate_mutated_pos(const Original_State_t& origin_state){
    //global lock for pos_mutexes and mutated_positions
    tbb::queuing_rw_mutex insert_lock;
    //mutex for each mutation_set at each position
    std::unordered_map<int,std::mutex*> pos_mutexes;
    tbb::parallel_for_each(origin_state.begin(),origin_state.end(),[&insert_lock,&pos_mutexes](const std::pair<std::string, Mutation_Set>& sample_mutations){
        for (const MAT::Mutation &m : sample_mutations.second) {
            //reader lock to find this position if it is already inserted
            tbb::queuing_rw_mutex::scoped_lock pos_lock(insert_lock,false);
            auto iter=mutated_positions.find(m);
            std::unordered_map<std::string, nuc_one_hot>* samples;
            std::mutex* sample_mutex;
            if (iter==mutated_positions.end()) {
                //not found, need to insert, with writer lock
                samples=new std::unordered_map<std::string, nuc_one_hot>;
                pos_lock.upgrade_to_writer();
                auto emplace_result=mutated_positions.emplace(m,samples);
                if (!emplace_result.second) {
                    delete samples;
                    samples=emplace_result.first->second;
                }
                //also insert the mutex for that position
                sample_mutex=new std::mutex();
                auto samp_mutex_iter=pos_mutexes.emplace(m.get_position(),sample_mutex);
                if (!samp_mutex_iter.second) {
                    delete sample_mutex;
                }
                sample_mutex=samp_mutex_iter.first->second;
                pos_lock.release();
            }else {
                samples=iter->second;
                sample_mutex=pos_mutexes[m.get_position()];
                pos_lock.release();
            }
            //add sample to mutation mapping
            sample_mutex->lock();
            samples->emplace(sample_mutations.first,
                                 m.get_all_major_allele());
            sample_mutex->unlock();
        }
    });
    //clean up all mutexes
    for(auto m:pos_mutexes){
        delete m.second;
    }
}
//Use Full fitch sankoff to reassign state from scratch
static void reassign_states(MAT::Tree& t, Original_State_t& origin_states){
    auto bfs_ordered_nodes = t.breadth_first_expansion();

    check_samples(t.root, origin_states, &t);
    populate_mutated_pos(origin_states);
    std::unordered_set<std::string> ignored;
    std::unordered_set<std::string> ignored2;
    clean_up_internal_nodes(t.root,t,ignored,ignored2);
    if (t.root->children.size()>1) {
        add_root(&t);
    }
    bfs_ordered_nodes = t.breadth_first_expansion();
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    //get mutation vector
    tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&bfs_ordered_nodes, &output](
            const std::pair<MAT::Mutation,
                            std::unordered_map<std::string, nuc_one_hot> *>
                &pos) {
            std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
            Fitch_Sankoff_Whole_Tree(bfs_ordered_nodes, pos.first, *mutated,
                                     output);
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
}
//load from usher compatible pb
Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states){
    Mutation_Annotated_Tree::Tree t =
        Mutation_Annotated_Tree::load_mutation_annotated_tree(path);
    reassign_states(t, origin_states);
    fprintf(stderr, "original score:%zu\n", t.get_parsimony_score());
    return t;
}