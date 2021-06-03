#include "tree_rearrangement_internal.hpp"
#include <cstdio>
#include <string>
#include "apply_move/apply_move.hpp"
#include <tbb/task.h>
struct clean_up_internal_nodes_single_child_or_no_mutation:public tbb::task{
MAT::Node* this_node;
tbb::concurrent_vector<std::string>& changed_nodes;
tbb::concurrent_vector<std::string>& node_with_inconsistent_state;
tbb::concurrent_vector<std::string>& removed_nodes;
std::vector<std::mutex>& node_mutexes;
clean_up_internal_nodes_single_child_or_no_mutation(
    MAT::Node* this_node,
tbb::concurrent_vector<std::string>& changed_nodes,
tbb::concurrent_vector<std::string>& node_with_inconsistent_state,
tbb::concurrent_vector<std::string>& removed_nodes,
std::vector<std::mutex>& node_mutexes
):this_node(this_node),changed_nodes(changed_nodes),node_with_inconsistent_state(node_with_inconsistent_state),removed_nodes(removed_nodes),node_mutexes(node_mutexes){}
tbb::task* execute() override{

    std::vector<MAT::Node *> this_node_ori_children = this_node->children;
    if (this_node->parent) {
    
    
    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::mutex& parent_mutex=node_mutexes[this_node->parent->dfs_index];

    if ((this_node->children.size()==1||((!this_node->is_leaf())&&this_node->no_valid_mutation()))) {
        parent_mutex.lock();
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        parent_mutex.unlock();
        changed_nodes.push_back(this_node->parent->identifier);
        for (MAT::Node *child : this_node_ori_children) {
            child->parent = this_node->parent;
            if (this_node_ori_children.size() == 1) {
                if (merge_mutation_single_child(child, this_node->mutations)) {
                    node_with_inconsistent_state.push_back(child->identifier);
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
            }
            parent_mutex.lock();
            parent_children.push_back(child);
            parent_mutex.unlock();
        }
        removed_nodes.push_back(this_node->identifier);
        delete this_node;
    }
    }
    tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
    bool spawned=false;
    this_node_ori_children.erase(std::remove_if(this_node_ori_children.begin(), this_node_ori_children.end(), [](MAT::Node* node){return node->is_leaf();}),this_node_ori_children.end());
    empty->set_ref_count(this_node_ori_children.size());
    for (MAT::Node *child : this_node_ori_children) {
        empty->spawn(*new (empty->allocate_child())clean_up_internal_nodes_single_child_or_no_mutation(child,changed_nodes,node_with_inconsistent_state,removed_nodes,node_mutexes));
        spawned=true;
    }
    return spawned?nullptr:empty;
}
};

static void clean_up_internal_nodes_single_child_or_no_mutation(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<std::string>& changed_nodes,std::unordered_set<std::string>& node_with_inconsistent_state){
    tbb::concurrent_vector<std::string> changed_nodes_rep;
    tbb::concurrent_vector<std::string> node_with_inconsistent_state_rep;
    tbb::concurrent_vector<std::string> removed_nodes;
    auto dfs=tree.depth_first_expansion();
    std::vector<std::mutex> node_mutexes(dfs.size());
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())struct clean_up_internal_nodes_single_child_or_no_mutation(this_node,changed_nodes_rep,node_with_inconsistent_state_rep,removed_nodes,node_mutexes));
    changed_nodes.insert(changed_nodes_rep.begin(),changed_nodes_rep.end());
    node_with_inconsistent_state.insert(node_with_inconsistent_state_rep.begin(),node_with_inconsistent_state_rep.end());
    for( const auto& node_id:removed_nodes){
        changed_nodes.erase(node_id);
        node_with_inconsistent_state.erase(node_id);
        tree.all_nodes.erase(node_id);
    }
}
char get_major_allele(MAT::Node* node, int position){
    auto iter=node->mutations.find(position);
    if (iter==node->mutations.end()) {
        return 0;
    }else {
        return iter->get_all_major_allele();
    }
}
static void clean_tree_load(MAT::Tree& t,std::unordered_set<std::string>& changed_nodes,Original_State_t& ori_state){
    std::unordered_set<std::string> node_with_inconsistent_states;
        clean_up_internal_nodes_single_child_or_no_mutation(t.root, t, changed_nodes,node_with_inconsistent_states);
        //check_samples(t.root, ori_state, &t);
#ifdef CHECK_STATE_REASSIGN
    MAT::Tree new_tree=reassign_state_full(t);
    new_tree.save_detailed_mutations("Reassigned_Ref.pb");
    //MAT::Tree new_tree;
    //new_tree.load_detatiled_mutations("Reassigned_Ref.pb");
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
char get_state(char* sample,int position){
    auto iter=mutated_positions.find(position);
    if (iter==mutated_positions.end()) {
        return MAT::Mutation::refs[position].get_nuc_no_check();
    }else {
        auto sample_iter=iter->second->find(sample);
        if (iter->second->end()==sample_iter) {
        return MAT::Mutation::refs[position].get_nuc_no_check();
        }else {
            return sample_iter->second.get_nuc_no_check();
        }
    }
}
MAT::Tree load_vcf_nh_directly(const std::string& nh_path,const std::string& vcf_path,Original_State_t& origin_states){
    MAT::Tree t=Mutation_Annotated_Tree::create_tree_from_newick(nh_path);
    VCF_input(vcf_path.c_str(),t);
    puts("Finished loading from VCF\n");
    std::unordered_set<std::string> changed_nodes;
    clean_tree_load(t, changed_nodes,origin_states);
    t.condense_leaves();
    check_samples(t.root, origin_states, &t);
    //populate_mutated_pos(origin_states);
    changed_nodes.clear();
    printf("%zu condensed nodes\n",t.condensed_nodes.size());
    for(const auto &condensed:t.condensed_nodes){
        changed_nodes.insert(t.get_node(condensed.first)->parent->identifier);
    }
    clean_tree_load(t, changed_nodes,origin_states);
    return t;
}