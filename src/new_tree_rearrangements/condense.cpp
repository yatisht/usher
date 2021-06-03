#include "mutation_annotated_tree.hpp"
#include "tbb/task.h"
#include <atomic>
static void intersect_allele(const Mutation_Annotated_Tree::Mutations_Collection& in1, const Mutation_Annotated_Tree::Mutations_Collection& in2, Mutation_Annotated_Tree::Mutations_Collection& out){
    if (in1.empty()||in2.empty()) {
        return;
    }
    auto in1_iter=in1.begin();
    auto in1_end=in1.end();
    for(const auto& mut:in2){
        while (in1_iter!=in1_end&&in1_iter->get_position()<mut.get_position()) {
            in1_iter++;
        }
        if (in1_iter!=in1_end&&in1_iter->get_position()==mut.get_position()){
            nuc_one_hot common=in1_iter->get_all_major_allele()&mut.get_all_major_allele();
            if (common!=mut.get_par_one_hot()) {
                out.push_back(mut);
                out.back().set_auxillary(common, 0xf&(~common));
            }
            in1_iter++;
        }
    }
}

static void get_common_mutations(std::vector<Mutation_Annotated_Tree::Node*>::const_iterator begin,std::vector<Mutation_Annotated_Tree::Node*>::const_iterator end,Mutation_Annotated_Tree::Mutations_Collection& out){
    if (end-begin==3) {
        Mutation_Annotated_Tree::Mutations_Collection temp;
        intersect_allele((*begin)->mutations, (*(begin+1))->mutations, temp);
        intersect_allele(temp, (*(begin+2))->mutations, out);
    } else if (end-begin==2) {
        intersect_allele((*begin)->mutations, (*(begin+1))->mutations, out);
    } else {
        Mutation_Annotated_Tree::Mutations_Collection temp1;
        Mutation_Annotated_Tree::Mutations_Collection temp2;
        auto middle=begin+(end-begin)/2;
        get_common_mutations(begin,middle , temp1);
        get_common_mutations(middle,end,temp2);
        intersect_allele(temp1,temp2, out);
    }
}

struct Node_Condenser:public tbb::task{
    Mutation_Annotated_Tree::Tree::condensed_node_t & condensed_nodes;
    const std::vector<std::string>& missing_samples;
    Mutation_Annotated_Tree::Node* root;
    std::atomic<size_t>& condensed_nodes_count;
    Node_Condenser(Mutation_Annotated_Tree::Tree::condensed_node_t & removed_nodes,const std::vector<std::string>& missing_samples,Mutation_Annotated_Tree::Node* root,std::atomic<size_t>& condensed_nodes_count):condensed_nodes(removed_nodes),missing_samples(missing_samples),root(root),condensed_nodes_count(condensed_nodes_count){}
    tbb::task* execute() override{
        std::vector<Mutation_Annotated_Tree::Node*> polytomy_nodes;
        std::vector<tbb::task*> to_spawn;
        tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
        std::vector<Mutation_Annotated_Tree::Node*> new_children;
        new_children.reserve(root->children.size());
        for(auto node:root->children){
            if (node->is_leaf()) {
                if (std::find(missing_samples.begin(), missing_samples.end(), node->identifier) == missing_samples.end()
                &&node->no_valid_mutation()) {
                    polytomy_nodes.push_back(node);
                }else {
                    new_children.push_back(node);
                }
            }else {
                to_spawn.push_back(new (
                    empty->allocate_child()) Node_Condenser(condensed_nodes,missing_samples,node,condensed_nodes_count));
                new_children.push_back(node);
            }
        }
        empty->set_ref_count(to_spawn.size());
        for(auto task:to_spawn){
            empty->spawn(*task);
        }
        if (polytomy_nodes.size()>1) {
            Mutation_Annotated_Tree::Mutations_Collection new_shared_muts;
            get_common_mutations(polytomy_nodes.begin(), polytomy_nodes.end(), new_shared_muts);
#ifdef CHECK_CONDENSE
            for(const auto node:polytomy_nodes){
                for(const auto& mut:new_shared_muts){
                    auto iter=node->mutations.find(mut);
                    assert(iter!=node->mutations.end());
                    assert(iter->get_all_major_allele()&mut.get_all_major_allele());
                }
            }
#endif
            std::vector<std::string> children_id{std::move(polytomy_nodes[0]->identifier)};
            children_id.reserve(polytomy_nodes.size());
            new_children.push_back(polytomy_nodes[0]);
            polytomy_nodes[0]->identifier="node_" + std::to_string(++condensed_nodes_count) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
            polytomy_nodes[0]->mutations.swap(new_shared_muts);
            for (size_t replaced_child_idx=1; replaced_child_idx<polytomy_nodes.size(); replaced_child_idx++) {
                children_id.push_back(std::move(polytomy_nodes[replaced_child_idx]->identifier));
                delete polytomy_nodes[replaced_child_idx];
            }
            //assert(res.second);
            root->children.swap(new_children);
            std::vector<std::string> temp;
            auto condensed_node_insert_result=condensed_nodes.emplace(polytomy_nodes[0]->identifier,temp);
            while (!condensed_node_insert_result.second) {
                polytomy_nodes[0]->identifier="node_" + std::to_string(++condensed_nodes_count) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
                condensed_node_insert_result=condensed_nodes.emplace(polytomy_nodes[0]->identifier,temp);
            }
            condensed_node_insert_result.first->second=std::move(children_id);
        }
        return to_spawn.empty()?empty:nullptr;
    }
};
void Mutation_Annotated_Tree::Tree::condense_leaves(std::vector<std::string> missing_samples) {
    if (condensed_nodes.size() > 0) {
        fprintf(stderr, "WARNING: tree contains condensed nodes. It may be condensed already!\n");
    }
    std::atomic<size_t> condensed_nodes_count(0);
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Node_Condenser(condensed_nodes,missing_samples,root,condensed_nodes_count));
    //assert(condensed_nodes_count.load()==condensed_nodes.size());
    for(const auto& condensed:condensed_nodes){
        Mutation_Annotated_Tree::Node* ori_node=get_node(condensed.second[0]);
        all_nodes.emplace(ori_node->identifier,ori_node);
        for(auto node_id:condensed.second){
            assert(condensed_nodes.count(node_id)==0);
            all_nodes.erase(node_id);
        }

    }
}