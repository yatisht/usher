#include "mutation_annotated_tree.hpp"
#include "tbb/task.h"
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <utility>
#include <vector>
//find mutated alleles shared by both of the input
static void intersect_allele(const Mutation_Annotated_Tree::Mutations_Collection& in1, const Mutation_Annotated_Tree::Mutations_Collection& in2, Mutation_Annotated_Tree::Mutations_Collection& out) {
    if (in1.empty()||in2.empty()) {
        return;
    }
    auto in1_iter=in1.begin();
    auto in1_end=in1.end();
    for(const auto& mut:in2) {
        while (in1_iter!=in1_end&&in1_iter->get_position()<mut.get_position()) {
            in1_iter++;
            //skip all muts present in only one of the inout
        }
        if (in1_iter!=in1_end&&in1_iter->get_position()==mut.get_position()) {
            //output common major allele
            nuc_one_hot common=in1_iter->get_all_major_allele()&mut.get_all_major_allele();
            if (common!=mut.get_par_one_hot()) {
                out.push_back(mut);
                out.back().set_auxillary(common, 0xf&(~common));
            }
            in1_iter++;
        }
    }
}
//merge sortish way to find major alleles shared by all nodes from begin to end
static void get_common_mutations(std::vector<Mutation_Annotated_Tree::Node*>::const_iterator begin,std::vector<Mutation_Annotated_Tree::Node*>::const_iterator end,Mutation_Annotated_Tree::Mutations_Collection& out) {
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
        get_common_mutations(begin,middle, temp1);
        get_common_mutations(middle,end,temp2);
        intersect_allele(temp1,temp2, out);
    }
}
//Condense nodes, this parallelization is more
//for getting around stack overflow than performance,
//as condensing nodes doesn't take too long anyway
struct Node_Condenser {
    const std::vector<Mutation_Annotated_Tree::Node*>& condensable_nodes;
    tbb::concurrent_vector<std::vector<size_t>>& condensed_ids;
    void condense(Mutation_Annotated_Tree::Node* root) const {
        //condensed children of root
        std::vector<Mutation_Annotated_Tree::Node*> polytomy_nodes;
        //instead of removing,just copy children over
        std::vector<Mutation_Annotated_Tree::Node*> new_children;
        new_children.reserve(root->children.size());
        for(auto node:root->children) {
            if (node->is_leaf()&&node->no_valid_mutation()) {
                //have no valid mutation, so condense it
                polytomy_nodes.push_back(node);
            } else {
                new_children.push_back(node);
            }
        }
        if (polytomy_nodes.size()>1) {
            Mutation_Annotated_Tree::Mutations_Collection new_shared_muts;
            //get all ambiguous mutation, to place on the condensed node
            get_common_mutations(polytomy_nodes.begin(), polytomy_nodes.end(), new_shared_muts);
#ifdef CHECK_CONDENSE
            for(const auto node:polytomy_nodes) {
                for(const auto& mut:new_shared_muts) {
                    auto iter=node->mutations.find(mut);
                    assert(iter!=node->mutations.end());
                    assert(iter->get_all_major_allele()&mut.get_all_major_allele());
                }
            }
#endif
            //replace the first node to condense with the condensed node in place
            std::vector<size_t> condensed_node_group_ids{polytomy_nodes[0]->node_id};
            condensed_node_group_ids.reserve(polytomy_nodes.size());
            new_children.push_back(polytomy_nodes[0]);
            //get a new name for condensed node
            polytomy_nodes[0]->mutations.swap(new_shared_muts);
            //record original sample name
            for (size_t replaced_child_idx=1; replaced_child_idx<polytomy_nodes.size(); replaced_child_idx++) {
                condensed_node_group_ids.push_back(polytomy_nodes[replaced_child_idx]->node_id);
                delete polytomy_nodes[replaced_child_idx];
            }
            //assert(res.second);
            root->children.swap(new_children);
            condensed_ids.emplace_back(std::move(condensed_node_group_ids));
        }
    }
    void operator()(tbb::blocked_range<size_t> r) const {
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            condense(condensable_nodes[idx]);
        }
    }
};
void find_condensable_nodes(Mutation_Annotated_Tree::Node* root,std::vector<Mutation_Annotated_Tree::Node*>& condensable_nodes) {
    size_t child_count=0;
    for (auto child : root->children) {
        if (child->is_leaf()) {
            child_count++;
        } else {
            find_condensable_nodes(child, condensable_nodes);
        }
    }
    if (child_count>1) {
        condensable_nodes.push_back(root);
    }
}
void Mutation_Annotated_Tree::Tree::condense_leaves(std::vector<std::string> missing_samples) {
    if (condensed_nodes.size() > 0) {
        fprintf(stderr, "WARNING: tree contains condensed nodes. It may be condensed already!\n");
    }
    std::vector<Mutation_Annotated_Tree::Node*> condensable_nodes;
    find_condensable_nodes(root, condensable_nodes);
    tbb::concurrent_vector<std::vector<size_t>> condensed_ids;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,condensable_nodes.size()),Node_Condenser{condensable_nodes,condensed_ids});
    //assert(condensed_nodes_count.load()==condensed_nodes.size());
    size_t condensed_nodes_count(0);
    for(const auto& condensed:condensed_ids) {
        auto new_id="node_" + std::to_string(++condensed_nodes_count) + "_condensed_" + std::to_string(condensed.size()) + "_leaves";
        auto preserved_id=condensed[0];
        std::vector<std::string> sample_names;
        sample_names.reserve(condensed.size());
        for (auto old_ids : condensed) {
            auto iter=node_names.find(old_ids);
            if(iter==node_names.end()) {
                fprintf(stderr,"old node id %zu not found, condensed vector addr %lx \n",old_ids,(unsigned long)&condensed);
                raise(SIGTRAP);
                continue;
            }
            sample_names.push_back(iter->second);
            auto back_iter=node_name_to_idx_map.find(iter->second);
            node_names.erase(iter);
            if(back_iter!=node_name_to_idx_map.end()) {
                node_name_to_idx_map.erase(back_iter);
            }
        }
        for (size_t idx=1; idx<condensed.size(); idx++) {
            all_nodes[condensed[idx]]=nullptr;
        }
        condensed_nodes.emplace(preserved_id,std::move(sample_names));
        node_names.emplace(preserved_id,new_id);
        node_name_to_idx_map.emplace(new_id,preserved_id);
    }
}
