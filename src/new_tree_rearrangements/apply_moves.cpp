#include "mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "Fitch_Sankoff.hpp"
#include "tree_rearrangement_internal.hpp"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include <algorithm>
#include <string>
namespace MAT = Mutation_Annotated_Tree;
    int get_parsimmony_score_only(MAT::Node* src, MAT::Node* dst,MAT::Node* start_node);
/*
static void acc_counts(std::array<int, 4> &count,nuc_one_hot majority_nuc){
    
}
bool get_children_alleles_count(const MAT::Node *this_node,
                                       const std::vector<MAT::Node *>overriden_nodes,
                                       std::array<int, 4> &count,
                                       int position,nuc_one_hot this_allele,int& mutation_count) {
    int par_allele_count = 0;
    for (MAT::Node *child : this_node->children) {
        auto mut_iter = child->mutations.find(position);
        auto end_iter=child->mutations.end();
                    if(mut_iter!=end_iter){
                if(mut_iter->is_valid()){
                    mutation_count++;
                }
        }
        if (std::find(overriden_nodes.begin(),overriden_nodes.end(),child)!=overriden_nodes.end()) {
            continue;
        }
        if (mut_iter == end_iter) {
            par_allele_count++;
        } else {
            nuc_one_hot majority_nuc =
                mut_iter->get_mut_one_hot() | mut_iter->get_tie_one_hot();
            for (int i = 0; i < 4; i++) {
                if (majority_nuc & (1 << i)) {
                    count[i]++;
                }
            }
        }
    }
    count[one_hot_to_two_bit(this_allele)] += par_allele_count;
}*/

static void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,tbb::concurrent_vector<MAT::Node*>& to_filter){

    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::vector<MAT::Node *> this_node_ori_children = this_node->children;

    if (this_node->children.size()==1&&tree.condensed_nodes.count(this_node->identifier)==0) {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        tree.all_nodes.erase(this_node->identifier);
        for (MAT::Node *child : this_node_ori_children) {
            child->parent = this_node->parent;
            parent_children.push_back(child);
        }
        
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        clean_up_internal_nodes(child, tree,to_filter);
    }
}

void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,tbb::concurrent_vector<MAT::Node*>& to_filter) {
    std::vector<MAT::Node*> removed;
    /*for(auto m:all_moves){
        fprintf(stderr, "src: %s, dst %s,LCA %s\n",m->get_src()->identifier.c_str(),m->get_dst()->identifier.c_str(),m->LCA->identifier.c_str());
        int score_ref=get_parsimmony_score_only(m->get_src(), m->get_dst(),m->dst_to_LCA.back());
        if(m->score_change!=score_ref){
            fprintf(stderr,"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nrecorded score_change %d, reference %d\n\n",m->score_change,score_ref);
            individual_move(m->get_src(),m->get_dst(),m->LCA);
            get_parsimmony_score_only(m->get_src(), m->get_dst(),m->dst_to_LCA.back());
        }
    }*/
    for (auto m : all_moves) {
        MAT::Node* src=m->get_src();
        MAT::Node* dst=m->get_dst();
        MAT::Node* LCA=m->LCA;

        Mutation_Annotated_Tree::remove_child(src,removed);

        if (dst==LCA) {
            assert(m->src_to_LCA.back()->parent==LCA);
            MAT::Node* src_branch_node=m->src_to_LCA.back();
            auto& LCA_children=LCA->children;
            auto iter=std::find(LCA_children.begin(),LCA_children.end(),src_branch_node);
            assert(iter!=LCA_children.end());
            LCA_children.erase(iter);
            MAT::Node* new_src_branch_node=t.create_node(std::to_string(++t.curr_internal_node),LCA);
            new_src_branch_node->children.push_back(src_branch_node);
            src_branch_node->parent=new_src_branch_node;
            new_src_branch_node->children.push_back(src);
            src->parent=new_src_branch_node;
        }else {
            MAT::Node* dst_parent=dst->parent;
            MAT::Node* new_dst_branch=t.create_node(std::to_string(++t.curr_internal_node),dst_parent);
            auto& dst_parent_children=dst_parent->children;
            auto iter=std::find(dst_parent_children.begin(),dst_parent_children.end(),dst);
            assert(iter!=dst_parent_children.end());
            dst_parent_children.erase(iter);
            new_dst_branch->children.push_back(dst);
            dst->parent=new_dst_branch;
            new_dst_branch->children.push_back(src);
            src->parent=new_dst_branch;
        }
        delete m;
    }
    if (!removed.empty()) {
        fprintf(stderr, "removed %zu nodes \n",removed.size());
        tbb::concurrent_vector<MAT::Node*> filtered;
        for(auto node:to_filter){
            if (std::find(removed.begin(),removed.end(),node)==removed.end()) {
                filtered.push_back(node);
            }
        }
        to_filter=std::move(filtered);
    }
    clean_up_internal_nodes(t.root->children.size()==1?t.root->children[0]:t.root, t, to_filter);
    bfs_ordered_nodes = t.breadth_first_expansion();

    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
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