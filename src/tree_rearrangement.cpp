#include "tree_rearrangement.hpp"
#include "check_samples.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iterator>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for_each.h>
#include <src/Twice_Bloom_Filter.hpp>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;

static void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output){
    Twice_Bloom_Filter filter;
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            filter.insert(m.position);
        }
    }
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            if(filter.query(m.position)){
                output.push_back(n);
                break;
            }
        }
    }
}
static void split_rounds(std::vector<MAT::Node *> &this_round,
                         std::vector<MAT::Node *> &next_round,
                         std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> idx;

    {
        std::vector<size_t> nidx;
        nidx.reserve(this_round.size());
        for (auto node : this_round) {
            auto parent = node->parent;
            if (parent)
                nidx.push_back(parent->index);
        }
        std::sort(nidx.begin(), nidx.end());

        std::vector<size_t> deferred_idx;
        deferred_idx.reserve(next_round.size());
        for (auto node : next_round) {
            deferred_idx.push_back(node->index);
        }
        idx.reserve(nidx.size() + deferred_idx.size());
        merge(nidx.begin(), nidx.end(), deferred_idx.begin(), deferred_idx.end(),
              std::back_inserter(idx));
    }
    this_round.clear();
    this_round.reserve(idx.size());
    next_round.clear();
    next_round.reserve(idx.size());
    if(idx.empty()) return;
    for (auto iter = idx.begin(); iter < idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (!check_grand_parent(dfs_ordered_nodes[next_ele]->parent,dfs_ordered_nodes[*iter]->parent)) {
            this_round.push_back(dfs_ordered_nodes[*iter]);
        } else {
            next_round.push_back(dfs_ordered_nodes[*iter]);
        }
    }
    this_round.push_back(dfs_ordered_nodes[idx.back()]);
}
void Tree_Rearrangement::refine_trees(
    std::vector<MAT::Tree> &optimal_trees) {

    for (auto this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Sample_Mut_Type ori;
        std::vector<MAT::Node*>& new_nodes=this_tree.new_nodes;
        while(!new_nodes.empty()){
        check_samples(this_tree.root, ori);
        std::vector<MAT::Node *> optimized;
        std::vector<MAT::Node *>& this_round = new_nodes;;
        std::vector<MAT::Node *> next_round;
        find_nodes_with_recurrent_mutations(dfs_ordered_nodes, new_nodes);
        split_rounds(this_round, next_round, dfs_ordered_nodes);
        tbb::mutex mutex;
        while (!this_round.empty()) {
            bool have_update = false;
            //fprintf(stderr, "Optimizing %zu nodes \n",this_round.size());
            tbb::parallel_for(tbb::blocked_range<size_t>(0,this_round.size(),10),[&](tbb::blocked_range<size_t> r){
                for (size_t i=r.begin();i<r.end() ; i++) {
                    if(Tree_Rearrangement::move_nearest(
                    this_round[i], dfs_ordered_nodes, this_tree)){
                        std::lock_guard<tbb::mutex> lock(mutex);
                        have_update=true;
                        optimized.push_back(this_round[i]);
                    }
                }
            });
            
            if (have_update) {
                Sample_Mut_Type copy(ori);
                check_samples(this_tree.root, copy);
                dfs_ordered_nodes = this_tree.depth_first_expansion();
            }
        split_rounds(this_round, next_round, dfs_ordered_nodes);
        }
        new_nodes=optimized;
        }

        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        this_tree.reassign_level();
    }
}
