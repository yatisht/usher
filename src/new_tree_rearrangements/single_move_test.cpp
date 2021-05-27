#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstdio>
#include <utility>
#include "tbb/parallel_for_each.h"
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst){
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src=src->parent;
        }
        else if (src->dfs_index<dst->dfs_index) {
            dst=dst->parent;
        }
    }
    return src;
}
int main(int argc,char** argv){
    MAT::Tree t;
    t.load_detatiled_mutations("2021-01-20-msa-ambiguious-tree-patch-itermediate14.pb");
    Original_State_t origin_states;
    check_samples(t.root, origin_states, &t);
    auto bfs_ordered_nodes = t.breadth_first_expansion();
    for (MAT::Node *node : bfs_ordered_nodes) {
        for (const MAT::Mutation &m : node->mutations) {
            mutated_positions.emplace(
                m, new std::unordered_map<std::string, nuc_one_hot>);
        }
        node->tree = &t;
    }
        tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&origin_states](
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
        });
    std::vector<Profitable_Moves_ptr_t> all_moves{};
    fprintf(stderr,"%zu",t.get_parsimony_score());
                        output_t out;
    tbb::concurrent_vector<MAT::Node*> nodes_to_search;
        
nodes_to_search.push_back(t.get_node("USA/WA-UW-5947/2020|MT412289.1|20-04-07"));
nodes_to_search.push_back(t.get_node("USA/WA-UW-5170/2020|MT375482.1|20-04-02"));
nodes_to_search.push_back(t.get_node("USA/WA-UW-5717/2020|MT412270.1|20-04-06"));
nodes_to_search.push_back(t.get_node("USA/WA-UW-1724/2020|MT326154.1|20-03-21"));
nodes_to_search.push_back(t.get_node("USA/WA-UW-4734/2020|MT375446.1|20-03-28"));
nodes_to_search.push_back(t.get_node("6191"));
nodes_to_search.push_back(t.get_node("25832"));
nodes_to_search.push_back(t.get_node("39946"));
nodes_to_search.push_back(t.get_node("England/PHEC-16C12/2020|20-03-30"));
nodes_to_search.push_back(t.get_node("TUN/TUN202055304/2020|MW279303.1|20-10-31"));
optimize_tree(bfs_ordered_nodes, nodes_to_search, t, origin_states,10
            #ifdef CONFLICT_RESOLVER_DEBUG
            ,stderr
            #endif
            );
}
/*
change of -1 @ 24389 
change of -1 @ 24390 
change of 1 @ 28881 
change of 1 @ 28882 
change of 1 @ 28883 
*/