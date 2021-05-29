#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "apply_move/apply_move.hpp"
#include "tbb/parallel_for_each.h"
#include <vector>

int main(int argc,char** argv){
    MAT::Tree t;
    t.load_detatiled_mutations(argv[1]);
    auto bfs_ordered_nodes = t.breadth_first_expansion();
    Original_State_t origin_states;
    check_samples(t.root, origin_states, &t);
    // save_final_tree(t,origin_states,"tttttt");
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
                    mutated->emplace(sample.first,
                                     iter->get_all_major_allele());
                }
            }
    });    
    std::unordered_set<std::string> changed_nodes;
    clean_tree(t, changed_nodes);
}