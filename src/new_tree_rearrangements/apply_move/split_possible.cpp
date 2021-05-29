#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <vector>
static void find_back_mut();
static void find_split(std::vector<MAT::Node *> &all_nodes,
                       std::vector<MAT::Node *> &to_split) {
    for (auto node : all_nodes) {
        if (node->parent) {
            for (const auto mut : node->mutations) {
                auto &par_mutations = node->parent->mutations;
                auto iter = par_mutations.find(mut.get_position());
                if (iter != par_mutations.end()) {
                    if (mut.is_valid()&&iter->get_par_one_hot()==mut.get_mut_one_hot()) {
                        to_split.push_back(node);
                    }
                }
            }
        }
    }
}
