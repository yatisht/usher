#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <algorithm>
#include <climits>
#include <mutex>
#include <tbb/task_group.h>
#include <unordered_map>
#include <vector>
#include <signal.h>
namespace MAT = Mutation_Annotated_Tree;
static void add_mut(std::unordered_map<int, uint8_t> &output,
                    const MAT::Mutation &mut) {
    auto result = output.emplace(mut.get_position(), mut.get_mut_one_hot());
    if (!result.second) {
        result.first->second |= mut.get_mut_one_hot();
    }
}
void assign_descendant_possible_muts_cont(
    std::unordered_map<int, uint8_t> &output,
    const std::vector<std::unordered_map<int, uint8_t>> &childen_out,
    MAT::Node *root) {
    auto reserve_size = root->mutations.size();
    for (const auto &child_mut : childen_out) {
        reserve_size += child_mut.size();
    }
    output.reserve(reserve_size);
    for (const auto &child_mut_vec : childen_out) {
        for (const auto &mut : child_mut_vec) {
            auto result = output.emplace(mut.first, mut.second);
            if (!result.second) {
                result.first->second |= mut.second;
            }
        }
    }
    for (auto &mut : root->mutations) {
        auto result = output.emplace(mut.get_position(), mut.get_mut_one_hot());
        if (!result.second) {
            result.first->second |= mut.get_mut_one_hot();
        }
        mut.set_descendant_mut(result.first->second);
    }
}

void assign_descendant_possible_muts_recursive(
    MAT::Node *root,
    std::unordered_map<int, uint8_t> &output,
    tbb::task_group &tg) {
    
    std::vector<std::unordered_map<int, uint8_t>> childen_out(root->children.size());
    
    for (size_t idx = 0; idx < root->children.size(); idx++) {
        auto this_child = root->children[idx];
        if (this_child->children.empty()) {
            childen_out[idx].reserve(this_child->mutations.size());
            for (auto &mut : this_child->mutations) {
                mut.set_descendant_mut(mut.get_mut_one_hot());
                childen_out[idx].emplace(mut.get_position(), mut.get_mut_one_hot());
            }
        } else {
            tg.run([this_child, &childen_out, idx, &tg]() {
                assign_descendant_possible_muts_recursive(this_child, childen_out[idx], tg);
            });
        }
    }
    
    tg.wait();
    assign_descendant_possible_muts_cont(output, childen_out, root);
}
void assign_descendant_muts(MAT::Tree &in) {
    std::unordered_map<int, uint8_t> ignore;
    tbb::task_group tg;
    assign_descendant_possible_muts_recursive(in.root, ignore, tg);
}
