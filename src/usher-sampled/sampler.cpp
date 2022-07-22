#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <algorithm>
#include <climits>
#include <mutex>
#include <tbb/task.h>
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
struct Assign_Descendant_Possible_Muts_Cont : public tbb::task {
    std::unordered_map<int, uint8_t> &output;
    std::vector<std::unordered_map<int, uint8_t>> childen_out;
    MAT::Node *root;
    Assign_Descendant_Possible_Muts_Cont(
        std::unordered_map<int, uint8_t> &output, size_t child_count,
        MAT::Node *root)
        : output(output), childen_out(child_count), root(root) {}

    tbb::task *execute() {
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
        return nullptr;
    }
};

struct Assign_Descendant_Possible_Muts : public tbb::task {
    MAT::Node *root;
    std::unordered_map<int, uint8_t> &output;
    Assign_Descendant_Possible_Muts(MAT::Node *root,
                                    std::unordered_map<int, uint8_t> &output)
        : root(root), output(output) {}
    tbb::task *execute() {
        auto cont = new (allocate_continuation())
        Assign_Descendant_Possible_Muts_Cont(output, root->children.size(),
                                             root);
        std::vector<Assign_Descendant_Possible_Muts *> children_tasks;
        children_tasks.reserve(root->children.size());
        for (size_t idx = 0; idx < root->children.size(); idx++) {
            auto this_child = root->children[idx];
            if (this_child->children.empty()) {
                cont->childen_out[idx].reserve(this_child->mutations.size());
                for (auto &mut : this_child->mutations) {
                    mut.set_descendant_mut(mut.get_mut_one_hot());
                    cont->childen_out[idx].emplace(mut.get_position(), mut.get_mut_one_hot());
                }
            } else {
                auto new_task = new (cont->allocate_child())
                Assign_Descendant_Possible_Muts(this_child,
                                                cont->childen_out[idx]);
                children_tasks.push_back(new_task);
            }
        }
        cont->set_ref_count(children_tasks.size());
        for (auto child_task : children_tasks) {
            cont->spawn(*child_task);
        }
        return children_tasks.empty()?cont:nullptr;
    }
};
void assign_descendant_muts(MAT::Tree &in) {
    std::unordered_map<int, uint8_t> ignore;
    tbb::task::spawn_root_and_wait(
        *new (tbb::task::allocate_root())
        Assign_Descendant_Possible_Muts(in.root,ignore));
}
