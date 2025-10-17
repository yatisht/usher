#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <algorithm>
#include <climits>
#include <signal.h>
#include <unordered_map>
#include <vector>
#include <taskflow/taskflow.hpp>

namespace MAT = Mutation_Annotated_Tree;

static void add_mut(std::unordered_map<int, uint8_t> &output,
                    const MAT::Mutation &mut) {
    auto result = output.emplace(mut.get_position(), mut.get_mut_one_hot());
    if (!result.second) {
        result.first->second |= mut.get_mut_one_hot();
    }
}

struct Assign_Descendant_Possible_Muts_Cont {
    std::unordered_map<int, uint8_t> &output;
    std::vector<std::unordered_map<int, uint8_t>> children_out;
    MAT::Node *root;
    Assign_Descendant_Possible_Muts_Cont(
        std::unordered_map<int, uint8_t> &output, size_t child_count,
        MAT::Node *root)
        : output(output), children_out(child_count), root(root) {}

    void execute() {
        auto reserve_size = root->mutations.size();
        for (const auto &child_mut : children_out) {
            reserve_size += child_mut.size();
        }
        output.reserve(reserve_size);
        for (const auto &child_mut_vec : children_out) {
            for (const auto &mut : child_mut_vec) {
                auto result = output.emplace(mut.first, mut.second);
                if (!result.second) {
                    result.first->second |= mut.second;
                }
            }
        }
        for (auto &mut : root->mutations) {
            auto result =
                output.emplace(mut.get_position(), mut.get_mut_one_hot());
            if (!result.second) {
                result.first->second |= mut.get_mut_one_hot();
            }
            mut.set_descendant_mut(result.first->second);
        }
    }
};

struct Assign_Descendant_Possible_Muts {
    tf::Taskflow& taskflow;
    MAT::Node *root;
    std::unordered_map<int, uint8_t> &output;
    Assign_Descendant_Possible_Muts(tf::Taskflow& t, MAT::Node *root,
                                    std::unordered_map<int, uint8_t> &output)
        : taskflow{t}, root(root), output(output) {}
    void execute() const {
        Assign_Descendant_Possible_Muts_Cont cont(output, root->children.size(),
                                                  root);
        std::vector<Assign_Descendant_Possible_Muts> children_tasks;
        children_tasks.reserve(root->children.size());
        for (size_t idx = 0; idx < root->children.size(); idx++) {
            auto this_child = root->children[idx];
            if (this_child->children.empty()) {
                cont.children_out[idx].reserve(this_child->mutations.size());
                for (auto &mut : this_child->mutations) {
                    mut.set_descendant_mut(mut.get_mut_one_hot());
                    cont.children_out[idx].emplace(mut.get_position(),
                                                   mut.get_mut_one_hot());
                }
            } else {
                children_tasks.emplace_back(taskflow,
                                            this_child, cont.children_out[idx]);
            }
        }
        if (children_tasks.empty()) {
            cont.execute();
        } else {
          taskflow.emplace([&children_tasks](tf::Subflow& subflow){
            for (auto&& child_task : std::move(children_tasks)) {
                subflow.emplace([child_task = std::move(child_task)] {
                    child_task.execute();
                });
            }
          });
        }
    }
};

void assign_descendant_muts(MAT::Tree &in) {
    std::unordered_map<int, uint8_t> ignore;
    tf::Executor executor;
    tf::Taskflow taskflow;
    taskflow.emplace([&] {
      Assign_Descendant_Possible_Muts(taskflow, in.root, ignore).execute();
    });
    executor.run(taskflow).wait(); 
}
