#include "usher.hpp"
#include <algorithm>
#include <climits>
#include <mutex>
#include <tbb/task.h>
#include <unordered_map>
#include <vector>
#include <signal.h>
namespace MAT = Mutation_Annotated_Tree;
int check_split(const MAT::Node *root, std::vector<bool>& do_split,
                int threshold,bool skip_leaf) {
    if (root->children.empty()) {
        return 0;
    }
    int max_dist = 0;
    const MAT::Node *max_dist_node = nullptr;
    const MAT::Node *second_max_dist_node = nullptr;
    int second_max_dist = 0;
    for (const MAT::Node *c : root->children) {
        auto child_dist =
            check_split(c, do_split, threshold,skip_leaf) + c->branch_length;
        if (child_dist>=threshold&&(!(skip_leaf&&c->children.empty()))) {
            do_split[c->dfs_index]=true;
            child_dist=c->branch_length;
        }
        if (child_dist >= max_dist) {
            second_max_dist = max_dist;
            second_max_dist_node=max_dist_node;
            max_dist = child_dist;
            max_dist_node = c;
        }else if (child_dist>second_max_dist) {
            second_max_dist=child_dist;
            second_max_dist_node=c;
        }
    }
    if (max_dist + second_max_dist >= threshold) {
        if (max_dist_node->branch_length+second_max_dist>=threshold||(skip_leaf&&max_dist_node->children.empty())) {
            do_split[root->dfs_index]=true;
            return 0;
        }
        do_split[max_dist_node->dfs_index] = true;
        return std::max((int)max_dist_node->branch_length, second_max_dist);
    }
    return max_dist;
}
struct Sampled_Mut_Cont : public tbb::task {
    std::vector<Sampled_Tree_Mutation> parent_muts;
    tbb::task *execute() override { return nullptr; }
};
struct Sample_Mut_Op : public tbb::task {
    const std::vector<Sampled_Tree_Mutation> &parent_mutations;
    const std::vector<bool> &is_sampled;
    const MAT::Node *root;
    std::shared_ptr<std::mutex> new_tree_mutex;
    Sampled_Tree_Node *new_tree_root;
    Sample_Mut_Op(const std::vector<Sampled_Tree_Mutation> &parent_mutations,
                  const std::vector<bool> &is_sampled, const MAT::Node *root,
                  std::shared_ptr<std::mutex> new_tree_mutex,
                  Sampled_Tree_Node *new_tree_root)
        : parent_mutations(parent_mutations), is_sampled(is_sampled),
          root(root), new_tree_mutex(new_tree_mutex),
          new_tree_root(new_tree_root) {}
    tbb::task *execute() override {
        auto cont = new (allocate_continuation()) Sampled_Mut_Cont;
        cont->set_ref_count(root->children.size());
        auto iter = parent_mutations.begin();
        for (const auto &mut : root->mutations) {
            while (iter->position < mut.get_position()) {
                cont->parent_muts.push_back(*iter);
                iter++;
            }
            if (iter->position == mut.get_position()) {
                if(iter->par_nuc!=mut.get_mut_one_hot()){
                    cont->parent_muts.push_back(*iter);
                    cont->parent_muts.back().mut_nuc = mut.get_mut_one_hot();
                    cont->parent_muts.back().descendent_possible_nuc = mut.get_mut_one_hot();
                }
                iter++;
            } else {
                Sampled_Tree_Mutation temp;
                temp.chrom_idx = mut.get_chromIdx();
                temp.position = mut.get_position();
                temp.descendent_possible_nuc = mut.get_mut_one_hot();
                temp.mut_nuc = mut.get_mut_one_hot();
                temp.par_nuc=mut.get_par_one_hot();
                cont->parent_muts.push_back(temp);
            }
//            assert(cont->parent_muts.empty()||cont->parent_muts.back().par_nuc!=cont->parent_muts.back().mut_nuc);
        }
        cont->parent_muts.insert(cont->parent_muts.end(), iter,
                                 parent_mutations.end());
        assert(cont->parent_muts.back().position==INT_MAX);
        std::shared_ptr<std::mutex> new_tree_mutex_next = new_tree_mutex;
        Sampled_Tree_Node *new_tree_root_next = new_tree_root;
        if (is_sampled[root->dfs_index]) {
            new_tree_root_next = new Sampled_Tree_Node;
            new_tree_root_next->corresponding_main_tree_node = root;
            new_tree_root_next->mutations.swap(cont->parent_muts);
            Sampled_Tree_Mutation temp;
            temp.position=INT_MAX;
            cont->parent_muts.push_back(temp);
            assert(new_tree_root_next->mutations.back().position==INT_MAX);
            new_tree_root_next->mutations.pop_back();
            new_tree_root_next->parent = new_tree_root;
            {
                std::lock_guard<std::mutex> lk(*new_tree_mutex);
                new_tree_root->children.push_back(new_tree_root_next);
            }
            new_tree_mutex_next.reset(new std::mutex);
        }
        for (const auto c : root->children) {
            cont->spawn(*(new (cont->allocate_child()) Sample_Mut_Op(
                cont->parent_muts, is_sampled, c, new_tree_mutex_next,
                new_tree_root_next)));
        }
        return root->children.empty() ? cont : nullptr;
    }
};
static void add_mut(std::unordered_map<int, uint8_t> &output,
                    const Sampled_Tree_Mutation &mut) {
    auto result = output.emplace(mut.position, mut.mut_nuc);
    if (!result.second) {
        result.first->second |= mut.mut_nuc;
    }
}
struct Assign_Descendant_Possible_Muts_Cont : public tbb::task {
    std::unordered_map<int, uint8_t> &output;
    std::vector<std::unordered_map<int, uint8_t>> childen_out;
    Sampled_Tree_Node *root;
    Assign_Descendant_Possible_Muts_Cont(
        std::unordered_map<int, uint8_t> &output, size_t child_count,
        Sampled_Tree_Node *root)
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
            auto result = output.emplace(mut.position, mut.mut_nuc);
            if (!result.second) {
                result.first->second |= mut.mut_nuc;
            }
            mut.descendent_possible_nuc = result.first->second;
        }
        return nullptr;
    }
};

struct Assign_Descendant_Possible_Muts : public tbb::task {
    Sampled_Tree_Node *root;
    std::unordered_map<int, uint8_t> &output;
    Assign_Descendant_Possible_Muts(Sampled_Tree_Node *root,
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
                for (const auto &mut : this_child->mutations) {
                    cont->childen_out[idx].emplace(mut.position, mut.mut_nuc);
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
Sampled_Tree_Node *sample_tree(MAT::Tree &in, int threshold,bool avoid_leaf) {
    auto dfs = in.depth_first_expansion();
    std::vector<bool> is_sampled(dfs.size(), false);
    check_split(in.root, is_sampled, threshold,avoid_leaf);
    Sampled_Tree_Node *sample_tree_root = new Sampled_Tree_Node;
    sample_tree_root->corresponding_main_tree_node = in.root;
    sample_tree_root->parent = nullptr;
    std::vector<Sampled_Tree_Mutation> parent_mutations;
    Sampled_Tree_Mutation temp;
    temp.position=INT_MAX;
    parent_mutations.push_back(temp);
    tbb::task::spawn_root_and_wait(
        *new (tbb::task::allocate_root()) Sample_Mut_Op(
            parent_mutations, is_sampled, in.root,
            std::shared_ptr<std::mutex>(new std::mutex), sample_tree_root));
    std::unordered_map<int, uint8_t> ignore;
    tbb::task::spawn_root_and_wait(
        *new (tbb::task::allocate_root())
            Assign_Descendant_Possible_Muts(sample_tree_root, ignore));
    return sample_tree_root;
}
void sample_tree_dfs(Sampled_Tree_Node *sampled_tree_root,std::vector<Sampled_Tree_Node *>& output){
    sampled_tree_root->dfs_idx=output.size();
    output.push_back(sampled_tree_root);
    for (auto child : sampled_tree_root->children) {
        sample_tree_dfs(child, output);
    }
}
void remove_sampled_tree(Sampled_Tree_Node *sampled_tree_root){
    for (auto child : sampled_tree_root->children) {
        remove_sampled_tree(child);
    }
    delete sampled_tree_root;
}