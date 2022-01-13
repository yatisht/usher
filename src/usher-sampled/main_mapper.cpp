#include "mapper.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <bits/floatn-common.h>
#include <cassert>
#include <climits>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <signal.h>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#include <unordered_map>
#include <utility>
#include <vector>
template <typename Hook1, typename Hook2> struct Combine_Hook {
    Hook1 hook1;
    Hook2 hook2;
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        hook1.reserve(sample_mutation_count, target_mutation_count);
        hook2.reserve(sample_mutation_count, target_mutation_count);
    }
    void sample_mut_only(const MAT::Mutation &mut) {
        hook1.sample_mut_only(mut);
        hook2.sample_mut_only(mut);
    }
    void target_node_only(const MAT::Mutation &mut) {
        hook1.target_node_only(mut);
        hook2.target_node_only(mut);
    }
    void both(const MAT::Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        hook1.both(sample_mut, target_mut);
        hook2.both(sample_mut, target_mut);
    }
};
struct Empty_Hook {
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {}
    void sample_mut_only(const MAT::Mutation &mut) {}
    void target_node_only(const MAT::Mutation &mut) {}
    void both(const MAT::Mutation &mut, const MAT::Mutation &mut2) {}
};
static void insert_split(const MAT::Mutation &sample_mut,
                         const MAT::Mutation &target_mut,
                         MAT::Mutations_Collection &sample_mutations,
                         MAT::Mutations_Collection &splitted_mutations,
                         MAT::Mutations_Collection &shared_mutations) {
    auto par_nuc = target_mut.get_par_one_hot();
    auto sample_nuc = sample_mut.get_mut_one_hot();
    auto target_nuc = target_mut.get_mut_one_hot();
    auto coincided_mut = sample_nuc & target_nuc;
    if (coincided_mut & par_nuc) {
        coincided_mut = coincided_mut & par_nuc;
    }
    coincided_mut = 1 << __builtin_ctz(coincided_mut);
    assert(coincided_mut & sample_nuc & target_nuc);
    shared_mutations.push_back(target_mut);
    shared_mutations.back().set_par_mut(par_nuc, coincided_mut);
    if (coincided_mut != target_nuc) {
        splitted_mutations.push_back(target_mut);
        splitted_mutations.back().set_par_mut(coincided_mut, target_nuc);
    }
    if (coincided_mut != sample_nuc) {
        sample_mutations.push_back(sample_mut);
        sample_mutations.back().set_par_mut(coincided_mut, sample_nuc);
    }
}
struct Down_Sibling_Hook {
    MAT::Mutations_Collection &sample_mutations;
    MAT::Mutations_Collection &splitted_mutations;
    MAT::Mutations_Collection &shared_mutations;
    int &parsimony_score;
    Down_Sibling_Hook(Main_Tree_Target &target, int &parsimony_score)
        : sample_mutations(target.sample_mutations),
          splitted_mutations(target.splited_mutations),
          shared_mutations(target.shared_mutations),
          parsimony_score(parsimony_score) {
        sample_mutations.clear();
        splitted_mutations.clear();
        shared_mutations.clear();
        parsimony_score = 0;
    }
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        sample_mutations.reserve(sample_mutation_count);
        splitted_mutations.reserve(target_mutation_count);
        shared_mutations.reserve(target_mutation_count);
    }
    void sample_mut_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        sample_mutations.push_back(mut);
        if (!(mut.get_par_one_hot() & mut.get_mut_one_hot())) {
            parsimony_score++;
        }
    }
    void target_node_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        splitted_mutations.push_back(mut);
    }
    void both(const MAT::Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        if (sample_mut.get_mut_one_hot() != target_mut.get_mut_one_hot()) {
            if (sample_mut.get_mut_one_hot() & target_mut.get_mut_one_hot()) {
                insert_split(sample_mut, target_mut, sample_mutations,
                             splitted_mutations, shared_mutations);
            } else {
                sample_mutations.push_back(sample_mut);
                splitted_mutations.push_back(target_mut);
                parsimony_score++;
            }
        } else {
            if (__builtin_popcount(target_mut.get_mut_one_hot())!=1) {
                auto common=target_mut.get_par_one_hot()&target_mut.get_mut_one_hot();
                if (!common) {
                    common=1 << __builtin_ctz(target_mut.get_mut_one_hot());
                    splitted_mutations.push_back(target_mut);
                    splitted_mutations.back().set_mut_one_hot(common);
                }
                sample_mutations.push_back(sample_mut);
                sample_mutations.back().set_par_one_hot(common);
                splitted_mutations.push_back(target_mut);
                splitted_mutations.back().set_par_one_hot(common);
            }else {
                shared_mutations.push_back(target_mut);            
            }
        }
    }
};
struct Down_Decendant_Hook {
    MAT::Mutations_Collection &muts;
    Down_Decendant_Hook(MAT::Mutations_Collection &muts) : muts(muts) {
        muts.clear();
    }
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        muts.reserve(sample_mutation_count + target_mutation_count);
    }
    void sample_mut_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        muts.push_back(mut);
    }
    void target_node_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        muts.push_back(mut);
        auto &last_mut = muts.back();
        // back mutation
        last_mut.set_par_mut(last_mut.get_mut_one_hot(),
                             last_mut.get_par_one_hot());
    }
    void both(const MAT::Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        if (sample_mut.get_mut_one_hot() != target_mut.get_mut_one_hot()) {
            muts.push_back(sample_mut);
            auto &last_mut = muts.back();
            last_mut.set_par_one_hot(target_mut.get_mut_one_hot());
        }
        assert(muts.empty()||muts.back().get_mut_one_hot() != muts.back().get_par_one_hot());
    }
};
struct Upward_Sibling_Hook {
    MAT::Mutations_Collection &sample_mutations;
    MAT::Mutations_Collection &splitted_mutations;
    MAT::Mutations_Collection &shared_mutations;
    int &parsimony_score;
    Upward_Sibling_Hook(Main_Tree_Target &target, int &parsimony_score)
        : sample_mutations(target.sample_mutations),
          splitted_mutations(target.splited_mutations),
          shared_mutations(target.shared_mutations),
          parsimony_score(parsimony_score) {
        sample_mutations.clear();
        splitted_mutations.clear();
        shared_mutations.clear();
        parsimony_score = 0;
    }
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        sample_mutations.reserve(sample_mutation_count);
        splitted_mutations.reserve(target_mutation_count);
    }
    void sample_mut_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        // same as downward, pass through
        sample_mutations.push_back(mut);
        if (!(mut.get_mut_one_hot() & mut.get_par_one_hot())) {
            parsimony_score++;
        }
    }
    void target_node_only(const MAT::Mutation &mut) {
        // shared mutation
        assert(!(mut.get_mut_one_hot() & mut.get_par_one_hot()));
        assert(__builtin_popcount(mut.get_mut_one_hot()) == 1);
        shared_mutations.push_back(mut);
    }
    void both(const MAT::Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        // Different need split
        if (sample_mut.get_mut_one_hot() != target_mut.get_mut_one_hot()) {
            if ((sample_mut.get_mut_one_hot() & target_mut.get_mut_one_hot())) {
                insert_split(sample_mut, target_mut, sample_mutations,
                             splitted_mutations, shared_mutations);
            } else {
                parsimony_score++;
                splitted_mutations.push_back(target_mut);
                sample_mutations.push_back(target_mut);
                sample_mutations.back().set_par_one_hot(
                    target_mut.get_par_one_hot());
            }
        } else {
            raise(SIGTRAP);
        }
    }
};

struct Upward_Descendant_Hook {
    MAT::Mutations_Collection &muts;
    Upward_Descendant_Hook(MAT::Mutations_Collection &muts) : muts(muts) {
        muts.clear();
    }
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        muts.reserve(sample_mutation_count + target_mutation_count);
    }
    void sample_mut_only(const MAT::Mutation &mut) {
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        // same as downward, pass through
        muts.push_back(mut);
    }
    void target_node_only(const MAT::Mutation &mut) {
        // shared mutation
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        muts.push_back(mut);
    }
    void both(const MAT::Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        if (target_mut.get_par_one_hot() != sample_mut.get_mut_one_hot()) {
            muts.push_back(sample_mut);
            muts.back().set_par_one_hot(target_mut.get_par_one_hot());
        }
    }
};

template <typename Hook>
static void generic_merge(const MAT::Node *node,
                          const MAT::Mutations_Collection &par_mutations,
                          Hook hook) {
    hook.reserve(par_mutations.size(), node->mutations.size());
    auto par_iter = par_mutations.begin();
    for (const auto &mut : node->mutations) {
        while (par_iter->get_position() < mut.get_position()) {
            hook.sample_mut_only(*par_iter);
            par_iter++;
        }
        if (par_iter->get_position() == mut.get_position()) {
            hook.both(*par_iter, mut);
            par_iter++;
        } else {
            hook.target_node_only(mut);
        }
    }
    while (par_iter < par_mutations.end()) {
        hook.sample_mut_only(*par_iter);
        par_iter++;
    }
}

struct Main_Tree_Searcher : public tbb::task {
    int radius_left;
    MAT::Mutations_Collection this_muts;
    const MAT::Node *node;
    const MAT::Node *exclude_node;
    Output<Main_Tree_Target> &output;
#ifndef NDEBUG
    Mutation_Set &sample_mutations;
#endif
    Main_Tree_Searcher(int radius_left, MAT::Node *node,
                       const MAT::Node *exclude_node,
                       Output<Main_Tree_Target> &output
#ifndef NDEBUG
                       ,
                       Mutation_Set &sample_mutations
#endif
                       )
        : radius_left(radius_left), node(node), exclude_node(exclude_node),
          output(output)
#ifndef NDEBUG
          ,
          sample_mutations(sample_mutations)
#endif
    {
    }
    void register_target(Main_Tree_Target &target, int this_score) {
        if (output.best_par_score >= this_score) {
            target.sample_mutations.mutations.pop_back();
            std::lock_guard<std::mutex> lk(output.mutex);
            if (output.best_par_score > this_score) {
                output.best_par_score = this_score;
                output.targets.clear();
            }
            if (output.best_par_score == this_score) {
                output.targets.push_back(std::move(target));
            }
        }
    }
    tbb::task *execute() {
        std::vector<Main_Tree_Searcher *> children_tasks;
        children_tasks.reserve(node->children.size() + 1);
        Main_Tree_Target target;
        auto cont = new (allocate_continuation()) tbb::empty_task();
        for (const auto child : node->children) {
            if (child == exclude_node) {
                continue;
            }
            target.target_node = child;
            target.parent_node = const_cast<MAT::Node *>(node);
            int child_radius_left = radius_left - child->mutations.size();
            int parsimony_score = 0;
            if (child->is_leaf() || (child_radius_left <= 0)) {
                generic_merge(child, this_muts,
                              Combine_Hook<Empty_Hook, Down_Sibling_Hook>{
                                  Empty_Hook(),
                                  Down_Sibling_Hook(target, parsimony_score)});
            } else {
                children_tasks.push_back(
                    new (cont->allocate_child()) Main_Tree_Searcher(
                        child_radius_left, child, node, output
#ifndef NDEBUG
                        ,
                        sample_mutations
#endif
                        ));
                generic_merge(
                    child, this_muts,
                    Combine_Hook<Down_Decendant_Hook, Down_Sibling_Hook>{
                        Down_Decendant_Hook(children_tasks.back()->this_muts),
                        Down_Sibling_Hook(target, parsimony_score)});
#ifndef NDEBUG
                check_continuation(children_tasks.back()->node,
                                   sample_mutations,
                                   children_tasks.back()->this_muts);
#endif
            }
#ifndef NDEBUG
            check_mutations(sample_mutations, target);
#endif
            register_target(target, parsimony_score);
        }
        auto parent = node->parent;
        int parsimony_score = 0;
        if (parent->dfs_index==54964) {
            //std::raise(SIGTRAP);
        }
        int radius_left_parent = radius_left - node->mutations.size();
        if (parent && parent != exclude_node) {
            target.target_node = const_cast<MAT::Node *>(node);
            target.parent_node = parent;
            if (radius_left_parent > 0) {
                children_tasks.push_back(
                    new (cont->allocate_child()) Main_Tree_Searcher(
                        radius_left_parent, parent, node, output
#ifndef NDEBUG
                        ,
                        sample_mutations
#endif
                        ));
                generic_merge(
                    node, this_muts,
                    Combine_Hook<Upward_Descendant_Hook, Upward_Sibling_Hook>{
                        Upward_Descendant_Hook(
                            children_tasks.back()->this_muts),
                        Upward_Sibling_Hook(target, parsimony_score)});
#ifndef NDEBUG
                check_continuation(children_tasks.back()->node,
                                   sample_mutations,
                                   children_tasks.back()->this_muts);
#endif
            } else {
                generic_merge(node, this_muts,
                              Combine_Hook<Empty_Hook, Upward_Sibling_Hook>{
                                  Empty_Hook(), Upward_Sibling_Hook(
                                                    target, parsimony_score)});
            }
            register_target(target, parsimony_score);
        }
        cont->set_ref_count(children_tasks.size());
        for (auto child : children_tasks) {
            cont->spawn(*child);
        }
        return children_tasks.empty() ? cont : nullptr;
    }
};
void set_parent_muts(std::vector<Sampled_Tree_Mutation> &mutations_to_set,
                     const MAT::Node *node) {
    std::unordered_map<int, uint8_t *> to_set;
    to_set.reserve(mutations_to_set.size());
    assert(mutations_to_set.back().position==INT_MAX);
    for (size_t idx=0; idx<(mutations_to_set.size()-1); idx++) {
        to_set.emplace(mutations_to_set[idx].position, &mutations_to_set[idx].par_nuc);
    }
    while (!to_set.empty() && node) {
        for (const auto &mut : node->mutations) {
            auto iter = to_set.find(mut.get_position());
            if (iter != to_set.end()) {
                *iter->second = mut.get_mut_one_hot();
                to_set.erase(iter);
            }
        }
        node = node->parent;
    }
    for (const auto &remaining : to_set) {
        *(remaining.second) = MAT::Mutation::refs[remaining.first];
    }
}
struct Main_Tree_Searcher_Par_Op {
    std::vector<Output<Main_Tree_Target>> &output;
    std::vector<Sampled_Place_Target> &sampled_tree_result;
    MAT::Tree &main_tree;
    int sampling_radius;
#ifndef NDEBUG
    Mutation_Set &sample_mutations;
#endif
    void operator()(tbb::blocked_range<size_t> &range) const {
        for (size_t idx = range.begin(); idx < range.end(); idx++) {
            auto sampled_root = sampled_tree_result[idx].target_node;
            auto &mutations_left_in = sampled_tree_result[idx].muts;
            auto main_tree_node = const_cast<MAT::Node *>(
                sampled_root->corresponding_main_tree_node);
            set_parent_muts(mutations_left_in, main_tree_node);
            // Mutation vector set as if descendent of main tree node
            MAT::Mutations_Collection mutations_left;
            mutations_left.reserve(mutations_left_in.size());
            for (const auto &mut : mutations_left_in) {
                mutations_left.push_back(MAT::Mutation(
                    mut.chrom_idx, mut.position, mut.par_nuc, mut.mut_nuc));
            }

            output[idx].best_par_score = INT_MAX;
            auto main_tree_task_root = new (tbb::task::allocate_root())
                Main_Tree_Searcher(sampling_radius, main_tree_node, nullptr,
                                   output[idx]
#ifndef NDEBUG
                                   ,
                                   sample_mutations
#endif
                );
            main_tree_task_root->this_muts = std::move(mutations_left);
            tbb::task::spawn_root_and_wait(*main_tree_task_root);
            assert(!output[idx].targets.empty());

        }
    }
};
std::tuple<Main_Tree_Target, int, int>
place_main_tree(std::vector<Sampled_Place_Target> &sampled_output,
                MAT::Tree &main_tree, int sampling_radius
#ifndef NDEBUG
                ,
                Mutation_Set &sample_mutations
#endif
) {
    std::vector<Output<Main_Tree_Target>> main_tree_out(sampled_output.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sampled_output.size()),
                      Main_Tree_Searcher_Par_Op{main_tree_out, sampled_output,
                                                main_tree, sampling_radius
#ifndef NDEBUG
                                                ,
                                                sample_mutations
#endif
                      });
    int min_score = main_tree_out[0].best_par_score;
    const auto *selected_target = &main_tree_out[0].targets[0];
    int best_idx = 0;
    for (size_t idx = 1; idx < sampled_output.size(); idx++) {
        if (main_tree_out[idx].best_par_score < min_score) {
            min_score = main_tree_out[idx].best_par_score;
            selected_target = &main_tree_out[idx].targets[0];
            best_idx = idx;
        }
    }
    return std::make_tuple(std::move(*selected_target), best_idx, min_score);
}
typedef std::unordered_map<const MAT::Node *, Sampled_Tree_Node *> node_map_t;
void map_sampled_nodes(
    Sampled_Tree_Node *sample_tree_root,
    std::unordered_map<const MAT::Node *, Sampled_Tree_Node *> &out) {
    out.emplace(sample_tree_root->corresponding_main_tree_node,
                sample_tree_root);
    for (auto child : sample_tree_root->children) {
        map_sampled_nodes(child, out);
    }
}
Sampled_Tree_Node *probe_neighbor(MAT::Node *node, MAT::Node *exclude_node,
                                  int radius, const node_map_t &node_map) {
    auto iter = node_map.find(node);
    if (iter != node_map.end()) {
        return iter->second;
    }
    auto parent = node->parent;
    if (parent && parent != exclude_node) {
        int par_rad = radius - parent->mutations.size();
        if (par_rad > 0) {
            auto ret = probe_neighbor(parent, node, par_rad, node_map);
            if (ret) {
                return ret;
            }
        }
    }
    for (auto child : node->children) {
        if (child == exclude_node) {
            continue;
        }
        int child_rad = radius - child->mutations.size();
        if (child_rad > 0) {
            auto ret = probe_neighbor(child, node, child_rad, node_map);
            if (ret) {
                return ret;
            }
        }
    }
    return nullptr;
}
int distance(const MAT::Node *node1, const MAT::Node *node2) {
    std::unordered_map<const MAT::Node *, int> node1_to_root_dists;
    auto node1_par = node1;
    int dist_accumulated = 0;
    while (node1_par) {
        node1_to_root_dists.emplace(node1_par, dist_accumulated);
        dist_accumulated += node1_par->mutations.size();
        node1_par = node1_par->parent;
    }
    int dist = 0;
    auto node2_par = node2;
    while (node2_par) {
        auto iter = node1_to_root_dists.find(node2_par);
        if (iter != node1_to_root_dists.end()) {
            return dist + iter->second;
        }
        dist += node2_par->mutations.size();
        node2_par = node2_par->parent;
    }
    raise(SIGTRAP);
    return INT_MAX;
}
#ifndef NDEBUG
void optimality_check(Mutation_Set &sample_mutations, int parsimony,
                      MAT::Node *main_tree_root, int sampling_radius,
                      Sampled_Tree_Node *sample_tree_root,
                      const std::vector<Sampled_Place_Target> &sampled_out) {
    Output<Main_Tree_Target> temp;
    MAT::Mutations_Collection ori_mutations;
    ori_mutations.mutations = std::vector<MAT::Mutation>(
        sample_mutations.begin(), sample_mutations.end());
    ori_mutations.reserve(sample_mutations.size());
    for (auto &mut : ori_mutations) {
        mut.set_par_one_hot(mut.get_ref_one_hot());
    }
    std::sort(ori_mutations.begin(), ori_mutations.end());
    ori_mutations.push_back(MAT::Mutation(0,INT_MAX,0x1,0xf));
    temp.best_par_score = INT_MAX;
    auto main_tree_task_root = new (tbb::task::allocate_root())
        Main_Tree_Searcher(sampling_radius, main_tree_root, nullptr, temp,
                           sample_mutations);
    main_tree_task_root->this_muts = std::move(ori_mutations);
    tbb::task::spawn_root_and_wait(*main_tree_task_root);
    if (parsimony != temp.best_par_score) {
        int min_dist = INT_MAX;
        auto closest_node_found = temp.targets[0].target_node;
        fprintf(stderr, "Suboptimal placement\n");
        for (const auto &sampled_target : sampled_out) {
            auto dist = distance(
                closest_node_found,
                sampled_target.target_node->corresponding_main_tree_node);
            if (dist < sampling_radius) {
                raise(SIGTRAP);
            }
        }
        if (min_dist > sampling_radius) {
            std::unordered_map<const MAT::Node *, Sampled_Tree_Node *> node_map;
            map_sampled_nodes(sample_tree_root, node_map);
            auto neighbor = probe_neighbor(closest_node_found, nullptr,
                                           sampling_radius, node_map);
            fprintf(stderr, "%s within radius but not found\n",
                    neighbor->corresponding_main_tree_node->identifier.c_str());
        }
        raise(SIGTRAP);
    }
}
#endif