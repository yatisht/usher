#include "mapper.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cassert>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <string>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#include <utility>
#include <vector>

struct Lower_Bound_Hook {
    int &lower_bound;
    int &mutation_count;
    std::vector<To_Place_Sample_Mutation> &output;
    Lower_Bound_Hook(int &lower_bound, int &mutation_count,
                     std::vector<To_Place_Sample_Mutation> &output)
        : lower_bound(lower_bound), mutation_count(mutation_count),
          output(output) {
        lower_bound = 0;
        mutation_count = 0;
        output.clear();
    }
    void add(To_Place_Sample_Mutation& mut, const Sampled_Tree_Node *node) {
        assert(output.empty()||output.back().get_end_range()<mut.position);
        if (mut.mut_nuc!=0xf&&(!(mut.par_nuc & mut.mut_nuc))) {
            if ((mut.mut_nuc & mut.descendent_possible_nuc) == 0) {
                lower_bound++;
            }
            mutation_count++;
        }
        output.push_back(mut);
    }
};
struct Add_Only_Hook {
    std::vector<To_Place_Sample_Mutation> &output;
    int &mutation_count;
    Add_Only_Hook(std::vector<To_Place_Sample_Mutation> &output,
                  int &mutation_count)
        : output(output), mutation_count(mutation_count) {}
    void add(To_Place_Sample_Mutation mut, const ::Sampled_Tree_Node *node) {
        assert(output.empty()||output.back().get_end_range()<mut.position);
        if (mut.mut_nuc!=0xf && ((mut.mut_nuc & mut.par_nuc) == 0)) {
            mutation_count++;
        }
        output.push_back(mut);
    }
};

template <typename Mut_Add_Hook>
static void merge_mutation(std::vector<To_Place_Sample_Mutation> &par_mutations,
                           Sampled_Tree_Node *node, Mut_Add_Hook &hook) {
    hook.output.reserve(par_mutations.size() + node->mutations.size());
    auto par_iter = par_mutations.begin();
    for (const auto &mut : node->mutations) {
        while (par_iter->get_end_range() < mut.position) {
            hook.add(*par_iter, node);
            par_iter++;
        }
        assert(mut.position<=par_iter->get_end_range());
        if (par_iter->mut_nuc==0xf&&mut.position>=par_iter->position) {
            continue;
        }
        assert(mut.mut_nuc!=0xf);
        //assert(par_iter->mut_nuc!=0xf);
        if (par_iter->position == mut.position) {
            if (par_iter->mut_nuc != mut.mut_nuc) {
                To_Place_Sample_Mutation temp(*par_iter);
                temp.par_nuc = mut.mut_nuc;
                temp.descendent_possible_nuc&=mut.descendent_possible_nuc;
                hook.add(temp, node);
            }
            par_iter++;
        } else {
            To_Place_Sample_Mutation temp(mut.position,mut.chrom_idx,mut.par_nuc,mut.mut_nuc);
            hook.add(temp, node);
        }
    }
    while (par_iter < par_mutations.end()) {
        hook.add(*par_iter, node);
        par_iter++;
    }
}

// TODO: Handle sampled tree branch spliting
struct Sampled_Tree_Placer_Task : public tbb::task {
    std::vector<To_Place_Sample_Mutation> this_muts;
    Output<Sampled_Place_Target> &output;
    const Sampled_Tree_Node *node;
    int lower_bound;
#ifndef NDEBUG
    Mutation_Set &sample_mutations;
#endif
    Sampled_Tree_Placer_Task(const Sampled_Tree_Node *node,
                             Output<Sampled_Place_Target> &output,
                             std::vector<To_Place_Sample_Mutation> &&mutations,
                             int lower_bound
#ifndef NDEBUG
                             ,
                             Mutation_Set &sample_mutations
#endif

                             )
        : this_muts(std::move(mutations)), output(output), node(node),
          lower_bound(lower_bound)
#ifndef NDEBUG
          ,
          sample_mutations(sample_mutations)
#endif
    {
    }
    void register_sample(const Sampled_Place_Target& target,
                         int par_score) {
#ifdef DETAILED_MERGER_CHECK
        check_sampled_mutations(sample_mutations, target);
#endif
        assert(target.muts.back().position==INT_MAX);
        if (par_score <= output.best_par_score) {
            std::lock_guard<std::mutex> lk(output.mutex);
            if (par_score < output.best_par_score) {
                output.best_par_score = par_score;
                output.targets.clear();
            }
            if (par_score == output.best_par_score) {
                output.targets.emplace_back(target);
            }
        }
    }
    tbb::task *execute() override {
        if (lower_bound > output.best_par_score) {
            return nullptr;
        }
        if (node->dfs_idx == 3026) {
            //std::raise(SIGTRAP);
        }
        std::vector<Sampled_Tree_Placer_Task *> children_tasks;
        children_tasks.reserve(node->children.size());
        auto cont = new (allocate_continuation()) tbb::empty_task();
        assert(this_muts.back().position==INT_MAX);
        for (const auto child : node->children) {
            Sampled_Place_Target target;
            auto& child_mutations=target.muts;
            target.target_node=child;
            int mutation_count = 0;
            if (child->children.empty()) {
                Add_Only_Hook hook(child_mutations, mutation_count);
                merge_mutation(this_muts, child, hook);
                register_sample(target, mutation_count);
            } else {
                int par_lower_bound = 0;
                Lower_Bound_Hook hook(par_lower_bound, mutation_count,
                                      child_mutations);
                merge_mutation(this_muts, child, hook);
                lower_bound = std::min(lower_bound, par_lower_bound);
                register_sample(target,mutation_count);
                if (lower_bound <= output.best_par_score) {
                    auto child_task = new (cont->allocate_child())
                        Sampled_Tree_Placer_Task(child, output,
                                                 std::move(child_mutations),
                                                 par_lower_bound
#ifndef NDEBUG
                                                 ,
                                                 sample_mutations
#endif
                        );
                    assert(child_task->this_muts.back().position==INT_MAX);
                    children_tasks.push_back(child_task);
                }
            }
        }
        cont->set_ref_count(children_tasks.size());
        if (lower_bound <= output.best_par_score) {
            for (auto task : children_tasks) {
                cont->spawn(*task);
            }
        } else {
            for (auto task : children_tasks) {
                cont->destroy(*task);
            }
            children_tasks.clear();
            cont->set_ref_count(0);
        }
        return children_tasks.empty() ? cont : nullptr;
    }
};

std::vector<Sampled_Place_Target>
place_on_sampled_tree(Sampled_Tree_Node *sampled_tree_root,
                      std::vector<To_Place_Sample_Mutation> &&sample_mutations,
                      int& parsimony_score
#ifndef NDEBUG
                      ,
                      Mutation_Set &sample_mutations_set
#endif
) {
    Output<Sampled_Place_Target> sampled_output;
    sampled_output.best_par_score = INT_MAX;
    To_Place_Sample_Mutation temp(INT_MAX,0,0xf);
    sample_mutations.push_back(temp);
    tbb::task::spawn_root_and_wait(
        *(new (tbb::task::allocate_root()) Sampled_Tree_Placer_Task(
            sampled_tree_root, sampled_output, std::move(sample_mutations), 0
#ifndef NDEBUG
            ,
            sample_mutations_set
#endif
            )));
    assert(!sampled_output.targets.empty());
    parsimony_score = sampled_output.best_par_score;
    return sampled_output.targets;
}