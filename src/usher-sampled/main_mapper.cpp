#include "mapper.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
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
    void sample_mut_only(const To_Place_Sample_Mutation &mut) {
        hook1.sample_mut_only(mut);
        hook2.sample_mut_only(mut);
    }
    void target_node_only(const MAT::Mutation &mut) {
        hook1.target_node_only(mut);
        hook2.target_node_only(mut);
    }
    void target_N_skiped(const MAT::Mutation &mut) {
        hook1.target_N_skiped(mut);
        hook2.target_N_skiped(mut);
    }
    void both(const To_Place_Sample_Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        hook1.both(sample_mut, target_mut);
        hook2.both(sample_mut, target_mut);
    }
};
struct Empty_Hook {
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {}
    void sample_mut_only(const To_Place_Sample_Mutation &mut) {}
    void target_node_only(const MAT::Mutation &mut) {}
    void target_N_skiped(const MAT::Mutation &mut) {}
    void both(const To_Place_Sample_Mutation &mut, const MAT::Mutation &mut2) {}
};
static void insert_split(const To_Place_Sample_Mutation &sample_mut,
                         const MAT::Mutation &target_mut,
                         std::vector<To_Place_Sample_Mutation> &sample_mutations,
                         MAT::Mutations_Collection &splitted_mutations,
                         MAT::Mutations_Collection &shared_mutations) {
    auto par_nuc = target_mut.get_par_one_hot();
    auto sample_nuc = sample_mut.mut_nuc;
    auto target_nuc = target_mut.get_mut_one_hot();
    auto coincided_mut = sample_nuc & target_nuc;
    if (coincided_mut & par_nuc) {
        coincided_mut = coincided_mut & par_nuc;
    }
    coincided_mut = 1 << __builtin_ctz(coincided_mut);
    assert(coincided_mut & sample_nuc & target_nuc);
    if (par_nuc!=coincided_mut) {
        shared_mutations.push_back(target_mut);
        shared_mutations.back().set_par_mut(par_nuc, coincided_mut);
        assert(shared_mutations.back().get_descendant_mut()&shared_mutations.back().get_mut_one_hot());
    }
    if (coincided_mut != target_nuc) {
        splitted_mutations.push_back(target_mut);
        splitted_mutations.back().set_par_mut(coincided_mut, target_nuc);
    }
    if (coincided_mut != sample_nuc) {
        sample_mutations.push_back(sample_mut);
        sample_mutations.back().par_nuc=coincided_mut;
        sample_mutations.back().mut_nuc= sample_nuc;
    }
}
static void n_skiped_sibling(
        MAT::Mutations_Collection &splitted_mutations,
    MAT::Mutations_Collection &shared_mutations,
    const MAT::Mutation &mut
){
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        if(__builtin_popcount(mut.get_mut_one_hot())==1){
            splitted_mutations.push_back(mut);
            return;
        }
        auto concensus_mut=mut.get_mut_one_hot()&mut.get_par_one_hot();
        if (!concensus_mut) {
            concensus_mut=1<<__builtin_ctz(mut.get_mut_one_hot());
            assert(shared_mutations.empty()||shared_mutations.back().get_position()<mut.get_position());
            shared_mutations.push_back(MAT::Mutation(mut.get_chromIdx(),mut.get_position(),mut.get_par_one_hot(),concensus_mut));
            shared_mutations.back().set_descendant_mut(mut.get_descendant_mut());
            assert(shared_mutations.back().get_descendant_mut()&shared_mutations.back().get_mut_one_hot());
        }
        assert(splitted_mutations.empty()||splitted_mutations.back().get_position()<mut.get_position());
        splitted_mutations.push_back(mut);
        splitted_mutations.back().set_par_one_hot(concensus_mut);
}
static void sample_check_mutation(
    const std::vector<To_Place_Sample_Mutation> &sample_mutations,
    const MAT::Mutations_Collection &splitted_mutations,
    const MAT::Mutations_Collection &shared_mutations,
    int position
){
#ifdef DETAILED_MERGER_CHECK
        assert(shared_mutations.empty()||shared_mutations.mutations.back().get_position()<position);
        assert(splitted_mutations.empty()||splitted_mutations.mutations.back().get_position()<position);
        assert(sample_mutations.empty()||sample_mutations.back().position<position);
#endif
}
static void sample_mut_check_mutation(
    const std::vector<To_Place_Sample_Mutation> &sample_mutations,
    const MAT::Mutations_Collection &splitted_mutations,
    const MAT::Mutations_Collection &shared_mutations,
    const To_Place_Sample_Mutation& mut
){
#ifdef DETAILED_MERGER_CHECK
        assert(shared_mutations.empty()||shared_mutations.mutations.back().get_position()<=mut.get_end_range());
        assert(splitted_mutations.empty()||splitted_mutations.mutations.back().get_position()<=mut.get_end_range());
        assert(sample_mutations.empty()||sample_mutations.back().position<mut.position);
#endif
}
struct Down_Sibling_Hook {
    std::vector<To_Place_Sample_Mutation> &sample_mutations;
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
        sample_mutations.reserve(sample_mutation_count+1);
        splitted_mutations.reserve(target_mutation_count+1);
        shared_mutations.reserve(target_mutation_count+1);
    }
    void sample_mut_only(const To_Place_Sample_Mutation &mut) {
        sample_mut_check_mutation(sample_mutations, splitted_mutations, shared_mutations, mut);
        sample_mutations.push_back(mut);
        if (mut.mut_nuc!=0xf) {
            assert(mut.mut_nuc != mut.par_nuc);
            if (!(mut.par_nuc & mut.mut_nuc)) {
                parsimony_score++;
            }
        }
    }
    void target_node_only(const MAT::Mutation &mut) {
        sample_check_mutation(sample_mutations, splitted_mutations, shared_mutations, mut.get_position());
        assert(mut.get_mut_one_hot() != mut.get_par_one_hot());
        splitted_mutations.push_back(mut);
    }
    void target_N_skiped(const MAT::Mutation &mut) {
        sample_check_mutation(sample_mutations, splitted_mutations, shared_mutations, mut.get_position());
        n_skiped_sibling(splitted_mutations, shared_mutations, mut);
    }
    void both(const To_Place_Sample_Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        /*if (sample_mut.get_position()==241) {
            fputc('a', stderr);
        }*/
        assert(sample_mut.position!=0xf);
        assert(sample_mut.position==target_mut.get_position());
        sample_check_mutation(sample_mutations, splitted_mutations, shared_mutations, sample_mut.position);
        if (sample_mut.mut_nuc != target_mut.get_mut_one_hot()) {
            if (sample_mut.mut_nuc & target_mut.get_mut_one_hot()) {
                if (target_mut.get_par_one_hot()&target_mut.get_mut_one_hot()) {
                    splitted_mutations.push_back(target_mut);
                    sample_mutations.push_back(sample_mut);
                    if (!(sample_mut.mut_nuc&sample_mut.par_nuc)) {
                        parsimony_score++;
                    }
                }
                else {
                    insert_split(sample_mut, target_mut, sample_mutations,
                             splitted_mutations, shared_mutations);
                }
            } else {
                sample_mutations.push_back(sample_mut);
                splitted_mutations.push_back(target_mut);
                    if (!(sample_mut.mut_nuc&sample_mut.par_nuc)) {
                        parsimony_score++;
                    }
            }
        } else {
            if (__builtin_popcount(target_mut.get_mut_one_hot())!=1) {
                auto common=target_mut.get_par_one_hot()&target_mut.get_mut_one_hot();
                if (!common) {
                    common=1 << __builtin_ctz(target_mut.get_mut_one_hot());
                    shared_mutations.push_back(target_mut);
                    shared_mutations.back().set_mut_one_hot(common);
                    assert(shared_mutations.back().get_descendant_mut()&shared_mutations.back().get_mut_one_hot());
                }
                sample_mutations.push_back(sample_mut);
                sample_mutations.back().par_nuc=common;
                if (!(sample_mut.mut_nuc&common)) {
                        parsimony_score++;
                }
                splitted_mutations.push_back(target_mut);
                splitted_mutations.back().set_par_one_hot(common);
            }else {
                shared_mutations.push_back(target_mut);
                assert(shared_mutations.back().get_descendant_mut()&shared_mutations.back().get_mut_one_hot());
            }
        }
        assert(shared_mutations.empty()||!(shared_mutations.back().get_par_one_hot()&shared_mutations.back().get_mut_one_hot()));
    }
};
struct Down_Decendant_Hook {
    std::vector<To_Place_Sample_Mutation> &muts;
    int& lower_bound;
    Down_Decendant_Hook(std::vector<To_Place_Sample_Mutation> &muts,int& lower_bound) : muts(muts),lower_bound(lower_bound) {
        lower_bound=0;
        muts.clear();
    }
    void reserve(size_t sample_mutation_count, size_t target_mutation_count) {
        muts.reserve(sample_mutation_count + target_mutation_count+1);
    }
    void sample_mut_only(const To_Place_Sample_Mutation &mut) {
        assert(mut.mut_nuc==0xf||mut.mut_nuc != mut.par_nuc);
        assert(muts.empty()||muts.back().position<mut.position);
        muts.push_back(mut);
        if (mut.mut_nuc!=0xf) {
        if (!(mut.descendent_possible_nuc&mut.mut_nuc)) {
            assert(!(mut.par_nuc&mut.mut_nuc));
            lower_bound++;
        }            
        }
    }
    void target_N_skiped(const MAT::Mutation &mut) {
    }
    void target_node_only(const MAT::Mutation &mut) {
        assert(!(mut.get_mut_one_hot() & mut.get_par_one_hot()));
        assert(muts.empty()||muts.back().position<mut.get_position());
        muts.push_back(To_Place_Sample_Mutation(mut.get_position(),mut.get_chromIdx(),mut.get_par_one_hot(),mut.get_mut_one_hot(),mut.get_descendant_mut()));
        if (!(muts.back().mut_nuc&muts.back().descendent_possible_nuc)) {
            assert(!(muts.back().mut_nuc&muts.back().par_nuc));
            lower_bound++;
        }
    }
    void both(const To_Place_Sample_Mutation &sample_mut,
              const MAT::Mutation &target_mut) {
        /*if (sample_mut.get_position()==241) {
            fputc('a', stderr);
        }*/
        assert(sample_mut.mut_nuc!=0xf);
        assert(sample_mut.position==target_mut.get_position());
        assert(muts.empty()||muts.back().position<sample_mut.position);
        if (sample_mut.mut_nuc != target_mut.get_mut_one_hot()) {
            muts.push_back(sample_mut);
            auto &last_mut = muts.back();
            last_mut.par_nuc=target_mut.get_mut_one_hot();
            last_mut.descendent_possible_nuc&=target_mut.get_descendant_mut();
            if (!(last_mut.descendent_possible_nuc&last_mut.mut_nuc)) {
                assert(!(last_mut.par_nuc&last_mut.mut_nuc));
                lower_bound++;
            }
        }
        //assert(muts.empty()||muts.back().get_mut_one_hot() != muts.back().get_par_one_hot());
    }
};

template <typename Hook>
static void generic_merge(const MAT::Node *node,
                          const std::vector<To_Place_Sample_Mutation> &par_mutations,
                          Hook hook) {
    hook.reserve(par_mutations.size(), node->mutations.size());
    auto par_iter = par_mutations.begin();
    for (const auto &mut : node->mutations) {
        while (par_iter->get_end_range() < mut.get_position()) {
            hook.sample_mut_only(*par_iter);
            par_iter++;
        }
        if (par_iter->mut_nuc==0xf&&par_iter->position<=mut.get_position()) {
            hook.target_N_skiped(mut);
            continue;
        }
        if (par_iter->position == mut.get_position()) {
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
    void register_target(Main_Tree_Target &target, int this_score,Output<Main_Tree_Target> &output) {
        #ifndef NDEBUG
        int initial_par_score=0;
        for (const auto & mut : target.target_node->mutations) {
            if (!(mut.get_par_one_hot()&mut.get_mut_one_hot())) {
                initial_par_score++;
            }
        }
        for (const auto& mut : target.splited_mutations) {
            if (!(mut.get_par_one_hot()&mut.get_mut_one_hot())) {
                initial_par_score--;
            }
        }
        initial_par_score-=target.shared_mutations.size();
        assert(initial_par_score==0);
        int sample_par_score=0;
        for (const auto& mut : target.sample_mutations) {
            if (mut.mut_nuc!=0xf&&!(mut.par_nuc&mut.mut_nuc)) {
                sample_par_score++;
            }
        }
        assert(sample_par_score==this_score);
        #endif
        if (output.best_par_score >= this_score) {
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
static void search_serial(const MAT::Node* node,std::vector<To_Place_Sample_Mutation>& this_muts,Output<Main_Tree_Target> &output){
    Main_Tree_Target target;
    for (const auto child : node->children) {
            target.target_node = child;
            target.parent_node = const_cast<MAT::Node *>(node);
            int parsimony_score = 0;
            if (child->is_leaf()) {
                generic_merge(child, this_muts,
                              Combine_Hook<Empty_Hook, Down_Sibling_Hook>{
                                  Empty_Hook(),
                                  Down_Sibling_Hook(target, parsimony_score)});
                register_target(target, parsimony_score,output);
            } else {
                int lower_bound = 0;
                std::vector<To_Place_Sample_Mutation> descendant_mutations;
                generic_merge(
                    child, this_muts,
                    Combine_Hook<Down_Decendant_Hook, Down_Sibling_Hook>{
                        Down_Decendant_Hook(descendant_mutations, lower_bound),
                        Down_Sibling_Hook(target, parsimony_score)});
#ifdef DETAILED_MERGER_CHECK
                check_continuation(child,
                                   sample_mutations, descendant_mutations);
#endif
                register_target(target, parsimony_score,output);
#ifndef BOUND_CHECK
                assert(parsimony_score>=curr_lower_bound);
                assert(curr_lower_bound<=lower_bound);
                if (lower_bound <= output.best_par_score) {
#endif
                    search_serial(child, descendant_mutations, output);
#ifndef BOUND_CHECK
                }
#endif
            }
#ifdef DETAILED_MERGER_CHECK
            check_mutations(sample_mutations, target);
#endif
        }
}
struct Main_Tree_Searcher : public tbb::task {
    int curr_lower_bound;
    std::vector<To_Place_Sample_Mutation> this_muts;
    const MAT::Node *node;
    Output<Main_Tree_Target> &output;
#ifdef DETAILED_MERGER_CHECK
    Mutation_Set &sample_mutations;
#endif
    Main_Tree_Searcher(int curr_lower_bound,MAT::Node *node,
                       Output<Main_Tree_Target> &output
#ifdef DETAILED_MERGER_CHECK
                       ,
                       Mutation_Set &sample_mutations
#endif
                       )
        : curr_lower_bound(curr_lower_bound),node(node),
          output(output)
#ifdef DETAILED_MERGER_CHECK
          ,
          sample_mutations(sample_mutations)
#endif
    {
    }
    tbb::task *execute() {
#ifndef BOUND_CHECK
        if(curr_lower_bound>output.best_par_score){
            return nullptr;
        }
#endif
        if(node->bfs_index<10){
            search_serial(node, this_muts, output);
            return nullptr;
        }
        std::vector<Main_Tree_Searcher *> children_tasks;
        children_tasks.reserve(node->children.size() + 1);
        Main_Tree_Target target;
        auto cont = new (allocate_continuation()) tbb::empty_task();
        for (const auto child : node->children) {
            target.target_node = child;
            target.parent_node = const_cast<MAT::Node *>(node);
            int parsimony_score = 0;
            if (child->is_leaf()) {
                generic_merge(child, this_muts,
                              Combine_Hook<Empty_Hook, Down_Sibling_Hook>{
                                  Empty_Hook(),
                                  Down_Sibling_Hook(target, parsimony_score)});
            } else {
                int lower_bound = 0;
                std::vector<To_Place_Sample_Mutation> descendant_mutations;
                generic_merge(
                    child, this_muts,
                    Combine_Hook<Down_Decendant_Hook, Down_Sibling_Hook>{
                        Down_Decendant_Hook(descendant_mutations, lower_bound),
                        Down_Sibling_Hook(target, parsimony_score)});
#ifdef DETAILED_MERGER_CHECK
                check_continuation(child,
                                   sample_mutations, descendant_mutations);
#endif
                assert(curr_lower_bound<=lower_bound);
#ifndef BOUND_CHECK
                if (lower_bound <= output.best_par_score) {
#endif
                    children_tasks.push_back(
                        new (cont->allocate_child())
                            Main_Tree_Searcher(lower_bound, child, output
#ifdef DETAILED_MERGER_CHECK
                                               ,
                                               sample_mutations
#endif
                                               ));
                    children_tasks.back()->this_muts =
                        std::move(descendant_mutations);
#ifndef BOUND_CHECK
                }
#endif
            }
#ifdef DETAILED_MERGER_CHECK
            check_mutations(sample_mutations, target);
#endif
            assert(parsimony_score>=curr_lower_bound);
            register_target(target, parsimony_score,output);
        }
        cont->set_ref_count(children_tasks.size());
        for (auto child : children_tasks) {
            cont->spawn(*child);
        }
        return children_tasks.empty() ? cont : nullptr;
    }
};
std::tuple<std::vector<Main_Tree_Target>, int>
place_main_tree(std::vector<To_Place_Sample_Mutation> &mutations,
                MAT::Tree &main_tree
#ifdef DETAILED_MERGER_CHECK
                ,
                Mutation_Set &sample_mutations
#endif
){
            Output<Main_Tree_Target> output;
            output.targets.reserve(1000);
            output.best_par_score = INT_MAX;
            auto main_tree_task_root = new (tbb::task::allocate_root())
                Main_Tree_Searcher(0,main_tree.root,
                                   output
#ifdef DETAILED_MERGER_CHECK
                                   ,
                                   sample_mutations
#endif
                );
            main_tree_task_root->this_muts = mutations;
            To_Place_Sample_Mutation temp(INT_MAX,0,0xf);
            main_tree_task_root->this_muts.push_back(temp);
            tbb::task::spawn_root_and_wait(*main_tree_task_root);
            assert(!output.targets.empty());

    return std::make_tuple(std::move(output.targets), output.best_par_score);
}
