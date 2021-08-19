#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "priority_conflict_resolver.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task.h>
#include <thread>
#include <unistd.h>
#include <utility>
#include "Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
size_t find_profitable_moves(MAT::Node *src, output_t &out,int radius,
                           stack_allocator<Mutation_Count_Change>& allocator,int starting_parsimony_score
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,MAT::Tree* tree
#endif
                          );
thread_local allocator_state<Mutation_Count_Change> FIFO_allocator_state;
extern tbb::task_group_context search_context;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst) {
    while (src!=dst) {
        //as dfs index of parent node will always smaller than its children's , so
        //substitute the node with larger dfs index with its parent and see whether
        //it will be the same as the other node
        if (src->dfs_index>dst->dfs_index) {
            src=src->parent;
        } else if (src->dfs_index<dst->dfs_index) {
            dst=dst->parent;
        }
    }
    return src;
}

bool check_not_ancestor(MAT::Node* dst,MAT::Node* src) {
    while (dst) {
        if (dst==src) {
            return false;
        }
        dst=dst->parent;
    }
    return true;
}
//The progress bar that does a linear extrapolation of time left for a round
static void print_progress(
    const std::atomic<int> *checked_nodes,
    std::chrono::time_point<
    std::chrono::steady_clock,
    std::chrono::duration<long, struct std::ratio<1, 1000000000>>>
    start_time,
    size_t total_nodes,
    const tbb::concurrent_vector<MAT::Node *> *deferred_nodes,
    bool *done,std::mutex* done_mutex) {
    while (true) {
        {
            if (*done) {
                return;
            }
            std::unique_lock<std::mutex> done_lock(*done_mutex);
            //I want it to print every second, but it somehow managed to print every minute...
            if (timed_print_progress) {
                progress_bar_cv.wait_for(done_lock,std::chrono::seconds(1));
            } else {
                progress_bar_cv.wait(done_lock);
            }
        }
        if (*done) {
            return;
        }
        int checked_nodes_temp = checked_nodes->load(std::memory_order_relaxed);
        std::chrono::duration<double> elpased_time =
            std::chrono::steady_clock::now() - start_time;
        double seconds_left = elpased_time.count() *
                              (total_nodes - checked_nodes_temp) /
                              checked_nodes_temp;
        if(((deferred_nodes->size())&&(std::chrono::steady_clock::now()-last_save_time)>=save_period)||deferred_nodes->size()>=max_queued_moves) {
            fprintf(stderr, "Timeout\n");
            return;
        }
        fprintf(stderr,"\rchecked %d nodes, estimated %f minutes left,found %zu nodes "
                "profitable",
                checked_nodes_temp, seconds_left / 60, deferred_nodes->size());
    }
}
void log_move_detail(const std::vector<Profitable_Moves_ptr_t> & moves, FILE* out,int iteration,int radius){
    for(const auto move:moves){
        fprintf(out, "%s\t%s\t%d\t%d\t%d\t%lu\n",move->src->identifier.c_str(),move->get_dst()->identifier.c_str()
            ,iteration,move->score_change,radius-move->radius_left,move->src->dfs_end_index-move->src->dfs_index);
    }
}
std::pair<size_t, size_t> optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
                     tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
                     MAT::Tree &t,int radius,FILE* log,bool allow_drift,int iteration
#ifndef NDEBUG
                     , Original_State_t origin_states
#endif
                    ) {
    std::atomic<size_t> node_searched_this_iter(0);
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", bfs_ordered_nodes.size());
    //for resolving conflicting moves
    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    Deferred_Move_t deferred_moves;
    Conflict_Resolver resolver(bfs_ordered_nodes.size(),deferred_moves);
    //progress bar
    std::atomic<int> checked_nodes(0);
    bool done=false;
    std::mutex done_mutex;
    auto search_start= std::chrono::steady_clock::now();
    std::thread progress_meter(print_progress,&checked_nodes,search_start, nodes_to_search.size(), &deferred_nodes,&done,&done_mutex);
    fputs("Start searching for profitable moves\n",stderr);
    //Actual search of profitable moves
    output_t out;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                      [&nodes_to_search, &resolver,
                                         &deferred_nodes,radius,&checked_nodes,&allow_drift,&node_searched_this_iter
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                         ,&t
#endif
                      ](tbb::blocked_range<size_t> r) {
        
        stack_allocator<Mutation_Count_Change> this_thread_FIFO_allocator(FIFO_allocator_state);
        size_t nodes_searched_this_batch=0;
        for (size_t i = r.begin(); i < r.end(); i++) {
            if (search_context.is_group_execution_cancelled()) {
                break;
            }
            if(((!deferred_nodes.size())||(std::chrono::steady_clock::now()-last_save_time)<save_period)&&deferred_nodes.size()<max_queued_moves) {
                output_t out;
                nodes_searched_this_batch+=find_profitable_moves(nodes_to_search[i], out, radius,this_thread_FIFO_allocator,allow_drift?-1:0
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                      ,&t
#endif
                                     );
                assert(this_thread_FIFO_allocator.empty());
                if (!out.moves.empty()) {
                    //resolve conflicts
                    if (allow_drift) {
                        std::shuffle(out.moves.begin(), out.moves.end(), rng);
                    } else {
                        deferred_nodes.push_back(
                            out.moves[0]->get_src());
                    }
                    resolver(out.moves);
                }
                checked_nodes.fetch_add(1,std::memory_order_relaxed);
            } else {
                deferred_nodes.push_back(nodes_to_search[i]);
            }
        }
        node_searched_this_iter.fetch_add(nodes_searched_this_batch,std::memory_order_relaxed);
    },search_context);
    {
        //stop the progress bar
        //std::unique_lock<std::mutex> done_lock(done_mutex);
        done=true;
        //done_lock.unlock();
        progress_bar_cv.notify_all();
    }
    auto searh_end=std::chrono::steady_clock::now();
    std::chrono::duration<double> elpased_time =searh_end-search_start;
    fprintf(stderr, "\nSearch took %f minutes\n",elpased_time.count()/60);
    //apply moves
    fputs("Start applying moves\n",stderr);
    std::vector<Profitable_Moves_ptr_t> all_moves;
    resolver.schedule_moves(all_moves);
    fprintf(stderr, "Applying %zu moves\n",all_moves.size());
    if (iteration>0) {
        log_move_detail(all_moves, log, iteration, radius);
    }
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
                ,origin_states
#endif
               );
    auto apply_end=std::chrono::steady_clock::now();
    elpased_time =apply_end-searh_end;
    fprintf(stderr, "apply moves took %f seconds\n",elpased_time.count());
    //recycle conflicting moves
    int init_deferred=deferred_moves.size();
    int recycled=0;
    fputs("Start recycling conflicting moves\n",stderr);
    while (!deferred_moves.empty()&&(!allow_drift)) {
        bfs_ordered_nodes=t.breadth_first_expansion();
        {
            Deferred_Move_t deferred_moves_next;
            Conflict_Resolver resolver(bfs_ordered_nodes.size(),deferred_moves_next);
            if (timed_print_progress) {
                fprintf(stderr,"\rrecycling conflicting moves, %zu left",deferred_moves.size());
            }
            tbb::parallel_for(tbb::blocked_range<size_t>(0,deferred_moves.size()),[&deferred_moves,&resolver,&t](const tbb::blocked_range<size_t>& r) {
                for (size_t i=r.begin(); i<r.end(); i++) {
                    MAT::Node* src=t.get_node(deferred_moves[i].first);
                    MAT::Node* dst=t.get_node(deferred_moves[i].second);
                    if (src&&dst) {
                        output_t out;
                        if (check_not_ancestor(dst, src)) {
                            individual_move(src,dst,get_LCA(src, dst),out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                            ,&t
#endif
                                           );
                        }
                        if (!out.moves.empty()) {
                            resolver(out.moves);
                        }
                    }

                }
            });
            all_moves.clear();
            resolver.schedule_moves(all_moves);
            if (iteration>0) {
                log_move_detail(all_moves, log, iteration, radius);
            }
            recycled+=all_moves.size();
            apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
                        ,origin_states
#endif
                       );
            deferred_moves=std::move(deferred_moves_next);
        }

    }
    fprintf(stderr, "Recycled %d moves\n",recycled);
    auto recycle_end=std::chrono::steady_clock::now();
    elpased_time =recycle_end-apply_end;
#ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
#endif
    nodes_to_search = std::move(deferred_nodes);
    progress_meter.join();
    fprintf(stderr, "recycled %f of conflicting moves \n",(double)recycled/(double)init_deferred);
    fprintf(stderr, "recycling moves took %f seconds\n",elpased_time.count());
    return std::make_pair(t.get_parsimony_score(),node_searched_this_iter.load());
}
tbb::concurrent_unordered_map<MAT::Mutation,
    tbb::concurrent_unordered_map<std::string, nuc_one_hot> *,
    Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
void save_final_tree(MAT::Tree &t, Original_State_t& origin_states,
                     const std::string &output_path) {
#ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
#endif
    std::vector<MAT::Node *> dfs = t.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
    [&dfs](tbb::blocked_range<size_t> r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            dfs[i]->mutations.remove_invalid();
        }
    });
    fix_condensed_nodes(&t);
    fprintf(stderr, "%zu condensed_nodes\n",t.condensed_nodes.size());
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, output_path);
}
