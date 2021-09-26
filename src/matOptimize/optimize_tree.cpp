#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "priority_conflict_resolver.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <memory>
#include <mpi.h>
#include <mutex>
#include <random>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task.h>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>
#include "Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
extern bool changing_radius;
#define MAX_MOVE_SIZE ((size_t)0x1000)
#define MAX_MOVE_MSG_SIZE (1+4*MAX_MOVE_SIZE)
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
        /*if(((deferred_nodes->size())&&(std::chrono::steady_clock::now()-last_save_time)>=save_period)||deferred_nodes->size()>=max_queued_moves) {
            fprintf(stderr, "Timeout\n");
            return;
        }*/
        fprintf(stderr,"\rchecked %d nodes, estimated %f minutes left,found %zu nodes "
                "profitable",
                checked_nodes_temp, seconds_left / 60, deferred_nodes->size());
    }
}
void log_move_detail(const std::vector<Profitable_Moves_ptr_t> & moves, FILE* out,int iteration,int radius) {
    for(const auto move:moves) {
        fprintf(out, "%s\t%s\t%d\t%d\t%d\t%lu\n",move->src->identifier.c_str(),move->dst->identifier.c_str()
                ,iteration,move->score_change,radius-move->radius_left,move->src->dfs_end_index-move->src->dfs_index);
    }
}

typedef tbb::flow::function_node<std::vector<Profitable_Moves_ptr_t>*> resolver_node_t;
static void MPI_recieve_move(const std::vector<MAT::Node*>& dfs_ordered_nodes,resolver_node_t& resover_node ){
    int* buffer=new int[MAX_MOVE_MSG_SIZE];
    while (true) {
    std::vector<Profitable_Moves_ptr_t>* out =new std::vector<Profitable_Moves_ptr_t>;
    MPI_Status stat;
    MPI_Recv(buffer,MAX_MOVE_MSG_SIZE , MPI_INT, MPI_ANY_SOURCE, MOVE_TAG, MPI_COMM_WORLD, &stat);
    auto msg_size=stat._ucount;
    if (msg_size==0) {
        fprintf(stderr, "+++++++++Move receiver exit\n");
        return;
    }
    int moves_size=(msg_size/4-1)/4;
    //fprintf(stderr, "Recieved %d moves\n",moves_size);
    MAT::Node* src=dfs_ordered_nodes[buffer[0]];
    out->reserve(moves_size);
    for (int move_idx=0;move_idx<moves_size;move_idx++){
        int LCA_idx=buffer[1+4*move_idx];
        int dst_idx=buffer[2+4*move_idx];
        int score_change=buffer[3+4*move_idx];
        int radius_left=buffer[4+4*move_idx];
        //fprintf(stderr, "Recevieved move with src %d, dst %d,LCA %d,score_change %d, radius %d \n",buffer[0],dst_idx,LCA_idx,score_change,radius_left);
        out->push_back(Profitable_Moves_ptr_t(new Profitable_Moves{score_change,src,dfs_ordered_nodes[dst_idx],dfs_ordered_nodes[LCA_idx],radius_left}));
    }
    resover_node.try_put(out);
    }
    delete[] buffer;
}
struct MPI_move_sender{
    void operator()(std::vector<Profitable_Moves_ptr_t>* to_send){
        size_t moves_to_send=std::min(to_send->size(),MAX_MOVE_SIZE);
        size_t msg_size=1+4*moves_to_send;
        int* buffer=new int[msg_size];
        buffer[0]=(*to_send)[0]->src->dfs_index;
        for (size_t move_idx=0; move_idx<moves_to_send; move_idx++) {
            const auto& this_move=(*(*to_send)[move_idx]);
            buffer[1+4*move_idx]=this_move.LCA->dfs_index;
            buffer[2+4*move_idx]=this_move.dst->dfs_index;
            buffer[3+4*move_idx]=this_move.score_change;
            buffer[4+4*move_idx]=this_move.radius_left;
            //fprintf(stderr, "sent move with src %d, dst %d,LCA %d,score_change %d, radius %d \n",buffer[0],buffer[2+4*move_idx],buffer[1+4*move_idx],buffer[3+4*move_idx],buffer[4+4*move_idx]);
        }
        //fprintf(stderr, "Sent %zu moves\n",moves_to_send);
        MPI_Send(buffer, msg_size, MPI_INT, 0, MOVE_TAG, MPI_COMM_WORLD);
        delete[] buffer;
        delete to_send;
    }
};
struct progress_meter_t{
    bool done;
    std::mutex done_mutex;
    std::thread progress_meter;
    progress_meter_t(size_t total_size, std::atomic<int>& checked_nodes,tbb::concurrent_vector<MAT::Node*>* deferred_nodes):done(false),progress_meter(print_progress,&checked_nodes,std::chrono::steady_clock::now(), total_size, deferred_nodes,&done,&done_mutex){
    }
    void stop(){
        //stop the progress bar
        std::unique_lock<std::mutex> done_lock(done_mutex);
        done=true;
        done_lock.unlock();
        progress_bar_cv.notify_all();
    }
    ~progress_meter_t(){
        progress_meter.join();
    }
};
template<typename T>
static void search(std::vector<MAT::Node *> &nodes_to_search,MAT::Tree &t,int radius,tbb::concurrent_vector<MAT::Node*>& deferred_nodes,T& resolver_node){
    std::atomic<int> checked_nodes(0);
    progress_meter_t progress_meter(nodes_to_search.size(),checked_nodes,&deferred_nodes);
    fputs("Start searching for profitable moves\n",stderr);
    //Actual search of profitable moves
    #ifdef CHECK_BOUND
std::atomic<size_t> saved(0);
std::atomic<size_t> total(0);
#endif
    t.populate_ignored_range();
/*
counters count;
output_t out;
find_moves_bounded(t.get_node("s14640s"), out,radius,count);
*/
    auto search_start=std::chrono::steady_clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                      [&nodes_to_search, &resolver_node,
                                         &deferred_nodes,radius,&checked_nodes
    #ifdef CHECK_BOUND
                                         ,&total,&saved
    #endif
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                         ,&t
#endif
                      ](tbb::blocked_range<size_t> r) {
                              #ifdef CHECK_BOUND
counters count;
#endif
        //stack_allocator<Mutation_Count_Change> this_thread_FIFO_allocator(FIFO_allocator_state);
        for (size_t i = r.begin(); i < r.end(); i++) {
        //for (size_t i = 0; i < nodes_to_search.size(); i++) {
            if (search_context.is_group_execution_cancelled()) {
                break;
            }
    //        if(((!deferred_nodes.size())||(std::chrono::steady_clock::now()-last_save_time)<save_period)&&deferred_nodes.size()<max_queued_moves) {
                output_t out;
                out.moves=new std::vector<Profitable_Moves_ptr_t>;
                find_moves_bounded(nodes_to_search[i], out,radius
                              #ifdef CHECK_BOUND
                ,count
                #endif
                );
                //assert(this_thread_FIFO_allocator.empty());
                if (!out.moves->empty()) {
                    //resolve conflicts
                        deferred_nodes.push_back(nodes_to_search[i]);
                    resolver_node.try_put(out.moves);
                }else {
                    delete out.moves;
                }
                checked_nodes.fetch_add(1,std::memory_order_relaxed);
            //} else {
               // deferred_nodes.push_back(nodes_to_search[i]);
            //}
        }
     #ifdef CHECK_BOUND
        total+=count.total;
        saved+=count.saved;
    #endif
    },search_context);
    progress_meter.stop();
    auto searh_end=std::chrono::steady_clock::now();
    std::chrono::duration<double> elpased_time =searh_end-search_start;
    fprintf(stderr, "\nSearch took %f minutes on %d\n",elpased_time.count()/60,this_rank);
}
void optimize_tree_main_thread(std::vector<MAT::Node *> &nodes_to_search,
                                        MAT::Tree &t,int radius,FILE* log,bool allow_drift,int iteration,
                                        std::vector<MAT::Node*>& deferred_nodes_out,bool MPI_involved
#ifndef NDEBUG
                     , Original_State_t& origin_states
#endif
                                       ){
    t.breadth_first_expansion();
    auto dfs_ordered_nodes=t.depth_first_expansion();
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", dfs_ordered_nodes.size());
    //for resolving conflicting moves
    Deferred_Move_t deferred_moves;
    Cross_t potential_crosses(dfs_ordered_nodes.size(),nullptr);
    tbb::flow::graph resolver_g;
    resolver_node_t resover_node(resolver_g,1,Conflict_Resolver(potential_crosses,deferred_moves),tbb::flow::queueing());
    std::thread move_reciever(MPI_recieve_move,std::ref(dfs_ordered_nodes),std::ref(resover_node));
    //progress bar
    tbb::concurrent_vector<MAT::Node*> deferred_nodes_this_thread;
    search(nodes_to_search,t,radius,deferred_nodes_this_thread,resover_node);
    std::vector<std::string> defered_node_identifier;
    if (MPI_involved) {
        MPI_Barrier(MPI_COMM_WORLD);
        int temp;    
        MPI_Send(&temp, 0, MPI_INT, 0, MOVE_TAG, MPI_COMM_WORLD);
        fprintf(stderr, "Sent finish msg\n");
        std::vector<int> other_process_deferred_size(process_count);
        size_t self_count=0;
        MPI_Gather(&self_count, 1, MPI_INT, other_process_deferred_size.data(), 1,MPI_INT, 0 ,MPI_COMM_WORLD);
        std::vector<int> recieve_offsets(process_count);
        size_t offset_so_far=0;
        for (int proc_idx=0; proc_idx<process_count; proc_idx++) {
            recieve_offsets[proc_idx]=offset_so_far;
            fprintf(stderr, "Expecting %d nodes form %d\n",other_process_deferred_size[proc_idx],proc_idx);
            offset_so_far+=other_process_deferred_size[proc_idx];
        }
        std::vector<size_t> other_defered_nodes_dfs_idx(offset_so_far);
        fprintf(stderr, "Expecting %zu nodes in total\n",offset_so_far);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(&offset_so_far, 0, MPI_UNSIGNED_LONG, other_defered_nodes_dfs_idx.data(),other_process_deferred_size.data(), recieve_offsets.data(),MPI_UNSIGNED_LONG, 0,MPI_COMM_WORLD);
        fprintf(stderr, "reserving %zu nodes in total\n",deferred_nodes_this_thread.size()+offset_so_far);
        defered_node_identifier.reserve(deferred_nodes_this_thread.size()+offset_so_far);
        for (const auto node : deferred_nodes_this_thread) {
            defered_node_identifier.emplace_back(node->identifier);
        }
        for(auto idx:other_defered_nodes_dfs_idx){
            defered_node_identifier.push_back(dfs_ordered_nodes[idx]->identifier);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    fprintf(stderr, "finished recieving defered nodes\n");

    move_reciever.join();
    resolver_g.wait_for_all();
    //apply moves
    #ifdef CHECK_BOUND
    fprintf(stderr, "Total %lu arcs, saved %lu arcs\n",total.load(),saved.load());
    #endif
    auto apply_start=std::chrono::steady_clock::now();
    fputs("Start applying moves\n",stderr);
    std::vector<Profitable_Moves_ptr_t> all_moves;
    schedule_moves(potential_crosses,all_moves);
    fprintf(stderr, "Applying %zu moves\n",all_moves.size());
    if (iteration>0) {
        log_move_detail(all_moves, log, iteration, radius);
    }
    apply_moves(all_moves, t
#ifdef CHECK_STATE_REASSIGN
                ,origin_states
#endif
               );
    auto apply_end=std::chrono::steady_clock::now();
    auto elpased_time =apply_end-apply_start;
    fprintf(stderr, "apply moves took %ld seconds\n",elpased_time.count());
    //recycle conflicting moves
    int init_deferred=deferred_moves.size();
    int recycled=0;
    fputs("Start recycling conflicting moves\n",stderr);
    while (!deferred_moves.empty()&&(!allow_drift)) {
        {
            Deferred_Move_t deferred_moves_next;
            potential_crosses.clear();
            auto bfs_ordered_nodes=t.breadth_first_expansion();
            potential_crosses.resize(bfs_ordered_nodes.size(),nullptr);
            tbb::flow::graph resolver_g;
            resolver_node_t resover_node(resolver_g,1,Conflict_Resolver(potential_crosses,deferred_moves_next),tbb::flow::queueing());
            if (timed_print_progress) {
                fprintf(stderr,"\rrecycling conflicting moves, %zu left",deferred_moves.size());
            }
            tbb::parallel_for(tbb::blocked_range<size_t>(0,deferred_moves.size()),[&deferred_moves,&resover_node,&t](const tbb::blocked_range<size_t>& r) {
                for (size_t i=r.begin(); i<r.end(); i++) {
                    MAT::Node* src=t.get_node(deferred_moves[i].first);
                    MAT::Node* dst=t.get_node(deferred_moves[i].second);
                    if (src&&dst) {
                        output_t out;
                        out.moves=new std::vector<Profitable_Moves_ptr_t>;
                        if (check_not_ancestor(dst, src)) {
                            individual_move(src,dst,get_LCA(src, dst),out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                            ,&t
#endif
                                           );
                        }
                        if (!out.moves->empty()) {
                            resover_node.try_put(out.moves);
                        }else{
                            delete out.moves;
                        }
                    }

                }
            });
            all_moves.clear();
            resolver_g.wait_for_all();
            schedule_moves(potential_crosses,all_moves);
            if (iteration>0) {
                log_move_detail(all_moves, log, iteration, radius);
            }
            recycled+=all_moves.size();
            apply_moves(all_moves, t
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
    //check_samples(t.root, origin_states, &t);
#endif
    deferred_nodes_out.clear();
    deferred_nodes_out.reserve(defered_node_identifier.size());
    for (const auto& id : defered_node_identifier) {
        auto iter=t.all_nodes.find(id);
        if (iter!=t.all_nodes.end()) {
            deferred_nodes_out.push_back(iter->second);
        }
    }
    fprintf(stderr, "recycled %f of conflicting moves \n",(double)recycled/(double)init_deferred);
    fprintf(stderr, "recycling moves took %ld seconds\n",elpased_time.count());
}
void optimize_tree_worker_thread(MAT::Tree &t,int radius,std::vector<MAT::Node *> &nodes_to_search){
    tbb::concurrent_vector<MAT::Node*> defered_nodes;
    t.depth_first_expansion();
    tbb::flow::graph sender_g;
    resolver_node_t resolver_node(sender_g,1,MPI_move_sender());
    search(nodes_to_search, t, radius, defered_nodes, resolver_node);
    sender_g.wait_for_all();
    MPI_Barrier(MPI_COMM_WORLD);
    int defered_size=defered_nodes.size();
    fprintf(stderr, "Sending %d deferend nodes form %d\n",defered_size,this_rank);
    MPI_Gather(&defered_size, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<size_t> deferred_idx;
    deferred_idx.reserve(defered_nodes.size());
    for (const auto node : defered_nodes) {
        deferred_idx.push_back(node->dfs_index);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(deferred_idx.data(), defered_size, MPI_UNSIGNED_LONG, NULL, NULL, NULL, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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
