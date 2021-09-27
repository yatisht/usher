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
size_t nodes_per_min_per_thread=100;
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
            progress_bar_cv.wait_for(done_lock,std::chrono::minutes(1));
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
        fprintf(stderr,"checked %d nodes, estimated %f minutes left,found %zu nodes "
                "profitable\n",
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
    int msg_size;
    MPI_Get_count(&stat, MPI_INT,&msg_size);
    if (msg_size==0) {
        fprintf(stderr, "+++++++++Move receiver exit\n");
        return;
    }
    int moves_size=(msg_size-1)/4;
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
    int* buffer;
    void init(){
        buffer=new int[1+4*MAX_MOVE_SIZE];
    }
    MPI_move_sender(){
        init();
    }
    MPI_move_sender(const MPI_move_sender&){
        init();
    }
    void operator()(std::vector<Profitable_Moves_ptr_t>* to_send){
        size_t moves_to_send=std::min(to_send->size(),MAX_MOVE_SIZE);
        size_t msg_size=1+4*moves_to_send;
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
        delete to_send;
    }
    ~MPI_move_sender(){
        delete[] buffer;
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
typedef tbb::flow::multifunction_node<std::vector<size_t>* , tbb::flow::tuple<std::vector<Profitable_Moves_ptr_t>*>> searcher_node_t;
struct move_searcher{
    const std::vector<MAT::Node*>& dfs_ordered_nodes;
    int radius;
    void operator()(std::vector<size_t>* to_search,searcher_node_t::output_ports_type& output)const {
        int r=radius;
        auto start_time=std::chrono::steady_clock::now();
        for (auto idx:*to_search) {
            auto node_to_search=dfs_ordered_nodes[idx];
            output_t out;
            out.moves=new std::vector<Profitable_Moves_ptr_t>;
            find_moves_bounded(node_to_search, out,r
                              #ifdef CHECK_BOUND
                ,count
                #endif
                );
                if (!out.moves->empty()) {
                    //resolve conflicts
                    std::get<0>(output).try_put(out.moves);
                }else {
                    delete out.moves;
                }
            //} else {
               // deferred_idx.push_back(nodes_to_search[i]);
            //}
        }
        float seconds_duration=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_time).count();
        float currate=to_search->size()/((seconds_duration+1.0)/60.0);
        nodes_per_min_per_thread=nodes_per_min_per_thread*0.99+0.01*currate;
        //fprintf(stderr, "Cur rate %f, %zu\n",currate,nodes_per_min_per_thread);
    }
};
static void node_distributor(const std::vector<size_t>& node_to_search_idx,std::atomic<bool>& done){
    size_t idx=0;
    while (idx<node_to_search_idx.size()) {
        size_t count_to_send;
        MPI_Status stat;
        MPI_Recv(&count_to_send, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &stat);
        count_to_send=std::min(count_to_send,node_to_search_idx.size()-idx);
        fprintf(stderr, "sending %zu nodes\n",count_to_send);
        MPI_Send(node_to_search_idx.data()+idx, count_to_send, MPI_UNSIGNED_LONG, stat.MPI_SOURCE, WORK_RES_TAG, MPI_COMM_WORLD);
        idx+=count_to_send;
    }
    for (int zero_sent=0;zero_sent<process_count;zero_sent++) {
        size_t count_to_send=0;
        MPI_Status stat;
        MPI_Recv(&count_to_send, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &stat);
        MPI_Send(node_to_search_idx.data(), 0, MPI_UNSIGNED_LONG, stat.MPI_SOURCE, WORK_RES_TAG, MPI_COMM_WORLD);
    }
    fprintf(stderr, "distributor exit\n");
}
struct fetcher{
    std::vector<size_t>& nodes_to_push;
    bool do_request;
    fetcher(std::vector<size_t>& nodes_to_push,bool use_MPI):nodes_to_push(nodes_to_push),do_request(use_MPI){}
    bool operator()(std::vector<size_t>*& out) const{
        auto nodes_count_this_call=nodes_to_push.size();
        auto call_start_time=std::chrono::steady_clock::now();
        if (nodes_to_push.empty()) {
            if (!do_request) {
                return false;
            }
            size_t req_size=num_threads*nodes_per_min_per_thread;
            MPI_Send(&req_size, 1, MPI_UNSIGNED_LONG, 0, WORK_REQ_TAG, MPI_COMM_WORLD);
            MPI_Status stat;
            nodes_to_push.resize(req_size);
            MPI_Recv(nodes_to_push.data(), req_size, MPI_UNSIGNED_LONG, 0, WORK_RES_TAG, MPI_COMM_WORLD, &stat);
            int recieve_count;
            MPI_Get_count(&stat, MPI_UNSIGNED_LONG, &recieve_count);
            fprintf(stderr, "requesting %zu nodes from %d, got %zu nodes \n",req_size,this_rank,recieve_count);
            if (recieve_count==0) {
                fprintf(stderr, "fetcher exit\n");
                return false;
            }
            nodes_to_push.resize(recieve_count);
        }
        auto nodes_to_release_this_round=std::min(nodes_to_push.size(),nodes_per_min_per_thread);
        //fprintf(stderr, "buf size %zu, releasing %lu nodes \n",nodes_to_push.size(),nodes_to_release_this_round);
        auto split_iter=nodes_to_push.end()-nodes_to_release_this_round;
        out=new std::vector<size_t>(split_iter,nodes_to_push.end());
        nodes_to_push.erase(split_iter,nodes_to_push.end());
        //fprintf(stderr, "left %zu nodes at %d \n",nodes_to_push.size(),this_rank);
        return true;
    }
};
void optimize_tree_main_thread(std::vector<size_t> &nodes_to_search,
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
    std::thread* distributor_thread;
    std::atomic<bool> done(false);
    if (MPI_involved) {
        distributor_thread=new std::thread(node_distributor,std::ref(nodes_to_search), std::ref(done));
    }
    //for resolving conflicting moves
    Deferred_Move_t deferred_moves;
    Cross_t potential_crosses(dfs_ordered_nodes.size(),nullptr);
    tbb::flow::graph g;
    std::vector<std::string> defered_node_identifier; 
    resolver_node_t resover_node(g, 1,
                                 Conflict_Resolver(potential_crosses,
                                                   deferred_moves,
                                                   &defered_node_identifier),
                                 tbb::flow::queueing());
    std::thread move_reciever(MPI_recieve_move,std::ref(dfs_ordered_nodes),std::ref(resover_node));
    //progress bar
    searcher_node_t searcher(g,num_threads+1,move_searcher{dfs_ordered_nodes,radius});
    tbb::flow::make_edge(std::get<0>(searcher.output_ports()),resover_node);
    std::vector<size_t> local_nodes_to_search;
    tbb::flow::source_node<std::vector<size_t>*> fetcher_node(g,fetcher(MPI_involved?local_nodes_to_search:nodes_to_search,MPI_involved));
    tbb::flow::make_edge(fetcher_node,searcher);
    if (MPI_involved) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int temp;
    MPI_Send(&temp, 0, MPI_INT, 0, MOVE_TAG, MPI_COMM_WORLD);
    fprintf(stderr, "Sent finish msg\n");
    move_reciever.join();
    g.wait_for_all();
    done.store(true);
    fprintf(stderr, "Waiting for distributor thread\n");
    if (MPI_involved) {
        distributor_thread->join();
        delete distributor_thread;
    }
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
            std::vector<MAT::Node*> ignored;
            resolver_node_t resover_node(resolver_g,1,Conflict_Resolver(potential_crosses,deferred_moves_next,nullptr),tbb::flow::queueing());
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
void optimize_tree_worker_thread(MAT::Tree &t,int radius){
    auto dfs_ordered_nodes=t.depth_first_expansion();
    tbb::flow::graph g;
    resolver_node_t resolver_node(g,1,MPI_move_sender(),tbb::flow::queueing());
    searcher_node_t searcher(g,num_threads+1,move_searcher{dfs_ordered_nodes,radius});
    tbb::flow::make_edge(std::get<0>(searcher.output_ports()),resolver_node);
    std::vector<size_t> nodes_to_search;
    tbb::flow::source_node<std::vector<size_t>*> fetcher_node(g,fetcher(nodes_to_search,true));
    tbb::flow::make_edge(fetcher_node,searcher);
    g.wait_for_all();
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "%d finished",this_rank);
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
