#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "priority_conflict_resolver.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
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
float update_rate=0.01;
#define MAX_MOVE_SIZE ((size_t)0x1000)
#define MAX_MOVE_MSG_SIZE (1+4*MAX_MOVE_SIZE)
thread_local allocator_state<Mutation_Count_Change> FIFO_allocator_state;
extern tbb::task_group_context search_context;
size_t nodes_per_min_per_thread=100;
float target_fetch_period=0.5;
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
void log_move_detail(const std::vector<Profitable_Moves_ptr_t> & moves, FILE* out,int iteration,int radius) {
    for(const auto& move:moves) {
        fprintf(out, "%s\t%s\t%d\t%d\t%d\t%lu\n",move->src->identifier.c_str(),move->dst->identifier.c_str()
                ,iteration,move->score_change,radius-move->radius_left,move->src->dfs_end_index-move->src->dfs_index);
    }
}

typedef tbb::flow::function_node<std::vector<Profitable_Moves_ptr_t>*,tbb::flow::continue_msg,tbb::flow::queueing> resolver_node_t;
static void MPI_recieve_move(const std::vector<MAT::Node*>& dfs_ordered_nodes,resolver_node_t& resover_node ) {
    int* buffer=new int[MAX_MOVE_MSG_SIZE];
    while (true) {
        std::vector<Profitable_Moves_ptr_t>* out =new std::vector<Profitable_Moves_ptr_t>;
        MPI_Status stat;
        MPI_Recv(buffer,MAX_MOVE_MSG_SIZE, MPI_INT, MPI_ANY_SOURCE, MOVE_TAG, MPI_COMM_WORLD, &stat);
        int msg_size;
        MPI_Get_count(&stat, MPI_INT,&msg_size);
        if (msg_size==0) {
            delete[] buffer;
            fprintf(stderr, "+++++++++Move receiver exit\n");
            return;
        }
        int moves_size=(msg_size-1)/4;
        //fprintf(stderr, "Recieved %d moves\n",moves_size);
        MAT::Node* src=dfs_ordered_nodes[buffer[0]];
        out->reserve(moves_size);
        for (int move_idx=0; move_idx<moves_size; move_idx++) {
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
struct MPI_move_sender {
    int* buffer;
    void init() {
        buffer=new int[1+4*MAX_MOVE_SIZE];
    }
    MPI_move_sender() {
        init();
    }
    MPI_move_sender(const MPI_move_sender&) {
        init();
    }
    void operator()(std::vector<Profitable_Moves_ptr_t>* to_send) {
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
    ~MPI_move_sender() {
        delete[] buffer;
    }
};
typedef tbb::flow::multifunction_node<std::vector<size_t>*, tbb::flow::tuple<std::vector<Profitable_Moves_ptr_t>*>,tbb::flow::rejecting> searcher_node_t;
struct move_searcher {
    const std::vector<MAT::Node*>& dfs_ordered_nodes;
    int radius;
    bool do_drift;
    Reachable reachable;
    void operator()(std::vector<size_t>* to_search,searcher_node_t::output_ports_type& output)const {
        int r=radius;
        auto start_time=std::chrono::steady_clock::now();
        for (auto idx:*to_search) {
            auto node_to_search=dfs_ordered_nodes[idx];
            output_t out;
            out.moves=new std::vector<Profitable_Moves_ptr_t>;
            find_moves_bounded(node_to_search, out,r,do_drift,reachable
#ifdef CHECK_BOUND
                               ,count
#endif
                              );
            if (!out.moves->empty()) {
                //resolve conflicts
                std::get<0>(output).try_put(out.moves);
            } else {
                delete out.moves;
            }
            //} else {
            // deferred_idx.push_back(nodes_to_search[i]);
            //}
        }
        float seconds_duration=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_time).count();
        float currate=to_search->size()/((seconds_duration+1.0)/60.0);
        nodes_per_min_per_thread=nodes_per_min_per_thread*(1-update_rate)+update_rate*currate;
        //fprintf(stderr, "Cur rate %f, %zu\n",currate,nodes_per_min_per_thread);
        delete to_search;
    }
};
static void node_distributor(const std::vector<size_t>& node_to_search_idx,std::atomic<bool>& done,std::vector<size_t>& nodes_not_searched,std::chrono::steady_clock::time_point stop_time,bool is_one_proc) {
    size_t idx=0;
    auto stop_in_min=std::chrono::duration_cast<std::chrono::minutes>(stop_time-std::chrono::steady_clock::now()).count();
    fprintf(stderr,"Will stop in %zu min\n",stop_in_min );
    std::vector<size_t> rates(process_count,1000);
    size_t total_rate=process_count*1000;
    while (idx<node_to_search_idx.size()) {
        size_t curr_proc_rate;
        MPI_Status stat;
        MPI_Recv(&curr_proc_rate, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &stat);
        total_rate+=(curr_proc_rate-rates[stat.MPI_SOURCE]);
        rates[stat.MPI_SOURCE]=curr_proc_rate;
        size_t remaining_nodes=node_to_search_idx.size()-idx;
        float time_left=(float)remaining_nodes/(float)total_rate;
        float release_time=std::max(0.1f,time_left/2);
        size_t release_node_count=std::min(size_t(1+release_time*curr_proc_rate),curr_proc_rate);
        size_t count_to_send=std::min(release_node_count,remaining_nodes);
        fprintf(stderr, " %zu nodes left, %0.1f min left\n",remaining_nodes,time_left);
        MPI_Send(node_to_search_idx.data()+idx, count_to_send, MPI_UNSIGNED_LONG, stat.MPI_SOURCE, WORK_RES_TAG, MPI_COMM_WORLD);
        idx+=count_to_send;
        if (std::chrono::steady_clock::now()>=stop_time) {
            fprintf(stderr, "================timeout=========\n");
            nodes_not_searched.insert(nodes_not_searched.end(),node_to_search_idx.begin()+idx,node_to_search_idx.end());
            break;
        }
        if (interrupted) {
            fprintf(stderr, "================interrupted=========\n");
            break;
        }
    }
    for (int zero_sent=0; zero_sent<(is_one_proc?1:process_count); zero_sent++) {
        size_t count_to_send=0;
        MPI_Status stat;
        MPI_Recv(&count_to_send, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &stat);
        MPI_Send(node_to_search_idx.data(), 0, MPI_UNSIGNED_LONG, stat.MPI_SOURCE, WORK_RES_TAG, MPI_COMM_WORLD);
    }
    fprintf(stderr, "distributor exit\n");
}
struct fetcher {
    std::vector<size_t>& nodes_to_push;
    mutable size_t release_rate;
    //mutable int is_longer_count;
    std::chrono::steady_clock::time_point& last_request_time;
    fetcher(std::vector<size_t>& nodes_to_push,std::chrono::steady_clock::time_point& last_request_time):nodes_to_push(nodes_to_push),last_request_time(last_request_time) {
        nodes_per_min_per_thread=100;
        release_rate=100;
        update_rate=0.1;
    }
    bool operator()(std::vector<size_t>*& out) const {
        if (nodes_to_push.empty()) {
            size_t req_size=num_threads*nodes_per_min_per_thread;
            auto this_request_time=std::chrono::steady_clock::now();
            MPI_Send(&req_size, 1, MPI_UNSIGNED_LONG, 0, WORK_REQ_TAG, MPI_COMM_WORLD);
            MPI_Status stat;
            nodes_to_push.resize(req_size);
            MPI_Recv(nodes_to_push.data(), req_size, MPI_UNSIGNED_LONG, 0, WORK_RES_TAG, MPI_COMM_WORLD, &stat);
            int recieve_count;
            MPI_Get_count(&stat, MPI_UNSIGNED_LONG, &recieve_count);
            auto request_period=std::chrono::duration_cast<std::chrono::seconds>(this_request_time-last_request_time).count();
            /*if (request_period>60) {
                is_longer_count++;
            }else {
                is_longer_count--;
            }
            if (abs(is_longer_count)<3) {
                update_rate=std::max(0.01,update_rate-0.01);
            }else {
                update_rate=std::min(0.5,update_rate+0.01);
            }*/
            fprintf(stderr, "requesting %zu nodes from %d after %ld seconds, got %d nodes \n",req_size,this_rank,request_period,recieve_count);
            release_rate=1+recieve_count/num_threads;
            last_request_time=this_request_time;
            if (recieve_count==0) {
                fprintf(stderr, "fetcher exit\n");
                return false;
            }
            nodes_to_push.resize(recieve_count);
        }
        auto nodes_to_release_this_round=std::min(nodes_to_push.size(),release_rate);
        //fprintf(stderr, "buf size %zu, releasing %lu nodes \n",nodes_to_push.size(),nodes_to_release_this_round);
        auto split_iter=nodes_to_push.end()-nodes_to_release_this_round;
        out=new std::vector<size_t>(split_iter,nodes_to_push.end());
        nodes_to_push.erase(split_iter,nodes_to_push.end());
        //fprintf(stderr, "left %zu nodes at %d \n",nodes_to_push.size(),this_rank);
        return true;
    }
};
Reachable set_reachable(int radius,MAT::Tree& t,bool search_all_dir) {
    //auto height=t.get_max_level();
    Reachable reachable;
    reachable.always_search=search_all_dir;
    //reachable.reachable_change=radius<(2*height);
    reachable.reachable_change=true;
    return reachable;
}
void optimize_tree_main_thread(std::vector<size_t> &nodes_to_search,
                               MAT::Tree &t,int radius,FILE* log,bool allow_drift,int iteration,
                               std::vector<MAT::Node*>& deferred_nodes_out,bool MPI_involved,std::chrono::steady_clock::time_point end_time,bool do_continue,bool search_all_dir,bool isfirst_this_iter
#ifndef NDEBUG
                               , Original_State_t& origin_states
#endif
                              ) {
    t.breadth_first_expansion();
    auto dfs_ordered_nodes=t.depth_first_expansion();
    auto start_time=std::chrono::steady_clock::now();
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", dfs_ordered_nodes.size());
    std::atomic<bool> done(false);
    std::vector<size_t> incomplete_idx;
    std::thread distributor_thread(node_distributor,std::ref(nodes_to_search), std::ref(done),std::ref(incomplete_idx),end_time,!MPI_involved);
    //for resolving conflicting moves
    Deferred_Move_t deferred_moves;
    Cross_t potential_crosses(dfs_ordered_nodes.size(),nullptr);
    tbb::flow::graph g;
    std::vector<std::string> defered_node_identifier;
    resolver_node_t resover_node(g, 1,
                                 Conflict_Resolver(potential_crosses,
                                         deferred_moves,
                                         &defered_node_identifier));
    std::thread move_reciever(MPI_recieve_move,std::ref(dfs_ordered_nodes),std::ref(resover_node));
    //progress bar
    searcher_node_t searcher(g,num_threads+1,move_searcher{dfs_ordered_nodes,radius,allow_drift,set_reachable(radius, t,search_all_dir)});
    tbb::flow::make_edge(std::get<0>(searcher.output_ports()),resover_node);
    std::vector<size_t> local_nodes_to_search;
    auto last_request_time=std::chrono::steady_clock::now();
    tbb::flow::source_node<std::vector<size_t>*> fetcher_node(g,fetcher(local_nodes_to_search,last_request_time));
    tbb::flow::make_edge(fetcher_node,searcher);
    if (MPI_involved) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int temp;
    MPI_Send(&temp, 0, MPI_INT, 0, MOVE_TAG, MPI_COMM_WORLD);
    //fprintf(stderr, "Sent finish msg\n");
    move_reciever.join();
    g.wait_for_all();
    done.store(true);
    //fprintf(stderr, "Waiting for distributor thread\n");
    distributor_thread.join();
    if(do_continue) {
        defered_node_identifier.reserve(defered_node_identifier.size()+incomplete_idx.size());
        for(auto idx:incomplete_idx) {
            defered_node_identifier.push_back(dfs_ordered_nodes[idx]->identifier);
        }
    }
    double search_min=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-start_time).count();
    fprintf(stderr, "Search took %f min \n",search_min/60000.0);
    //apply moves
#ifdef CHECK_BOUND
    fprintf(stderr, "Total %lu arcs, saved %lu arcs\n",total.load(),saved.load());
#endif
    if (isfirst_this_iter) {
        for (auto node : dfs_ordered_nodes) {
            node->clear_changed();
        }
    }
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
    auto elpased_time =std::chrono::duration_cast<std::chrono::seconds>(apply_end-apply_start);
    fprintf(stderr, "apply moves took %ld seconds\n",elpased_time.count());
    //recycle conflicting moves
    int init_deferred=deferred_moves.size();
    int recycled=0;
    fputs("Start recycling conflicting moves\n",stderr);
    while (!deferred_moves.empty()&&(!allow_drift)) {
        {
            fprintf(stderr, "\r %zu nodes left",deferred_moves.size());
            Deferred_Move_t deferred_moves_next;
            potential_crosses.clear();
            auto bfs_ordered_nodes=t.breadth_first_expansion();
            potential_crosses.resize(bfs_ordered_nodes.size(),nullptr);
            tbb::flow::graph resolver_g;
            std::vector<MAT::Node*> ignored;
            resolver_node_t resover_node(resolver_g,1,Conflict_Resolver(potential_crosses,deferred_moves_next,nullptr));
            tbb::parallel_for(tbb::blocked_range<size_t>(0,deferred_moves.size()),[&deferred_moves,&resover_node,&t,allow_drift](const tbb::blocked_range<size_t>& r) {
                for (size_t i=r.begin(); i<r.end(); i++) {
                    MAT::Node* src=t.get_node(deferred_moves[i].first);
                    MAT::Node* dst=t.get_node(deferred_moves[i].second);
                    if (src&&dst) {
                        output_t out;
                        out.moves=new std::vector<Profitable_Moves_ptr_t>;
                        if (check_not_ancestor(dst, src)) {
                            individual_move(src,dst,get_LCA(src, dst),out,allow_drift
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                                            ,&t
#endif
                                           );
                        }
                        if (!out.moves->empty()) {
                            resover_node.try_put(out.moves);
                        } else {
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
    elpased_time =std::chrono::duration_cast<std::chrono::seconds>(recycle_end-apply_end);
#ifndef NDEBUG
    //check_samples(t.root, origin_states, &t);
#endif
    if(!do_continue) {
        return;
    }
    clean_tree(t);
    deferred_nodes_out.clear();
    fprintf(stderr, "First stage %zu deferred node \n",defered_node_identifier.size());
    deferred_nodes_out.reserve(defered_node_identifier.size());
    for (const auto& id : defered_node_identifier) {
        auto iter=t.all_nodes.find(id);
        if (iter!=t.all_nodes.end()) {
            deferred_nodes_out.push_back(iter->second);
        }
    }
    fprintf(stderr, "recycled %f of conflicting moves \n",(double)recycled/(double)init_deferred);
    fprintf(stderr, "recycling moves took %ld seconds\n",elpased_time.count());
    t.populate_ignored_range();
}
void optimize_tree_worker_thread(MAT::Tree &t,int radius,bool do_drift,bool search_all_dir) {
    auto dfs_ordered_nodes=t.depth_first_expansion();
    tbb::flow::graph g;
    resolver_node_t resolver_node(g,1,MPI_move_sender());
    searcher_node_t searcher(g,num_threads+1,move_searcher{dfs_ordered_nodes,radius,do_drift,set_reachable(radius, t,search_all_dir)});
    tbb::flow::make_edge(std::get<0>(searcher.output_ports()),resolver_node);
    std::vector<size_t> nodes_to_search;
    auto last_request_time=std::chrono::steady_clock::now();
    tbb::flow::source_node<std::vector<size_t>*> fetcher_node(g,fetcher(nodes_to_search,last_request_time));
    tbb::flow::make_edge(fetcher_node,searcher);
    g.wait_for_all();
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "%d finished",this_rank);
}

tbb::concurrent_unordered_map<MAT::Mutation,
    tbb::concurrent_unordered_map<std::string, nuc_one_hot> *,
    Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
