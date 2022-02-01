#include "place_sample.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <tbb/concurrent_queue.h>
#include <tbb/flow_graph.h>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <mpi.h>

static void gather_par_mutation_step(std::unordered_map<int, int>& to_find,const MAT::Mutations_Collection& upstream, MAT::Mutations_Collection& output){
    for (const auto& mut : upstream) {
        auto iter=to_find.find(mut.get_position());
        if (iter!=to_find.end()) {
            output[iter->second].set_par_one_hot(mut.get_mut_one_hot());
            to_find.erase(iter);
        }
    }
}
static void gather_par_mut(std::unordered_map<int, int>& to_find,MAT::Node*& node, MAT::Mutations_Collection& output){
    while (node&&(!to_find.empty())) {
        gather_par_mutation_step(to_find,node->mutations,output);
        node=node->parent;
    }
    for(auto& temp: to_find){
        output[temp.second].set_par_one_hot(output[temp.second].get_ref_one_hot());
    }
}

static void discretize_mutations(const std::vector<To_Place_Sample_Mutation> &in,
                                 const MAT::Mutations_Collection &shared_mutations,
                                 MAT::Node *parent_node,
                                 MAT::Mutations_Collection &out) {
    out.reserve(in.size());
    std::unordered_map<int, int> par_nuc_idx;
    assert(in.back().position==INT_MAX);
    for (size_t idx=0;idx<(in.size()-1) ; idx++) {
        const auto & mut=in[idx];
        if (mut.mut_nuc == 0xf) {
            for (int pos = mut.position; pos <= mut.get_end_range(); pos++) {
                par_nuc_idx.emplace(pos, out.size());
                out.push_back(MAT::Mutation(mut.chrom_idx, pos, 0, 0xf));
                out.back().set_descendant_mut(0xf);
            }
        } else {
            out.push_back(MAT::Mutation(mut.chrom_idx, mut.position,
                                        mut.par_nuc, mut.mut_nuc));
            out.back().set_descendant_mut(mut.mut_nuc);
        }
    }
    gather_par_mutation_step(par_nuc_idx, shared_mutations,out);
    gather_par_mut(par_nuc_idx, parent_node, out);
}

static std::tuple<std::vector<Main_Tree_Target>, size_t,bool> *
deserialize_move(char *in, int size, MAT::Tree &tree) {
    Mutation_Detailed::search_result result;
    bool parse_result=result.ParseFromArray((void *)in, size);
    if (!parse_result) {
        fprintf(stderr, "parse error");
        raise(SIGTRAP);
    }
    auto out = new std::tuple<std::vector<Main_Tree_Target>, size_t,bool>;
    std::get<2>(*out)=false;
    auto &targets = std::get<0>(*out);
    std::get<1>(*out) = result.sample_id();
    if (tree.get_node_name(std::get<1>(*out))=="") {
        fprintf(stderr, "node %zu no name\n",std::get<1>(*out));
        raise(SIGTRAP);
    }
    auto target_count = result.place_targets_size();
    targets.reserve(target_count);
    for (int count = 0; count < target_count; count++) {
        Main_Tree_Target target;
        auto &parsed_target = result.place_targets(count);
        target.target_node = tree.get_node(parsed_target.target_node_id());
        if (target.target_node->parent==nullptr) {
            fprintf(stderr, "unset patetn\n");
            raise(SIGTRAP);
        }
        target.parent_node = tree.get_node(parsed_target.parent_node_id());
        load_mutations(parsed_target.sample_mutation_positions(),
                       parsed_target.sample_mutation_other_fields(),
                       target.sample_mutations);
        load_mutations(parsed_target.split_mutation_positions(),
                       parsed_target.split_mutation_other_fields(),
                       target.splited_mutations.mutations);
        load_mutations(parsed_target.shared_mutation_positions(),
                       parsed_target.shared_mutation_other_fields(),
                       target.shared_mutations.mutations);
        if (parsed_target.sample_mutation_positions_size()==0) {
        raise(SIGTRAP);
        }
        if (parsed_target.sample_mutation_positions_size()!=target.sample_mutations.size()) {
        fprintf(stderr, "sample mut mismatch %d, loaded %d\n",parsed_target.sample_mutation_positions_size(),target.sample_mutations.size());
        raise(SIGTRAP);
        }
        targets.emplace_back(std::move(target));
    }
    return out;
}
typedef tbb::flow::function_node<move_type *,size_t> placer_node_t;
typedef tbb::flow::buffer_node<std::tuple<std::vector<Main_Tree_Target>, size_t,bool> *> Other_Proc_Buffer_T;
static void recieve_place(Other_Proc_Buffer_T &handler,
                          MAT::Tree &tree,int processes_left) {
    MPI_Status status;
    while (true) {
        MPI_Probe(MPI_ANY_SOURCE, PROPOSED_PLACE, MPI_COMM_WORLD, &status);
        int msg_size;
        MPI_Get_count(&status, MPI_BYTE, &msg_size);
        if (msg_size == 0) {
            processes_left--;
            if (processes_left) {
                continue;
            }else {
                break;            
            }
        }
        auto temp = new char[msg_size];
        mpi_trace_print("recieving proposed place from %d\n",status.MPI_SOURCE);
        MPI_Recv(temp, msg_size, MPI_BYTE, status.MPI_SOURCE, PROPOSED_PLACE, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        auto recieved = deserialize_move(temp, msg_size, tree);
        delete[] temp;
        mpi_trace_print("recieved proposed place\n");
        handler.try_put(recieved);
    }
}
static std::string *
serialize_placed_sample(const MAT::Mutations_Collection &sample_mutations,
                        const MAT::Mutations_Collection &splited_mutations,
                        const MAT::Mutations_Collection &shared_mutations,
                        size_t target_id, size_t split_id, size_t sample_id) {
    Mutation_Detailed::placed_target target;
    target.set_sample_id(sample_id);
    target.set_target_node_id(target_id);
    target.set_split_node_id(split_id);
    fill_mutation_vect(target.mutable_sample_mutation_positions(),
                       target.mutable_sample_mutation_other_fields(),
                       sample_mutations);
    fill_mutation_vect(target.mutable_split_mutation_positions(),
                       target.mutable_split_mutation_other_fields(),
                       splited_mutations);
    fill_mutation_vect(target.mutable_shared_mutation_positions(),
                       target.mutable_shared_mutation_other_fields(),
                       shared_mutations);
    auto out = new std::string;
    *out = target.SerializeAsString();
    return out;
}

static void send_paced_move(tbb::concurrent_bounded_queue<std::string*>& send_queue){
    while (true) {
        std::string* in;
        send_queue.pop(in);
        if (in==nullptr) {
            fprintf(stderr, "send queue exit\n");
            break;
        }
        auto size = in->size();
        MPI_Bcast((void *)&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        mpi_trace_print("Main sending move of size %zu\n",size);
        MPI_Bcast((void *)in->c_str(), size, MPI_BYTE, 0, MPI_COMM_WORLD);
        mpi_trace_print("Main sent move\n");
        //fprintf(stderr, "send queue size %zu\n",send_queue.size());
        delete in;
    }
}
float redo;
float total;
struct main_tree_new_place_handler {
    MAT::Tree &main_tree;
    std::vector<MAT::Node *> &deleted_nodes;
    tbb::concurrent_bounded_queue<std::string*>& send_queue;
    size_t
    operator()(move_type *in) {
	auto &search_result = std::get<0>(*in);
        auto idx = std::get<1>(*in);
        for (const auto &placement : search_result) {
            if (placement.parent_node == nullptr ||
                placement.target_node == nullptr ||
                placement.parent_node != placement.target_node->parent) {
                if (placement.parent_node&&placement.target_node) {
                    auto actual_par=placement.target_node->parent;
                    fprintf(stderr, "par mismatch: recieved %zu, actual %zu, target %zu",placement.parent_node->node_id,actual_par?actual_par->node_id:-1,placement.target_node->node_id);
                    if (std::get<2>(*in)) {
                        fprintf(stderr, " from self\n");
                    }else {
                        fprintf(stderr, " from other\n");
                    }
                }else {
                    fprintf(stderr, "Node not found\n");
                }
                fprintf(stderr, "Redoing %s \n",
                        main_tree.get_node_name(idx).c_str());
                delete in;
                redo++;
		float redo_rate=redo/total;
		fprintf(stderr,"%f%% redon, %f redos, %f total\n",redo_rate*100,redo,total);
		return idx;
            }
        }
        total++;
        auto smallest_idx = search_result[0].target_node->node_id;
        auto arg_small_idx = 0;
        for (int target_idx = 1; target_idx < search_result.size();
             target_idx++) {
            if (smallest_idx > search_result[0].target_node->node_id) {
                smallest_idx = search_result[target_idx].target_node->node_id;
                arg_small_idx = target_idx;
            }
        }
        auto target = search_result[arg_small_idx];
        auto target_idx = target.target_node->node_id;
        MAT::Mutations_Collection sample_mutations;
        discretize_mutations(target.sample_mutations, target.shared_mutations,
                             target.parent_node, sample_mutations);
        auto out = update_main_tree(sample_mutations, target.splited_mutations,
                                    target.shared_mutations, target.target_node,
                                    idx, main_tree, 0);
        if (out.deleted_nodes) {
            deleted_nodes.push_back(out.deleted_nodes);
        }
        auto split_id = out.splitted_node ? out.splitted_node->node_id : 0;
        send_queue.push(serialize_placed_sample(
            sample_mutations, target.splited_mutations, target.shared_mutations,
            target_idx, split_id, idx));
        delete in;
        return SIZE_MAX;
    }
};
static void main_tree_distribute_samples(std::atomic_uint64_t& curr_idx,std::vector<Sample_Muts>& to_place,int processes_left){
    int ignore;
    while (true) {
        MPI_Status status;
        MPI_Recv(&ignore, 0, MPI_BYTE, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &status);
        mpi_trace_print("main recieved work req \n");
        int reply_dest=status.MPI_SOURCE;
        Mutation_Detailed::sample_to_place temp;
        auto send_idx=curr_idx++;
        mpi_trace_print("cur send idx %zu \n",send_idx);
        if (send_idx>=to_place.size()) { 
            MPI_Send(&ignore, 0, MPI_BYTE, reply_dest, WORK_RES_TAG, MPI_COMM_WORLD);
            processes_left--;
            break;
        }
        auto to_send=to_place[send_idx];
        temp.set_sample_id(to_send.sample_idx);
        fill_mutation_vect(temp.mutable_sample_mutation_positions(), temp.mutable_sample_mutation_other_fields(), to_send.muts);
        auto buffer=temp.SerializeAsString();
        mpi_trace_print( "main sending work res \n");
        MPI_Send(buffer.c_str(), buffer.size(), MPI_BYTE, reply_dest, WORK_RES_TAG, MPI_COMM_WORLD);
        mpi_trace_print("main sent work res \n");
    }
    for(;processes_left>0;processes_left--){
        MPI_Status status;
        MPI_Recv(&ignore, 0, MPI_BYTE, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &status);
        int reply_dest=status.MPI_SOURCE;
        MPI_Send(&ignore, 0, MPI_BYTE, reply_dest, WORK_RES_TAG, MPI_COMM_WORLD);    
    }
}

typedef tbb::flow::multifunction_node<size_t,tbb::flow::tuple<Sample_Muts*>> Pusher_Node_T;
struct Pusher{
    std::atomic_uint64_t& curr_idx;
    std::vector<Sample_Muts>& to_place;

    void operator()(size_t idx_in,Pusher_Node_T::output_ports_type& out ){
        if (idx_in==SIZE_MAX) {
            auto send_idx=curr_idx++;
            if (send_idx<to_place.size()) {
                std::get<0>(out).try_put(&to_place[send_idx]);
            }
            return;
        }
        std::get<0>(out).try_put(&to_place[idx_in-to_place[0].sample_idx]);
    }
};
void place_sample_leader(std::vector<Sample_Muts> &sample_to_place, MAT::Tree &main_tree,int batch_size,int proc_count) {
    redo=0;
    total=0;
    std::vector<MAT::Node*> deleted_nodes;
    deleted_nodes.reserve(sample_to_place.size());
    auto start_time=std::chrono::steady_clock::now();
    std::thread *listen_thread;
    std::thread *distribute_thread;
    std::thread *sender_thread;
    std::atomic_uint64_t curr_idx(0);
    tbb::concurrent_bounded_queue<std::string*> send_queue;
        tbb::flow::graph g;
        Pusher_Node_T init(g, 1, Pusher{curr_idx, sample_to_place});
        tbb::flow::function_node<
            Sample_Muts*, move_type *>
            searcher(g, tbb::flow::unlimited,
                     Finder{main_tree,false});
        tbb::flow::make_edge(std::get<0>(init.output_ports()),searcher);
        placer_node_t placer(g, 1,
                     main_tree_new_place_handler{main_tree,deleted_nodes,send_queue});
        tbb::flow::make_edge(searcher,placer);
        tbb::flow::make_edge(placer,init);
        Other_Proc_Buffer_T buffer(g);
        tbb::flow::make_edge(buffer,placer);
        if (proc_count>0) {
            sender_thread= new std::thread(send_paced_move,std::ref(send_queue));
            MPI_Barrier(MPI_COMM_WORLD);
            listen_thread=new std::thread (recieve_place,std::ref(buffer),std::ref(main_tree),proc_count-1);
            distribute_thread= new std::thread(main_tree_distribute_samples,std::ref(curr_idx), std::ref(sample_to_place),proc_count-1);
            
        }
        for(int temp=0;temp<batch_size;temp++){
            init.try_put(SIZE_MAX);
        }
        g.wait_for_all();
        send_queue.push(nullptr);
        for(auto node:deleted_nodes){
            delete node;
        }
        deleted_nodes.clear();
        if (proc_count>0) {
            listen_thread->join();
            distribute_thread->join();    
            delete listen_thread;
            delete distribute_thread;
        }       
}
