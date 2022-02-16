#include "place_sample.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <sched.h>
#include <string>
#include <tbb/concurrent_queue.h>
#include <tbb/flow_graph.h>
#include <tbb/task.h>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <mpi.h>
std::atomic_size_t backlog;
enum Status{OK,DONE,NOTHING};
//std::atomic_size_t backlog_prep;
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

void discretize_mutations(const std::vector<To_Place_Sample_Mutation> &in,
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
typedef std::pair<char*,int> pair_str;
move_type* deser_other_thread_move(char* in, int size,MAT::Tree &tree,std::vector<Sample_Muts>& to_place,size_t first_sample_idx,const std::vector<size_t>* idx_map)  {
    Mutation_Detailed::search_result result;
    bool parse_result=result.ParseFromArray((void *)in, size);
    if (!parse_result) {
        fprintf(stderr, "parse error");
        raise(SIGTRAP);
    }
    auto out = new move_type;
    std::get<2>(*out)=false;
    auto &targets = std::get<0>(*out);
    auto sample_idx=result.sample_id()-first_sample_idx;
    if (idx_map) {
        sample_idx=(*idx_map)[sample_idx];
    }
    /*if (sample_idx>to_place.size()) {
        fprintf(stderr, "Sample id out of range\n");
        raise(SIGTRAP);
    }*/
    std::get<1>(*out) = &to_place[sample_idx];
    /*if (tree.get_node_name(std::get<1>(*out))=="") {
        fprintf(stderr, "node %zu no name\n",std::get<1>(*out));
        raise(SIGTRAP);
    }*/
    auto target_count = result.place_targets_size();
    targets.reserve(target_count);
    for (int count = 0; count < target_count; count++) {
        Main_Tree_Target target;
        auto &parsed_target = result.place_targets(count);
        target.target_node = tree.get_node(parsed_target.target_node_id());
        /*if (target.target_node->parent==nullptr) {
            fprintf(stderr, "unset patetn\n");
            raise(SIGTRAP);
        }*/
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
        /*if (parsed_target.sample_mutation_positions_size()==0) {
        raise(SIGTRAP);
        }
        if (parsed_target.sample_mutation_positions_size()!=target.sample_mutations.size()) {
        fprintf(stderr, "sample mut mismatch %d, loaded %zu\n",parsed_target.sample_mutation_positions_size(),target.sample_mutations.size());
        raise(SIGTRAP);
        }*/
        targets.emplace_back(std::move(target));
    }
    delete[] in;
    return out;
    }
typedef tbb::concurrent_bounded_queue<move_type*> found_place_t;
struct Preped_Sample_To_Place;
struct Recieve_Place_State{
    found_place_t &handler;
    MAT::Tree &tree;
    int processes_left;
    std::vector<Sample_Muts>& to_place;
    const std::vector<size_t>* idx_map;
};
static Status recieve_place(Recieve_Place_State& state,size_t sample_start_idx) {
    MPI_Status status;
        int recieved;
        MPI_Iprobe(MPI_ANY_SOURCE, PROPOSED_PLACE, MPI_COMM_WORLD,&recieved, &status);
        if (!recieved) {
            return NOTHING;
        }
        //fprintf(stderr, "recieved\n");
        int msg_size;
        MPI_Get_count(&status, MPI_BYTE, &msg_size);
        auto temp = new char[msg_size];
        MPI_Recv(temp, msg_size, MPI_BYTE, status.MPI_SOURCE, PROPOSED_PLACE, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        if (msg_size == 0) {
 	    state.processes_left--;
	    fprintf(stderr, "reciever exiting %d proc left\n",state.processes_left);
            if (state.processes_left) {
                return OK;
            }else {
                return DONE;            
            }
        }
        //fprintf(stderr, "recieved proposed place from %d of size %d\n",status.MPI_SOURCE,msg_size);
        state.handler.push(deser_other_thread_move(temp, msg_size,state.tree,state.to_place,sample_start_idx,state.idx_map));
        return OK;
}
static std::string
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
    return target.SerializeAsString();
}
struct Preped_Sample_To_Place{
    size_t sample_id;
    MAT::Mutations_Collection shared_mutations;
    MAT::Mutations_Collection splited_mutations;
    MAT::Mutations_Collection sample_mutations;
    size_t target_id;
    size_t split_id;
};
typedef tbb::concurrent_bounded_queue<Preped_Sample_To_Place*> Placed_move_sended_state;
static Status send_paced_move(Placed_move_sended_state& send_queue){
        Preped_Sample_To_Place* in;
        auto did_pop=send_queue.try_pop(in);
        if (!did_pop) {
            return NOTHING;
        }
        if (in==nullptr) {
            fprintf(stderr, "send queue exit\n");
            return DONE;;
        }
        auto serialized=serialize_placed_sample(
            in->sample_mutations, in->splited_mutations, in->shared_mutations,
            in->target_id, in->split_id, in->sample_id);
        auto size = serialized.size();
        if (send_queue.size()>5) {
            fprintf(stderr, "queue size %ld\n",send_queue.size());        
        }
        MPI_Bcast((void *)&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        mpi_trace_print("Main sending move of size %zu\n",size);
        MPI_Bcast((void *)serialized.c_str(), size, MPI_BYTE, 0, MPI_COMM_WORLD);
        mpi_trace_print("Main sent move\n");
        //fprintf(stderr, "send queue size %zu\n",send_queue.size());
        delete in;
        return OK;
}
struct print_format{
    int parsimony_score;
    move_type* placement_info;
};

struct Print_Thread{
    MAT::Tree& tree;
    FILE* placement_stats_file;
    int max_parsimony;
    size_t max_uncertainty;
    size_t& node_count;
    std::vector<std::string>& low_confidence_samples;
    //std::vector<Sample_Muts*>& rejected_samples;
    std::vector<Clade_info>& samples_clade;
    size_t start_idx;
    void operator()(const print_format& in){
        auto sample_idx=std::get<1>(*in.placement_info)->sample_idx;
        auto sample_name = tree.get_node_name(sample_idx);
        auto sample_vec_idx=sample_idx-start_idx;
        const auto& search_result=std::get<0>(*in.placement_info);
        auto num_epps=search_result.size();
        fprintf(stderr,
                "Current tree size (#nodes): %zu\tSample name: %s\tParsimony "
                "score: %d\tNumber of parsimony-optimal placements: %zu\n",
                sample_idx, sample_name.c_str(), in.parsimony_score,
                num_epps);
        if (placement_stats_file) {
            fprintf(placement_stats_file, "%s\t%d\t%zu\t", sample_name.c_str(),
                    in.parsimony_score, num_epps);
        }
        if (num_epps > 1) {
            low_confidence_samples.emplace_back(sample_name);
            if (num_epps > max_uncertainty) {
                fprintf(
                    stderr,
                    "WARNING: Number of parsimony-optimal placements exceeds "
                    "maximum allowed value (%zu). Ignoring sample %s.\n",
                    max_uncertainty, sample_name.c_str());
                //rejected_samples.push_back(std::get<1>(*in.placement_info));
                delete in.placement_info;
                return;
            } else if (in.parsimony_score <= max_parsimony) {
                fprintf(stderr,
                        "WARNING: Multiple parsimony-optimal placements found. "
                        "Placement done without high confidence.\n");
            }
        }
        if (in.parsimony_score > max_parsimony) {
            fprintf(
                stderr,
                "WARNING: Parsimony score of the most parsimonious placement "
                "exceeds the maximum allowed value (%u). Ignoring sample %s.\n",
                max_parsimony, sample_name.c_str());
            //rejected_samples.push_back(std::get<1>(*in.placement_info));
            delete in.placement_info;
            return;
        }
        node_count++;
        auto& this_sample_clade=samples_clade[sample_vec_idx];
        this_sample_clade.valid=true;
        this_sample_clade.clade_assignments.clear();
        this_sample_clade.clade_assignments.resize(tree.get_num_annotations());
        this_sample_clade.best_clade_assignment.clear();
        this_sample_clade.best_clade_assignment.resize(tree.get_num_annotations());
        for (size_t c=0; c < tree.get_num_annotations(); c++) {
            this_sample_clade.clade_assignments[c].resize(search_result.size());
            //TODO: can be parallelized
            for (size_t k=0; k < search_result.size(); k++) {
                bool include_self = !search_result[k].target_node->is_leaf() && !search_result[k].splited_mutations.empty();
                auto clade_assignment = tree.get_clade_assignment(search_result[k].target_node, c, include_self);
                this_sample_clade.clade_assignments[c][k] = clade_assignment;
                if (k==0) {
                    this_sample_clade.best_clade_assignment[c] = clade_assignment;
                }
            }
            std::sort(this_sample_clade.clade_assignments[c].begin(), this_sample_clade.clade_assignments[c].end());
        }

        bool is_first=true;
        auto sample_node=tree.get_node(sample_idx);
        for (const auto& mut : sample_node->mutations) {
            if (__builtin_popcount(mut.get_mut_one_hot())!=1) {
                if(is_first){
                    fprintf (stderr, "Imputed mutations:\t");
                    is_first=false;
                }else {
                    fputc(';', stderr);
                    if(placement_stats_file){
                        fputc(';', placement_stats_file);
                    }
                }
                auto imputed_nuc=mut.get_par_one_hot()&mut.get_mut_one_hot();
                if (!imputed_nuc) {
                    imputed_nuc=1<<(__builtin_ctz(mut.get_mut_one_hot()));
                }
                fprintf (stderr, "%i:%c", mut.get_position(), MAT::get_nuc(imputed_nuc));
                if (placement_stats_file) {
                    fprintf (placement_stats_file, "%i:%c;", mut.get_position(), MAT::get_nuc(imputed_nuc));
                }

            }
        }
        if(!is_first){
            fprintf(stderr, "\n");
        }
        if (placement_stats_file) {
            fputc('\n', placement_stats_file);
        }
        delete in.placement_info;
        return;
    }
};
int count_mutation(const std::vector<To_Place_Sample_Mutation>& mutations){
    int mut_count=0;
    for (const auto& mut : mutations) {
        if (mut.mut_nuc!=0xf&&(!(mut.mut_nuc&mut.par_nuc))) {
            mut_count++;
        }
    }
    return mut_count;
};
typedef tbb::flow::function_node<Preped_Sample_To_Place *,size_t> placer_node_t;
typedef tbb::concurrent_bounded_queue<Sample_Muts*> retry_place_t;
typedef tbb::flow::multifunction_node<Sample_Muts*,tbb::flow::tuple<Sample_Muts*>> Pusher_Node_T;
typedef tbb::flow::function_node<print_format> Printer_Node_t;
static void place_sample_thread( MAT::Tree &main_tree,std::vector<MAT::Node *> &deleted_nodes,
    Placed_move_sended_state& send_queue,found_place_t& found_place_queue,Pusher_Node_T& retry_queue
    ,int all_size,std::atomic_bool& stop,std::atomic_size_t& curr_idx, 
    const int parsimony_increase_threshold, size_t sample_start_idx,bool dry_run,
    int max_parsimony,size_t max_uncertainty,bool multi_processing,Printer_Node_t& printer_node){
    int redo=0;
    int total=0;
    int parsimony_increase=0;
    int start_idx=curr_idx;
    int stop_count=all_size-start_idx;
    while (total<stop_count) {
        move_type* in;
        found_place_queue.pop(in);
        auto& search_result=std::get<0>(*in);
        bool skip=false;
        if (dry_run) {
            std::get<1>(*in)->sorting_key1=count_mutation(search_result[0].sample_mutations);
            std::get<1>(*in)->sorting_key2=search_result.size();
            fprintf(stderr, "Sample %zu, mutation count %d,curr_count %d, total %d\n",std::get<1>(*in)->sample_idx,std::get<1>(*in)->sorting_key1,total,stop_count);
            if (std::get<2>(*in)) {
                retry_queue.try_put(nullptr);
            }
            total++;
            delete in;
            continue;
        }
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
                retry_queue.try_put(std::get<1>(*in));
                delete in;
                /*if (!out->is_self) {
                    fprintf(stderr, "backlog Prep%zu \n",backlog_prep--);            
                }*/
                redo++;
                fprintf(stderr, "%f%% redone, %d redo, %d total\n",(float)redo/(float)total,redo,total);
                skip=true;
                break;
            }
        }
        if (skip) {
            continue;
        }
        total++;
        auto mut_size=count_mutation(search_result[0].sample_mutations);
        if (mut_size>max_parsimony||search_result.size()>max_uncertainty) {
            if (std::get<2>(*in)) {
                retry_queue.try_put(nullptr);
            }
            printer_node.try_put(print_format{mut_size,in});
            continue;
        }
        auto smallest_idx = search_result[0].target_node->node_id;
        auto arg_small_idx = 0;
        for (size_t target_idx = 1; target_idx < search_result.size();
             target_idx++) {
            if (smallest_idx > search_result[0].target_node->node_id) {
                smallest_idx = search_result[target_idx].target_node->node_id;
                arg_small_idx = target_idx;
            }
        }
        auto target = search_result[arg_small_idx];
        MAT::Mutations_Collection sample_mutations;
        discretize_mutations(target.sample_mutations, target.shared_mutations,
                             target.parent_node, sample_mutations);
        Preped_Sample_To_Place* to_ser;
        if (multi_processing) {
            to_ser=new Preped_Sample_To_Place;
            to_ser->target_id=target.target_node->node_id;
        }
        auto out = update_main_tree(sample_mutations, target.splited_mutations,
                                    target.shared_mutations, target.target_node,
                                    std::get<1>(*in)->sample_idx, main_tree, 0);
        if (out.deleted_nodes) {
            deleted_nodes.push_back(out.deleted_nodes);
            /*auto res=deleted_nodes.insert(out.deleted_nodes);
            if (!res.second) {
                fprintf(stderr, "repeated deleteion\n");
                raise(SIGTRAP);
            }*/
        }
        parsimony_increase+=mut_size;
        if (parsimony_increase>parsimony_increase_threshold) {
            stop.store(true);
        }
        //fprintf(stderr,"%zu:%s\t%d\t%zu\n",std::get<1>(*in)->sample_idx-sample_start_idx,main_tree.get_node_name(std::get<1>(*in)->sample_idx).c_str(), mut_size,search_result.size());
        if(multi_processing){
        to_ser->split_id = out.splitted_node ? out.splitted_node->node_id : 0;
        to_ser->splited_mutations=std::move(target.splited_mutations);
        to_ser->shared_mutations=std::move(target.shared_mutations);
        to_ser->sample_mutations=std::move(sample_mutations);
        to_ser->sample_id=std::get<1>(*in)->sample_idx;
        send_queue.push(to_ser);
        }
        //check_parent(main_tree.root,main_tree); 
        if (std::get<2>(*in)&&(parsimony_increase<=parsimony_increase_threshold)) {
            //fprintf(stderr, "self out\n");
            retry_queue.try_put(nullptr);
        }
        if (parsimony_increase>parsimony_increase_threshold) {
            stop_count=curr_idx-start_idx;
            fprintf(stderr, "curr_idx: %zu,stoped parsimpny score %d\n",curr_idx.load(),parsimony_increase);
        }
        printer_node.try_put(print_format{mut_size,in});
    }
}
struct Dist_sample_state{
    std::atomic_uint64_t& curr_idx;
    std::vector<Sample_Muts>& to_place;
    int processes_left;
    std::atomic_bool& stop;
};
static Status main_tree_distribute_samples(Dist_sample_state& state){
    int ignore;
        MPI_Status status;
        int recieved;
        MPI_Iprobe(MPI_ANY_SOURCE, PLACEMENT_WORK_REQ_TAG, MPI_COMM_WORLD, &recieved, &status);
        if (!recieved) {
            return NOTHING;
        }
        MPI_Recv(&ignore, 0, MPI_BYTE, MPI_ANY_SOURCE, PLACEMENT_WORK_REQ_TAG, MPI_COMM_WORLD, &status);
        mpi_trace_print("main recieved work req \n");
        int reply_dest=status.MPI_SOURCE;
        Mutation_Detailed::sample_to_place temp;
        size_t send_idx=SIZE_MAX;
        if (!state.stop) {
            send_idx=state.curr_idx++;    
        }
        //fprintf(stderr,"cur send idx %zu \n",send_idx);
        if (send_idx>=state.to_place.size()) { 
            MPI_Send(&ignore, 0, MPI_BYTE, reply_dest, PLACEMENT_WORK_RES_TAG, MPI_COMM_WORLD);
            state.processes_left--;
            if (state.processes_left) {
                return OK;
            }else {
                return DONE;
            }
        }
        auto& to_send=state.to_place[send_idx];
        temp.set_sample_id(to_send.sample_idx);
        fill_mutation_vect(temp.mutable_sample_mutation_positions(), temp.mutable_sample_mutation_other_fields(), to_send.muts);
        auto buffer=temp.SerializeAsString();
        mpi_trace_print( "main sending work res \n");
        MPI_Send(buffer.c_str(), buffer.size(), MPI_BYTE, reply_dest, PLACEMENT_WORK_RES_TAG, MPI_COMM_WORLD);
        mpi_trace_print("main sent work res \n");
        return OK;
}
static void mpi_loop(Dist_sample_state dist_sample,Recieve_Place_State recieve_place_state,Placed_move_sended_state& send_move_state,size_t sample_start_idx){
    bool dist_sample_done=false;
    bool recieve_place_done=false;
    bool send_move_done=false;
    while (!(dist_sample_done&recieve_place_done&send_move_done)) {
        bool not_idle=false;
        if (!dist_sample_done) {
            auto ret=main_tree_distribute_samples(dist_sample);
            if (ret==OK) {
                not_idle=true;
            }else if (ret==DONE) {
                dist_sample_done=true;
                fprintf(stderr, "Distributor exit \n");
            }
        }
        if(!recieve_place_done){
            auto ret=recieve_place(recieve_place_state,sample_start_idx);
            if (ret==OK) {
                not_idle=true;
            }else if (ret==DONE) {
                recieve_place_done=true;
                fprintf(stderr, "reciever exit \n");
            }
        }
        if (!send_move_done) {
            auto ret=send_paced_move(send_move_state);
            if (ret==OK) {
                not_idle=true;
            }else if (ret==DONE) {
                send_move_done=true;
            }
        }
        if (!not_idle) {
            sched_yield();
        }
    }
}
struct Pusher{
    std::atomic_uint64_t& curr_idx;
    std::vector<Sample_Muts>& to_place;
    std::atomic_bool& stop;
    void operator()(Sample_Muts* in,Pusher_Node_T::output_ports_type& out )const {
        if (!in) {
            size_t send_idx=SIZE_MAX;
            if (!stop) {
                send_idx=curr_idx++;        
            }
            if (send_idx<to_place.size()) {
                std::get<0>(out).try_put(&to_place[send_idx]);
                //fprintf(stderr, "main place %zu\n",send_idx);
            }
            return;
        }
        std::get<0>(out).try_put(in);
    }
};
struct Finder{
    MAT::Tree& tree;
    found_place_t& found_queue;
    void operator()(Sample_Muts* to_search) const{
        auto res=find_place(tree, to_search);
        std::get<2>(*res)=true;
            //fprintf(stderr, "self in\n");
        found_queue.push(res);
    }
};
void place_sample_leader(std::vector<Sample_Muts> &sample_to_place,
                         MAT::Tree &main_tree, int batch_size,
                         std::atomic_size_t &curr_idx,
                         int parsimony_increase_threshold,bool dry_run,
                         FILE *placement_stats_file,
                         int max_parsimony,size_t max_uncertainty,
                         std::vector<std::string>& low_confidence_samples,
                         std::vector<Clade_info>& samples_clade,
                         size_t sample_start_idx,std::vector<size_t>* idx_map
                         ) {
    std::vector<MAT::Node *> deleted_nodes;
    size_t node_count=main_tree.depth_first_expansion().size();
    deleted_nodes.reserve(sample_to_place.size());
    std::thread *mpi_thread=nullptr;
    tbb::concurrent_bounded_queue<Preped_Sample_To_Place*> send_queue;
    {   
        std::atomic_bool stop(false);
        tbb::flow::graph g;
        found_place_t found_queue;
        Pusher_Node_T init(g, tbb::flow::serial, Pusher{curr_idx, sample_to_place,stop});
        tbb::flow::function_node<
            Sample_Muts*>
            searcher(g, tbb::flow::unlimited,
                     Finder{main_tree,found_queue});
        tbb::flow::make_edge(std::get<0>(init.output_ports()),searcher);
        Printer_Node_t printer_node(g,tbb::flow::serial,
            Print_Thread{main_tree,placement_stats_file,max_parsimony,
            max_uncertainty,node_count,low_confidence_samples,samples_clade,sample_start_idx});
        if (process_count>1) {
            mpi_thread=new std::thread(mpi_loop,
            Dist_sample_state{curr_idx,sample_to_place,process_count-1,stop},
            Recieve_Place_State {found_queue,main_tree,process_count-1,sample_to_place,idx_map},
            std::ref(send_queue),sample_start_idx);
        }
        for(int temp=0;temp<batch_size;temp++){
            init.try_put(nullptr);
        }
        if (dry_run) {
            send_queue.try_push(nullptr);
        }
        place_sample_thread(main_tree, deleted_nodes, send_queue, found_queue,
                            init, sample_to_place.size(), stop, curr_idx,
                            parsimony_increase_threshold,
                            sample_start_idx, dry_run,max_parsimony,max_uncertainty,process_count>1,printer_node);
        g.wait_for_all();
        //placer_thread.join();
    }
     	send_queue.push(nullptr);
	fprintf(stderr,"main done\n");
        for(auto node:deleted_nodes){
            delete node;
        }
        deleted_nodes.clear();
        if (process_count>1) {
            mpi_thread->join();
            delete mpi_thread;
            size_t size=0;
            MPI_Bcast((void *)&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        }       
	fprintf(stderr,"main exit\n");
}
