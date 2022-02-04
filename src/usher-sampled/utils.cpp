#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <mpi.h>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for.h>
#include <tbb/pipeline.h>
#include <thread>
#include <unordered_set>
#include "tbb/parallel_for_each.h"
#include <utility>
#include <vector>
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "mutation_detailed.pb.h"
void convert_mut_type(const std::vector<MAT::Mutation> &in,
                      std::vector<To_Place_Sample_Mutation> &out) {
    out.reserve(in.size());
    for (const auto &mut : in) {
        if (mut.get_mut_one_hot() != 0xf) {
            out.push_back(To_Place_Sample_Mutation(mut.get_position(),mut.get_chromIdx(),mut.get_mut_one_hot(),mut.get_par_one_hot()));
        } else {
            if (out.empty() || // first
                out.back().mut_nuc!=0xf|| //last mutation not N
                (out.back().get_end_range() + 1) !=
                    mut.get_position() ||           // not contiguous
                out.back().range == UINT8_MAX // overflow
            ) {
                out.push_back(To_Place_Sample_Mutation(
                    mut.get_position(), mut.get_chromIdx(), mut.get_mut_one_hot()));
            } else {
                out.back().range++;
            }
        }
    }
}
void assign_levels(MAT::Node* root){
    if (root->parent) {
        root->level=root->parent->level+1;
    }else {
        root->level=0;
    }
    for(auto child:root->children){
        assign_levels(child);
    }
}
static void add_neighbor(MAT::Node* center,MAT::Node* exclude, int radius_left, std::unordered_set<size_t>& to_search_node_idx){
    to_search_node_idx.insert(center->dfs_index);
    if (radius_left<0) {
        return;
    }
    if (center->parent&&center->parent!=exclude) {
        add_neighbor(center->parent, center, radius_left-1, to_search_node_idx);
    }
    for (auto child : center->children) {
        add_neighbor(child, center, radius_left-1, to_search_node_idx);
    }
}
void find_moved_node_neighbors(int radius,size_t start_idx, MAT::Tree& tree, size_t cur_idx,std::vector<size_t>& node_to_search_idx){
    tree.depth_first_expansion();
    std::unordered_set<size_t> to_search_node_idx_dict;
    for (size_t ori_idx=0; ori_idx<cur_idx;ori_idx++ ) {
        add_neighbor(tree.get_node(ori_idx+start_idx), nullptr, radius, to_search_node_idx_dict);
    }
    node_to_search_idx.insert(node_to_search_idx.begin(),to_search_node_idx_dict.begin(),to_search_node_idx_dict.end());
    fprintf(stderr, "%zu nodes to search \n",node_to_search_idx.size());
}
static void send_positions(const std::vector<mutated_t>& to_send,int start_position,int end_position, int target_rank){
    //Send start position first
    auto length=end_position-start_position;
    MPI_Send(&length, 1, MPI_INT, target_rank, POSITION_TAG, MPI_COMM_WORLD);
    MPI_Send(&start_position, 1, MPI_INT, target_rank, POSITION_TAG, MPI_COMM_WORLD);
    for (size_t idx=start_position; idx<end_position; idx++) {
        Mutation_Detailed::mutation_at_each_pos to_serailize;
        to_serailize.mutable_mut()->Reserve(to_send[idx].size());
        to_serailize.mutable_node_id()->Reserve(to_send[idx].size());
        for (const auto& temp : to_send[idx]) {
            to_serailize.add_node_id(temp.first);
            to_serailize.add_mut(temp.second);
        }
        auto to_send_serialized=to_serailize.SerializeAsString();
        MPI_Send(to_send_serialized.data(), to_send_serialized.size(), MPI_BYTE, target_rank, POSITION_TAG,MPI_COMM_WORLD);
    }
    MPI_Send(&start_position, 0, MPI_INT, target_rank, POSITION_TAG, MPI_COMM_WORLD);
}
int follower_recieve_positions( std::vector<mutated_t>& to_recieve){
    int length;
    int start_position;
    MPI_Recv(&length, 1, MPI_INT, 0, POSITION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    MPI_Recv(&start_position, 1, MPI_INT, 0, POSITION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    to_recieve.reserve(length);
    while (length) {
        MPI_Status status;
        MPI_Probe(0, POSITION_TAG, MPI_COMM_WORLD, &status);
        int count;
        MPI_Get_count(&status, MPI_BYTE, &count);
        auto buffer = new char[count];
        MPI_Recv(buffer, count, MPI_INT, 0, POSITION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        Mutation_Detailed::mutation_at_each_pos parsed;
        parsed.ParseFromArray(buffer, count);
        to_recieve.emplace_back();
        auto& last=to_recieve.back();
        auto total_count=parsed.node_id_size();
        if (total_count==0) {
            //fprintf(stderr, "rank %d recieving empty position %d\n",this_rank,to_recieve.size());
        }
        last.reserve(total_count);
        for (int idx=0; idx<total_count; idx++) {
            last.emplace_back(parsed.node_id(idx),parsed.mut(idx));
        }
        delete[] buffer;
        length--;
    }
    return start_position;
}
void get_pos_samples_old_tree(MAT::Tree& tree,std::vector<mutated_t>& output){
    Original_State_t origin_states;
    check_samples(tree.root,origin_states,&tree);
    std::vector<tbb::concurrent_vector<std::pair<size_t,char>>> pos_mutated(MAT::Mutation::refs.size());
    tbb::parallel_for_each(origin_states.begin(),origin_states.end(),[&pos_mutated,&tree](const std::pair<size_t, Mutation_Set>& sample_mutations) {
        for (const MAT::Mutation &m : sample_mutations.second) {
            pos_mutated[m.get_position()].emplace_back(tree.get_node(sample_mutations.first)->node_id,m.get_all_major_allele());
        }
    });
    fprintf(stderr, "output size %zu, nuc size %zu, pos_mutated size %zu\n",output.size(),MAT::Mutation::refs.size(),pos_mutated.size());
    for (auto idx=0; idx<MAT::Mutation::refs.size(); idx++) {
        output[idx].insert(output[idx].end(),pos_mutated[idx].begin(),pos_mutated[idx].end());
    }
    auto positions_per_process=(output.size()/process_count)+1;
    auto start_idx=output.size()-(process_count-1)*positions_per_process;
    std::vector<std::thread> threads;
    for (int proc_idx=0; proc_idx<process_count-1; proc_idx++) {
        threads.emplace_back(send_positions,std::ref(output), start_idx+proc_idx*positions_per_process, start_idx+(proc_idx+1)*positions_per_process, proc_idx+1);
    }
    for (auto& t : threads) {
        t.join();
    }
    output.resize(start_idx);
}
void clean_up_internal_nodes(MAT::Node* this_node,MAT::Tree& tree,std::unordered_set<size_t>& changed_nodes_local,std::unordered_set<size_t>& node_with_inconsistent_state);
template<typename ACC_type>
static void acc_mutations(FS_result_per_thread_t& FS_result,ACC_type& accumulator){
    int total_size=0;
    for (auto& one_t_result : FS_result) {
        for (size_t idx=0;idx<one_t_result.output.size();idx++ ) {
            total_size+=one_t_result.output[idx].size();
        }
    }
    fprintf(stderr, "accumulated %d\n",total_size);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,FS_result.local().output.size()),[&FS_result, &accumulator](tbb::blocked_range<size_t> range){
            for (auto& one_t_result : FS_result) {
                for (size_t idx=range.begin();idx<range.end();idx++ ) {
                    accumulator(one_t_result.output[idx],idx);
                }
            }
        });
}
struct Recieve_Feeder{
    int &thread_left;
    bool operator()(std::pair<char* ,size_t>& out) const{
        int count;
        MPI_Status status;
        while (true) {
            MPI_Probe(MPI_ANY_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT,&count);
            fprintf(stderr, "recieve message of size %d\n",count);
            if (count<0) {
                raise(SIGTRAP);
            }
            if (count==0) {
                thread_left--;
                if (!thread_left) {
                    return false;
                }
            }else {
                break;
            }
        }
        out.second=count;
        out.first=new char[count];
        MPI_Recv(out.first, count, MPI_BYTE, status.MPI_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return true;
    }
};

struct Parse_result{
    FS_result_per_thread_t& FS_result;
    void operator()(std::pair<char* ,size_t> in){
        Mutation_Detailed::mutation_collection parsed;
        parsed.ParseFromArray(in.first, in.second);
        //fprintf(stderr, "adding to node %zu",parsed.node_idx());
        auto& to_add=FS_result.local().output[parsed.node_idx()];
        for (size_t idx=0; idx<parsed.positions_size(); idx++) {
            MAT::Mutation temp;
            *((int*)&temp)=parsed.positions(idx);
            *((int*)&temp+1)=parsed.other_fields(idx);
            to_add.push_back(temp);
        }
        delete [] in.first;
    }
};

void MPI_reassign_states(MAT::Tree& tree,const std::vector<mutated_t>& mutations,int start_position){
    if (this_rank==0) {
        std::unordered_set<size_t> ignored;
        std::unordered_set<size_t> ignored2;
        clean_up_internal_nodes(tree.root,tree,ignored,ignored2);
        fprintf(stderr, "parsiomony score %zu\n",tree.get_parsimony_score());
        for (auto node : tree.depth_first_expansion()) {
            node->mutations.clear();
        }
        tree.MPI_send_tree();
    }else {
        tree.delete_nodes();
        tree.MPI_receive_tree();
    }
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    //get mutation vector
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    FS_result_per_thread_t FS_result;
    fprintf(stderr, "rand %d assigning %zu nuc\n",this_rank,mutations.size());
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0,mutations.size()),
    [&FS_result,&child_idx_range,&parent_idx,&mutations,&tree,start_position](const tbb::blocked_range<size_t>& in) {
        auto& this_result=FS_result.local();
        this_result.init(child_idx_range.size());
        for (size_t idx=in.begin(); idx<in.end(); idx++) {
            if (mutations[idx].empty()) {
                //fprintf(stderr, "rank %d skipping empty position %zu \n",this_rank,idx);
                continue;
            }
            mutated_t translated;
            translated.reserve(mutations[idx].size());
            for (auto& node_p :mutations[idx] ) {
                auto node=tree.get_node(node_p.first);
                if (!node) {
                    continue;
                }
                translated.emplace_back(node->bfs_index,node_p.second);
            }
            std::sort(translated.begin(),translated.end(),mutated_t_comparator());
            translated.emplace_back(0,0xf);
            MAT::Mutation temp(0,idx+start_position,0,0);
            Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx, temp,translated,
                                     this_result);

        }
    });
    if (this_rank==0) {
        tbb::flow::graph g;
        int processes_left=process_count-1;
        tbb::flow::function_node<std::pair<char*, size_t>> proc(g,tbb::flow::unlimited,Parse_result{FS_result});
        bool is_first=true;
        char ignore;
        while (processes_left) {
            int count;
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &count);
            //fprintf(stderr, "recieved message of size %d from %d, status %d\n",count,status.MPI_SOURCE,status.MPI_ERROR);
            if (count==0) {
                if (!is_first) {
                    processes_left--;
                }
                MPI_Recv(&ignore, count, MPI_BYTE, status.MPI_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                continue;
            }
            is_first=false;
            auto buffer=new char[count];
            MPI_Recv(buffer, count, MPI_BYTE, status.MPI_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            proc.try_put(std::make_pair(buffer, count));
        }
        g.wait_for_all();
        tbb::parallel_for(tbb::blocked_range<size_t>(0,bfs_ordered_nodes.size()),[&FS_result,&bfs_ordered_nodes](tbb::blocked_range<size_t> range){
            for (auto & temp : FS_result) {
                for (size_t idx=range.begin(); idx<range.end(); idx++) {
                    auto& to_fill=bfs_ordered_nodes[idx]->mutations.mutations;
                    to_fill.insert(to_fill.end(),temp.output[idx].begin(),temp.output[idx].end());
                }
            }
            for (size_t idx=range.begin(); idx<range.end(); idx++) {
                    auto& to_fill=bfs_ordered_nodes[idx]->mutations.mutations;
                    std::sort(to_fill.begin(),to_fill.end());
            }
        });
        tree.populate_ignored_range();
        fprintf(stderr, "parsiomony score %zu\n",tree.get_parsimony_score());
        tree.MPI_send_tree();
    }else {
        std::vector<Mutation_Detailed::mutation_collection> mutation_collection_to_serialize(bfs_ordered_nodes.size());
        struct {
            std::vector<Mutation_Detailed::mutation_collection>& mutation_collection_to_serialize;
            void operator()(std::vector<MAT::Mutation>& to_output, size_t idx){
                for (const auto& mut : to_output) {
                        mutation_collection_to_serialize[idx].add_positions(mut.get_position());
                        mutation_collection_to_serialize[idx].add_other_fields(*((int*)&mut+1));
                    }
            }
        } other_T_acc{mutation_collection_to_serialize};
        acc_mutations(FS_result, other_T_acc);
        for (size_t idx=0; idx<bfs_ordered_nodes.size(); idx++) {
            if (mutation_collection_to_serialize[idx].positions_size()!=0) {
                mutation_collection_to_serialize[idx].set_node_idx(idx);
                auto serialized=mutation_collection_to_serialize[idx].SerializeAsString();
                //fprintf(stderr, "sending assigned of size %zu\n",serialized.size());
                MPI_Send(serialized.c_str(), serialized.size(), MPI_BYTE, 0, FS_RESULT_TAG, MPI_COMM_WORLD);
            }
        }
        MPI_Send(&this_rank, 0, MPI_BYTE, 0, FS_RESULT_TAG, MPI_COMM_WORLD);
        tree.delete_nodes();
        tree.MPI_receive_tree();
    }
}