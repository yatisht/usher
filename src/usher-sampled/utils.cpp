#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <atomic>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <sys/wait.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for.h>
#include <tbb/pipeline.h>
#include <thread>
#include <unistd.h>
#include <unordered_set>
#include "tbb/parallel_for_each.h"
#include <utility>
#include <vector>
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "mutation_detailed.pb.h"
#include <fstream>
#include <sstream>
void Min_Back_Fitch_Sankoff(MAT::Node* root_node,const MAT::Mutation& mut_template,
                            std::vector<std::vector<MAT::Mutation>>& mutation_output,mutated_t& positions,size_t dfs_size);
int set_descendant_count(MAT::Node* root) {
    size_t child_count=0;
    for (auto child : root->children) {
        child_count+=set_descendant_count(child);
    }
    root->bfs_index=child_count;
    return child_count;
}

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
void assign_levels(MAT::Node* root) {
    if (root->parent) {
        root->level=root->parent->level+1;
    } else {
        root->level=0;
    }
    for(auto child:root->children) {
        assign_levels(child);
    }
}
static void add_neighbor(MAT::Node* center,MAT::Node* exclude, int radius_left, tbb::concurrent_unordered_map<size_t,int>& to_search_node_idx) {
    auto res=to_search_node_idx.emplace(center->dfs_index,radius_left);
    if (!res.second) {
        if (res.first->second>=radius_left) {
            return;
        } else {
            res.first->second=radius_left;
        }
    }
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
void find_moved_node_neighbors(int radius,size_t start_idx, MAT::Tree& tree, size_t cur_idx,std::vector<size_t>& node_to_search_idx) {
    tree.depth_first_expansion();
    tbb::concurrent_unordered_map<size_t,int> to_search_node_idx_dict;
    to_search_node_idx_dict.rehash(tree.get_size_upper());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,cur_idx),[&](tbb::blocked_range<size_t> r) {
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            auto node=tree.get_node(idx+start_idx);
            if (!node) {
                continue;
            }
            add_neighbor(tree.get_node(idx+start_idx), nullptr, radius, to_search_node_idx_dict);
        }
    });
    node_to_search_idx.reserve(to_search_node_idx_dict.size());
    for(const auto& target:to_search_node_idx_dict) {
        node_to_search_idx.push_back(target.first);
    }
    fprintf(stderr, "%zu nodes to search \n",node_to_search_idx.size());
}
static void send_positions(const std::vector<mutated_t>& to_send,int start_position,int end_position, int target_rank) {
    //Send start position first
    auto length=end_position-start_position;
    MPI_Send(&length, 1, MPI_INT, target_rank, POSITION_TAG, MPI_COMM_WORLD);
    MPI_Send(&start_position, 1, MPI_INT, target_rank, POSITION_TAG, MPI_COMM_WORLD);
    for (int idx=start_position; idx<end_position; idx++) {
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
int follower_recieve_positions( std::vector<mutated_t>& to_recieve) {
    int length=0;
    int start_position;
    while (length==0) {
        MPI_Recv(&length, 1, MPI_INT, 0, POSITION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        fprintf(stderr, "follower recieving %d positions\n",length);
    }
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
void get_pos_samples_old_tree(MAT::Tree& tree,std::vector<mutated_t>& output) {
    Original_State_t origin_states;
    check_samples(tree.root,origin_states,&tree);
    std::vector<tbb::concurrent_vector<std::pair<size_t,char>>> pos_mutated(MAT::Mutation::refs.size());
    tbb::parallel_for_each(origin_states.begin(),origin_states.end(),[&pos_mutated,&tree](const std::pair<size_t, Mutation_Set>& sample_mutations) {
        for (const MAT::Mutation &m : sample_mutations.second) {
            pos_mutated[m.get_position()].emplace_back(tree.get_node(sample_mutations.first)->node_id,m.get_all_major_allele());
        }
    });
    fprintf(stderr, "output size %zu, nuc size %zu, pos_mutated size %zu\n",output.size(),MAT::Mutation::refs.size(),pos_mutated.size());
    for (size_t idx=0; idx<MAT::Mutation::refs.size(); idx++) {
        output[idx].insert(output[idx].end(),pos_mutated[idx].begin(),pos_mutated[idx].end());
    }
    if (process_count>1) {
        distribute_positions(output);
    }
}
void distribute_positions(std::vector<mutated_t>& output) {
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
template<typename ACC_type,typename result_t>
static void acc_mutations(result_t& FS_result,ACC_type& accumulator) {
    int total_size=0;
    for (auto& one_t_result : FS_result) {
        for (size_t idx=0; idx<one_t_result.output.size(); idx++ ) {
            total_size+=one_t_result.output[idx].size();
        }
    }
    fprintf(stderr, "accumulated %d\n",total_size);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,FS_result.local().output.size()),[&FS_result, &accumulator](tbb::blocked_range<size_t> range) {
        for (auto& one_t_result : FS_result) {
            for (size_t idx=range.begin(); idx<range.end(); idx++ ) {
                accumulator(one_t_result.output[idx],idx);
            }
        }
    });
}
template<typename result_t>
struct Parse_result {
    result_t& FS_result;
    size_t nodes_size;
    void operator()(std::pair<char*,size_t> in) {
        Mutation_Detailed::mutation_collection parsed;
        parsed.ParseFromArray(in.first, in.second);
        auto& local_res=FS_result.local().output;
        //fprintf(stderr, "adding to node %zu",parsed.node_idx());
        if (nodes_size>local_res.size()) {
            local_res.resize(nodes_size);
            fprintf(stderr, "Node idx %ld, total length %zu \n",parsed.node_idx(),FS_result.local().output.size());
        }
        auto& to_add=local_res[parsed.node_idx()];
        for (int idx=0; idx<parsed.positions_size(); idx++) {
            MAT::Mutation temp;
            *((int*)&temp)=parsed.positions(idx);
            *((int*)&temp+1)=parsed.other_fields(idx);
            to_add.push_back(temp);
        }
        delete [] in.first;
    }
};
static void reassign_state_preprocessing(MAT::Tree& tree) {
    std::unordered_set<size_t> ignored;
    std::unordered_set<size_t> ignored2;
    clean_up_internal_nodes(tree.root,tree,ignored,ignored2);
    fprintf(stderr, "parsiomony score %zu\n",tree.get_parsimony_score());
    for (auto node : tree.depth_first_expansion()) {
        node->mutations.clear();
    }
}
static void reassign_state_kernel(MAT::Tree& tree,const std::vector<mutated_t>& mutations,int start_position,std::vector<MAT::Node*>& bfs_ordered_nodes,FS_result_per_thread_t& FS_result) {
    //get mutation vector
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
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
            translated.reserve(mutations[idx].size()+1);
            for (auto& node_p :mutations[idx] ) {
                auto node=tree.get_node(node_p.first);
                if (!node) {
                    //fprintf(stderr, "%zu node not found from rank %d\n",node_p.first,this_rank);
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
}

struct mutated_t_inorder_comparator {
    bool operator()(const std::pair<long, nuc_one_hot>& lhs,const std::pair<long, nuc_one_hot>& rhs) const {
        return lhs.first<rhs.first;
    }
};
struct Min_Back_FS_Result_Container {
    std::vector<std::vector<Mutation_Annotated_Tree::Mutation>> output;
};
typedef tbb::enumerable_thread_specific<Min_Back_FS_Result_Container> Min_Back_FS_result;
static void min_backreassign_state_kernel(MAT::Tree& tree,const std::vector<mutated_t>& mutations,int start_position,std::vector<MAT::Node*>& dfs,Min_Back_FS_result& FS_result) {
    //get mutation vector
    fprintf(stderr, "rank %d min back assigning %zu nuc\n",this_rank,mutations.size());
    auto dfs_size=dfs.size();
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0,mutations.size()),
    [&FS_result,&mutations,&tree,start_position,dfs_size](const tbb::blocked_range<size_t>& in) {
        auto& this_result=FS_result.local().output;
        this_result.resize(dfs_size);
        for (size_t idx=in.begin(); idx<in.end(); idx++) {
            if (mutations[idx].empty()) {
                //fprintf(stderr, "rank %d skipping empty position %zu \n",this_rank,idx);
                continue;
            }
            mutated_t translated;
            translated.reserve(mutations[idx].size()+1);
            for (auto& node_p :mutations[idx] ) {
                auto node=tree.get_node(node_p.first);
                if (!node) {
                    //fprintf(stderr, "%zu node not found from rank %d\n",node_p.first,this_rank);
                    continue;
                }
                translated.emplace_back(node->dfs_index,node_p.second);
            }
            std::sort(translated.begin(),translated.end(),mutated_t_inorder_comparator());
            translated.emplace_back(INT_MAX,0xf);
            MAT::Mutation temp(0,idx+start_position,0,0);
            Min_Back_Fitch_Sankoff(tree.root, temp, this_result, translated, dfs_size);
        }
    });
}


template<typename result_t>
static void output_mutations(std::vector<MAT::Node*>& bfs_ordered_nodes,result_t& FS_result) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0,bfs_ordered_nodes.size()),[&FS_result,&bfs_ordered_nodes](tbb::blocked_range<size_t> range) {
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
}
void reassign_state_local(MAT::Tree& tree,const std::vector<mutated_t>& mutations,bool initial) {
    if (!initial) {
        reassign_state_preprocessing(tree);
    }
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    FS_result_per_thread_t FS_result;
    reassign_state_kernel(tree, mutations, 0, bfs_ordered_nodes, FS_result);
    output_mutations(bfs_ordered_nodes, FS_result);
}
void min_back_reassign_state_local(MAT::Tree& tree,const std::vector<mutated_t>& mutations) {
    reassign_state_preprocessing(tree);
    auto dfs_ordered_nodes = tree.depth_first_expansion();
    Min_Back_FS_result FS_result;
    min_backreassign_state_kernel(tree, mutations, 0, dfs_ordered_nodes, FS_result);
    output_mutations(dfs_ordered_nodes, FS_result);
}

template <typename result_t>
void FS_gather_mut(MAT::Tree &tree,
                   std::vector<Mutation_Annotated_Tree::Node *> &bfs_ordered_nodes,
                   result_t &FS_result) {
    if (this_rank == 0) {
        tbb::flow::graph g;
        int processes_left = process_count - 1;
        tbb::flow::function_node<std::pair<char *, size_t>> proc(
                    g, tbb::flow::unlimited,
                    Parse_result<result_t> {FS_result, bfs_ordered_nodes.size()});
        bool is_first = true;
        char ignore;
        fprintf(stderr, "Start recieving assignment message\n");
        size_t before_recieving_count = 0;
        for (const auto &res : FS_result) {
            for (const auto &temp : res.output) {
                before_recieving_count += temp.size();
            }
        }
        fprintf(stderr, "%zu mutations before recieving\n",
                before_recieving_count);
        while (processes_left) {
            int count;
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, FS_RESULT_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &count);
            // fprintf(stderr, "recieved message of size %d from %d, status
            // %d\n",count,status.MPI_SOURCE,status.MPI_ERROR);
            if (count == 0) {
                if (!is_first) {
                    processes_left--;
                }
                MPI_Recv(&ignore, count, MPI_BYTE, status.MPI_SOURCE,
                         FS_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                continue;
            }
            is_first = false;
            auto buffer = new char[count];
            MPI_Recv(buffer, count, MPI_BYTE, status.MPI_SOURCE, FS_RESULT_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            proc.try_put(std::make_pair(buffer, count));
        }
        fprintf(stderr, "Recieved last assignment message\n");
        g.wait_for_all();
        fprintf(stderr, "Finished parsing assignment messages\n");
        size_t after_recieving_count = 0;
        for (const auto &res : FS_result) {
            for (const auto &temp : res.output) {
                after_recieving_count += temp.size();
            }
        }
        fprintf(stderr, "%zu mutations after recieving\n",
                after_recieving_count);
        output_mutations(bfs_ordered_nodes, FS_result);
        fprintf(stderr, "Finished filling mutations\n");
        tree.populate_ignored_range();
        fprintf(stderr, "parsiomony score %zu\n", tree.get_parsimony_score());
        tree.MPI_send_tree();
    } else {
        std::vector<Mutation_Detailed::mutation_collection>
        mutation_collection_to_serialize(bfs_ordered_nodes.size());
        struct {
            std::vector<Mutation_Detailed::mutation_collection>
            &mutation_collection_to_serialize;
            void operator()(std::vector<MAT::Mutation> &to_output,
                            size_t idx) {
                for (const auto &mut : to_output) {
                    mutation_collection_to_serialize[idx].add_positions(
                        mut.get_position());
                    mutation_collection_to_serialize[idx].add_other_fields(
                        *((int *)&mut + 1));
                }
            }
        } other_T_acc{mutation_collection_to_serialize};
        acc_mutations(FS_result, other_T_acc);
        fprintf(stderr, "Start sending assignment messages\n");
        for (size_t idx = 0; idx < bfs_ordered_nodes.size(); idx++) {
            if (mutation_collection_to_serialize[idx].positions_size() != 0) {
                mutation_collection_to_serialize[idx].set_node_idx(idx);
                auto serialized =
                    mutation_collection_to_serialize[idx].SerializeAsString();
                // fprintf(stderr, "sending assigned of size
                // %zu\n",serialized.size());
                if (serialized.size() == 0) {
                    fprintf(stderr, "Sending assignment message of size 0\n");
                    raise(SIGTRAP);
                }
                MPI_Send(serialized.c_str(), serialized.size(), MPI_BYTE, 0,
                         FS_RESULT_TAG, MPI_COMM_WORLD);
            }
        }
        MPI_Send(&this_rank, 0, MPI_BYTE, 0, FS_RESULT_TAG, MPI_COMM_WORLD);
        fprintf(stderr, "Finish sending assignment messages\n");
        tree.delete_nodes();
        tree.MPI_receive_tree();
    }
}
void FS_scatter_tree(MAT::Tree &tree, bool initial) {
    if (this_rank == 0) {
        if (!initial) {
            reassign_state_preprocessing(tree);
        }
        tree.MPI_send_tree();
    } else {
        tree.delete_nodes();
        tree.MPI_receive_tree();
    }
}
void MPI_reassign_states(MAT::Tree &tree,
                         const std::vector<mutated_t> &mutations,
                         int start_position, bool initial) {
    FS_scatter_tree(tree, initial);
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    FS_result_per_thread_t FS_result;
    reassign_state_kernel(tree, mutations, start_position, bfs_ordered_nodes,
                          FS_result);
    FS_gather_mut(tree, bfs_ordered_nodes, FS_result);
}

void MPI_min_back_reassign_states(MAT::Tree &tree,
                                  const std::vector<mutated_t> &mutations,
                                  int start_position) {
    FS_scatter_tree(tree, false);
    auto dfs_ordered_nodes = tree.depth_first_expansion();
    Min_Back_FS_result FS_result;
    min_backreassign_state_kernel(tree, mutations, start_position, dfs_ordered_nodes,
                                  FS_result);
    FS_gather_mut(tree, dfs_ordered_nodes, FS_result);
}
static void remove_node_helper(MAT::Node* to_remove,MAT::Tree& tree) {
    tree.erase_node(to_remove->node_id);
    auto& par_children=to_remove->parent->children;
    auto iter=find(par_children.begin(),par_children.end(),to_remove);
    par_children.erase(iter);

}
static void remove_absent_leaves(MAT::Node* node,MAT::Tree& tree,std::unordered_set<std::string>& present) {
    if (node->children.empty()) {
        auto this_name=tree.get_node_name(node->node_id);
        auto iter=present.find(this_name);
        if (iter==present.end()&&(!node->is_root())) {
            fprintf(stderr,"sample %s absent in VCF,removed\n", this_name.c_str());
            //absent
            remove_node_helper(node, tree);
            delete node;
        }
        return;
    }
    auto old_children=node->children;
    for (auto child : old_children) {
        remove_absent_leaves(child, tree, present);
    }
    if (node->is_root()) {
        return;
    }
    if (node->children.empty()) {
        remove_node_helper(node, tree);
        delete node;
    }
    if (node->children.size()==1) {
        remove_node_helper(node, tree);
        node->parent->children.push_back(node->children[0]);
        node->children[0]->parent=node->parent;
        delete node;
    }
}

void remove_absent_leaves(MAT::Tree& tree,std::unordered_set<std::string>& present) {
    remove_absent_leaves(tree.root,tree,present);
}
static int output_newick(MAT::Tree& T,const output_options& options,int t_idx) {
    std::string uncondensed_string=options.print_uncondensed_tree?"uncondensed":"";
    auto final_tree_filename = options.outdir + "/";
    if (options.print_uncondensed_tree) {
        T.uncondense_leaves();
        final_tree_filename+="uncondensed-";
    }
    final_tree_filename+="final-tree";
    if (options.only_one_tree) {
        fprintf(stderr, "Writing uncondensed final tree to file %s \n", final_tree_filename.c_str());
        final_tree_filename+=".nh";
    } else {
        final_tree_filename += std::to_string(t_idx+1) + ".nh";
        fprintf(stderr, "Writing uncondensed final tree %d to file %s \n", (t_idx+1), final_tree_filename.c_str());
    }
    //auto parsimony_score = T.get_parsimony_score();
    //fprintf(stderr, "The parsimony score for this tree is: %zu \n", parsimony_score);
    int pid=fork();
    if (pid==0) {
        std::ofstream final_tree_file(final_tree_filename.c_str(), std::ofstream::out);
        T.write_newick_string(final_tree_file,T.root, true, true, options.retain_original_branch_len,true);
        final_tree_file.close();
        return 0;
    } else {
        return pid;
    }

}
void check_leaves(const MAT::Tree& T) {
    for (const auto node :T.depth_first_expansion() ) {
        if (node!=T.get_node(node->node_id)) {
            fprintf(stderr, "Node mismatch\n");
            raise(SIGTRAP);
        }
        if (node->children.empty()) {
            if (T.get_node_name(node->node_id)=="") {
                fprintf(stderr,"Leaf node without name %zu",node->node_id);
                raise(SIGTRAP);
            }
        }
    }
}
void print_annotation(const MAT::Tree &T, const output_options &options,
                      const std::vector<Clade_info> &assigned_clades,
                      size_t sample_start_idx, size_t sample_end_idx,
                      size_t num_annotations) {
    auto annotations_filename = options.outdir + "/clades.txt";

    FILE *annotations_file = fopen(annotations_filename.c_str(), "w");

    for (size_t s = 0; s < sample_end_idx - sample_start_idx; s++) {
        if (!assigned_clades[s].valid) {
            continue;
        }
        if (assigned_clades[s].best_clade_assignment.size() == 0) {
            // Sample was not placed (e.g. exceeded max EPPs) so no clades
            // assigned
            continue;
        }
        auto sample = T.get_node_name(sample_start_idx + s);

        fprintf(annotations_file, "%s\t", sample.c_str());
        for (size_t k = 0; k < num_annotations; k++) {
            fprintf(annotations_file, "%s",
                    assigned_clades[s].best_clade_assignment[k].c_str());
            // TODO
            fprintf(annotations_file, "*|");
            if (options.detailed_clades) {
                std::string curr_clade = "";
                int curr_count = 0;
                for (auto clade : assigned_clades[s].clade_assignments[k]) {
                    if (clade == curr_clade) {
                        curr_count++;
                    } else {
                        if (curr_count > 0) {
                            fprintf(
                                annotations_file, "%s(%i/%zu),",
                                curr_clade.c_str(), curr_count,
                                assigned_clades[s].clade_assignments[k].size());
                        }
                        curr_clade = clade;
                        curr_count = 1;
                    }
                }
                if (curr_count > 0) {
                    fprintf(annotations_file, "%s(%i/%zu)", curr_clade.c_str(),
                            curr_count,
                            assigned_clades[s].clade_assignments[k].size());
                }
            }
            if (k + 1 < num_annotations) {
                fprintf(annotations_file, "\t");
            }
        }
        fprintf(annotations_file, "\n");
    }

    fclose(annotations_file);
}
bool final_output(MAT::Tree &T, const output_options &options, int t_idx,
                  std::vector<Clade_info> &assigned_clades,
                  size_t sample_start_idx, size_t sample_end_idx,
                  std::vector<std::string> &low_confidence_samples,std::vector<mutated_t>& position_wise_out) {
    // If user need uncondensed tree output, write uncondensed tree(s) to
    // file(s)
    //check_leaves(T);
    {
        if (options.redo_FS_Min_Back_Mutations) {
            fprintf(stderr, "Parsimony score before %zu\n",T.get_parsimony_score());
            fprintf(stderr, "Back mutation count before %d\n",count_back_mutation(T));
            if (process_count>1) {
                MPI_min_back_reassign_states(T, position_wise_out, 0);
            } else {
                min_back_reassign_state_local(T, position_wise_out);
            }
            fprintf(stderr, "Parsimony score after %zu\n",T.get_parsimony_score());
            fprintf(stderr, "Back mutation count after %d\n",count_back_mutation(T));
        }
        std::unordered_set<size_t> ignored1;
        std::unordered_set<size_t> ignored2;
        clean_up_internal_nodes(T.root,T,ignored1,ignored2);
    }
    fix_condensed_nodes(&T);
    //check_leaves(T);
    //MAT::save_mutation_annotated_tree(T, "before_post_processing.pb");
    MPI_Finalize();
    int pid=output_newick(T, options,t_idx);
    if(!pid){
        _exit(EXIT_SUCCESS);
    }
    // For each final tree write the path of mutations from tree root to the
    // sample for each newly placed sample
    bool use_tree_idx = !options.only_one_tree;
    std::vector<MAT::Node *> targets;
    fprintf(stderr, "Sample start idx %zu, end index %zu\n",sample_start_idx,sample_end_idx);
    targets.reserve(sample_end_idx - sample_start_idx);

    for (size_t idx = sample_start_idx; idx < sample_end_idx; idx++) {
        auto node = T.get_node(idx);
        if (node) {
            targets.push_back(node);
        }
    }
    auto mutation_paths_filename = options.outdir + "/mutation-paths.txt";
    if (use_tree_idx) {
        mutation_paths_filename = options.outdir + "/mutation-paths-" +
                                  std::to_string(t_idx + 1) + ".txt";
        fprintf(stderr, "Writing mutation paths for tree %d to file %s \n",
                t_idx + 1, mutation_paths_filename.c_str());
    } else {
        fprintf(stderr, "Writing mutation paths to file %s \n",
                mutation_paths_filename.c_str());
    }
    MAT::get_sample_mutation_paths(&T, targets, mutation_paths_filename);
    // fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    // For each final tree write the annotations for each sample

    size_t num_annotations = T.get_num_annotations();
    //check_leaves(T);

    if (num_annotations > 0&&options.only_one_tree) {
        // timer.Start();

        print_annotation(T, options, assigned_clades, sample_start_idx, sample_end_idx,
                         num_annotations);

        // fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //check_leaves(T);
    if ((options.print_subtrees_single > 1) ) {
        fprintf(stderr, "Computing the single subtree for added samples with %zu random leaves. \n\n", options.print_subtrees_single);
        //timer.Start();
        // For each final tree, write a subtree of user-specified size around
        // each newly placed sample in newick format
        MAT::get_random_single_subtree(T, targets, options.outdir, options.print_subtrees_single, t_idx, use_tree_idx, options.retain_original_branch_len);
        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //check_leaves(T);
    if ((options.print_subtrees_size > 1)) {
        fprintf(stderr, "Computing subtrees for added samples. \n\n");

        // For each final tree, write a subtree of user-specified size around
        // each newly placed sample in newick format
        //timer.Start();
        MAT::get_random_sample_subtrees(T, targets, options.outdir, options.print_subtrees_size, t_idx, use_tree_idx, options.retain_original_branch_len);
        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    // Print warning message with a list of all samples placed with low
    // confidence (>=2 parsimony-optimal placements)
    if (low_confidence_samples.size() > 0) {
        fprintf(stderr, "WARNING: Following samples had multiple possibilities of parsimony-optimal placements:\n");
        for (auto lcs: low_confidence_samples) {
            fprintf(stderr, "%s\n", lcs.c_str());
        }
    }

    // Store mutation-annotated tree to a protobuf file if user has asked for it
    if (options.dout_filename != "") {

        //timer.Start();
        auto this_out_name=options.dout_filename;
        if (!options.only_one_tree) {
            this_out_name+="."+std::to_string(t_idx);
        }
        fprintf(stderr, "Saving mutation-annotated tree object to file (after condensing identical sequences) %s\n", options.dout_filename.c_str());
        // Recondense tree with new samples
        check_leaves(T);
        if (T.condensed_nodes.size() > 0) {
            T.uncondense_leaves();
        }
        check_leaves(T);
        T.condense_leaves();
        MAT::save_mutation_annotated_tree(T, options.dout_filename);

        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    wait(&pid);
    return true;
}
/*static void check_repeats(const std::vector<Sample_Muts>& samples_to_place,size_t sample_start_idx){
    std::vector<bool> encountered(samples_to_place.size());
    for(const auto& to_place:samples_to_place){
        auto array_idx=to_place.sample_idx-sample_start_idx;
        if (encountered[array_idx]) {
            fprintf(stderr, "Repeated\n");
        }
        encountered[array_idx]=true;
    }
    for (size_t idx=0; idx<samples_to_place.size(); idx++) {
        if (!encountered[idx]) {
            fprintf(stderr, "%zu missing\n",idx+sample_start_idx);
            raise(SIGTRAP);
        }
    }
}*/
