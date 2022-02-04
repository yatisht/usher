#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "usher.hpp"
#include <algorithm>
#include <atomic>
#include <boost/program_options.hpp>
#include <chrono>
#include <climits>
#include <complex>
#include <csignal>
#include <cstdio>
#include <mpi.h>
#include <random>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <utility>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
bool use_bound;
int process_count;
int this_rank;
unsigned int num_threads;
std::atomic_bool interrupted(false);
void fix_condensed_nodes(MAT::Tree *tree);
namespace po = boost::program_options;
static void clean_up_leaf(std::vector<MAT::Node*>& dfs){
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dfs.size()),[&dfs](tbb::blocked_range<size_t> range){
        for (size_t node_idx=range.begin(); node_idx<range.end(); node_idx++) {
            auto node=dfs[node_idx];
                /*if (node->identifier=="s1433144s") {
        raise(SIGTRAP);
    }
    if (node->identifier=="s2886812s") {
        raise(SIGTRAP);
    }
    if (node->identifier=="s2749940s") {
        raise(SIGTRAP);
    }*/
            if (node->is_leaf()) {
                auto& muts=node->mutations.mutations;
                muts.erase(std::remove_if(muts.begin(), muts.end(), [](const auto& mut){
                    return mut.get_par_one_hot()&mut.get_mut_one_hot();
                }),muts.end());
                for (auto &mut : muts) {
                    /*if (mut.get_par_one_hot()&mut.get_mut_one_hot()) {
                        raise(SIGTRAP);
                    }*/
                    mut.set_mut_one_hot(1<<__builtin_ctz(mut.get_mut_one_hot()));
                }
            }
        }
    });
}
static int set_descendant_count(MAT::Node* root){
    size_t child_count=0;
    for (auto child : root->children) {
        child_count+=set_descendant_count(child);
    }
    root->bfs_index=child_count;
    return child_count;
}

static int leader_thread(std::string& vcf_filename,
    std::string& protobuf_in,
    std::string& protobuf_out,int proc_count){
    auto start_time=std::chrono::steady_clock::now();
    MAT::Tree tree=MAT::load_mutation_annotated_tree(protobuf_in);
    tree.uncondense_leaves();
    std::vector<Sample_Muts> samples_to_place;
    std::vector<mutated_t> position_wise_out;
    Sample_Input(vcf_filename.c_str(),samples_to_place,tree,position_wise_out);
    get_pos_samples_old_tree(tree, position_wise_out);
    fprintf(stderr, "Placing %zu samples \n",samples_to_place.size());
    tree.condense_leaves();
    
    fix_parent(tree.root);
    #ifndef NDEBUG
    Original_State_t ori_state;
    check_samples(tree.root, ori_state, &tree);
    fprintf(stderr, "\n------\n%zu samples\n",ori_state.size());
    #endif
    auto dfs=tree.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,dfs.size()),[&dfs](tbb::blocked_range<size_t> r){
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            auto & mut=dfs[idx]->mutations.mutations;
            mut.erase(std::remove_if(mut.begin(), mut.end(), [](const MAT::Mutation& mut){
                return mut.get_mut_one_hot()==mut.get_par_one_hot();
            }),mut.end());
        }
    });
    for (auto node : dfs) {
        node->branch_length=node->mutations.size();
        #ifdef NDEBUG
        node->children.reserve(SIZE_MULT*node->children.size());
        #endif
    }
    fprintf(stderr, "main dfs size %zu\n",dfs.size());
    std::atomic_size_t curr_idx(0);
    FILE* ignored_file=fopen("/dev/null", "w");
    int optimization_radius=8;
    auto end_time=std::chrono::steady_clock::time_point::max();
    use_bound=true;
    //samples_to_place.resize(1000);
    while (true) {
        assign_descendant_muts(tree);
        assign_levels(tree.root);
        set_descendant_count(tree.root);
        if (proc_count>1) {
            fprintf(stderr, "Main sending tree\n");
            tree.MPI_send_tree();
        }
        place_sample_leader(samples_to_place, tree, 2, proc_count,curr_idx,10000);
        fprintf(stderr, "Main sent optimization prep\n");
        if (process_count==1) {
            Original_State_t origin_states;
            reassign_states(tree, origin_states);
            tree.populate_ignored_range();
        }
        MPI_reassign_states(tree, position_wise_out, 0);
        fprintf(stderr, "Main sent optimization prep done\n");
        std::vector<size_t> node_to_search_idx;
        if (curr_idx>samples_to_place.size()) {
            curr_idx.store(samples_to_place.size());
        }
        find_moved_node_neighbors(optimization_radius, samples_to_place[0].sample_idx, tree, curr_idx.load(), node_to_search_idx);
        fprintf(stderr, "Main found nodes to move\n");
        bool distributed=process_count>1;
        while (!node_to_search_idx.empty()) {
            std::vector<MAT::Node *> deferred_nodes_out;
            adjust_all(tree);
            fprintf(stderr, "Main sent tree_optimizing\n");
            optimize_tree_main_thread(node_to_search_idx, tree, optimization_radius, ignored_file, false, 1, deferred_nodes_out, distributed, end_time, true, true,true);
            node_to_search_idx.clear();
            dfs=tree.depth_first_expansion();
            for (auto node : deferred_nodes_out) {
                node_to_search_idx.push_back(node->dfs_index);
            }
            clean_tree(tree);
            distributed=false;
        }
        dfs=tree.depth_first_expansion();
        for (auto node : dfs) {
            node->mutations.remove_invalid();
        }
        if (curr_idx<samples_to_place.size()) {
            int to_board_cast=0;
            if (process_count>1) {            
                MPI_Bcast(&to_board_cast, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            for (auto node : dfs) {
                if(node->is_leaf()){
                    for (auto& mut : node->mutations) {
                        mut.set_mut_one_hot(mut.get_all_major_allele());
                    }
                }
            }
        }else {
            int to_board_cast=1;
            if (process_count>1) {
                MPI_Bcast(&to_board_cast, 1, MPI_INT, 0, MPI_COMM_WORLD);            
            }
            break;
        }
    }
    fprintf(stderr, "Main finised place\n");
    dfs=tree.depth_first_expansion();
    clean_up_leaf(dfs);
    fix_condensed_nodes(&tree);
    MAT::save_mutation_annotated_tree(tree, protobuf_out);
    auto duration=std::chrono::steady_clock::now()-start_time;
    fprintf(stderr, "Took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
    MPI_Finalize();
    return 0;
}
void wait_debug();
int main(int argc, char **argv) {
    
    //signal(SIGSEGV, handler);
    int ignored;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    std::string vcf_filename;
    std::string protobuf_in;
    std::string protobuf_out;
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    int first_n_sample;
    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_cores) +
                                      " detected on this machine]";
    desc.add_options()(
        "vcf,v", po::value<std::string>(&vcf_filename)->required(),
        "Input VCF file (in uncompressed or gzip-compressed .gz format) "
        "[REQUIRED]")
        ("load-mutation-annotated-tree,i",
                      po::value<std::string>(&protobuf_in)->default_value(""),
                      "Load mutation-annotated tree object")
        ("first_n_sample,n",po::value(&first_n_sample)->default_value(INT_MAX))
        ("save-mutation-annotated-tree,o",
        po::value<std::string>(&protobuf_out)->default_value(""),
        "Save output mutation-annotated tree object to the specified filename")(
        "threads,T",
        po::value<uint32_t>(&num_threads)->default_value(num_cores),
        num_threads_message.c_str())("version", "Print version number")(
        "help,h", "Print help messages");
    po::variables_map vm;
    //wait_debug();
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        // Return with error code 1 unless the user specifies help
        if(vm.count("help")) {
            return 0;
        } else
            return 1;
    }
    

    if (this_rank==0) {
    tbb::task_scheduler_init init(num_threads);
        return leader_thread(vcf_filename,protobuf_in,protobuf_out, process_count);
    }else {
    tbb::task_scheduler_init init(num_threads-1);
        MAT::Tree tree;
        std::vector<mutated_t> position_wise_mutations;
        int start_idx=follower_recieve_positions(position_wise_mutations);
        while (true) {
            fprintf(stderr, "follwer recieving tree\n");
            tree.MPI_receive_tree();
            assign_levels(tree.root);
            set_descendant_count(tree.root);
            auto dfs=tree.depth_first_expansion();
            for (auto node : dfs) {
                check_order(node->mutations);
                #ifdef NDEBUG
                node->children.reserve(SIZE_MULT*node->children.size());
                #endif
            }
            follower_place_sample(tree,2);
            fprintf(stderr, "follower recieving trees\n");
            MPI_reassign_states(tree, position_wise_mutations, start_idx);
            adjust_all(tree);
            use_bound=true;
            fprintf(stderr, "follower start optimizing\n");
            optimize_tree_worker_thread(tree, 8, false, true);
            int stop;
            MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (stop) {
                break;
            }
            tree.delete_nodes();
        }
	    MPI_Finalize();
    }
}
