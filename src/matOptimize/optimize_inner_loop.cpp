#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "tree_rearrangement_internal.hpp"
#include <mpi.h>
#include <cstdlib>
#include <unistd.h>
void make_output_path(std::string& path_template) {
    auto fd=mkstemps(const_cast<char*>(path_template.c_str()),3);
    close(fd);
}

size_t optimize_inner_loop(std::vector<MAT::Node*>& nodes_to_search,MAT::Tree& t,int radius,
    Move_Found_Callback& callback,
    bool allow_drift,
    bool search_all_dir,
    int minutes_between_save,
    bool no_write_intermediate,
    std::chrono::steady_clock::time_point search_end_time,
    std::chrono::steady_clock::time_point start_time,
    bool log_moves,
    int iteration,
    std::string intermediate_template,
    std::string intermediate_pb_base_name,
    std::string intermediate_nwk_out
    ){
    static std::chrono::steady_clock::time_point last_save_time=std::chrono::steady_clock::now();
    auto save_period=std::chrono::minutes(minutes_between_save);
    bool isfirst_this_iter=true;
    size_t new_score;
     while (!nodes_to_search.empty()) {
                auto dfs_ordered_nodes=t.depth_first_expansion();
                std::mt19937_64 rng;
                std::shuffle(nodes_to_search.begin(), nodes_to_search.end(),rng);
                bool distribute=(process_count>1)&&(nodes_to_search.size()>1000);
                if (distribute) {
                    MPI_Request req;
                    int radius_to_boardcast=abs(radius);
                    if (allow_drift) {
                        radius_to_boardcast|=DRIFT_MASK;
                    }
                    if(search_all_dir) {
                        radius_to_boardcast|=ALL_DIR_MASK;
                        fprintf(stderr, "Search all directions\n");
                    }
                    MPI_Ibcast(&radius_to_boardcast, 1, MPI_INT, 0, MPI_COMM_WORLD, &req);
                    fprintf(stderr, "Sent radius\n");
                    MPI_Wait(&req, MPI_STATUS_IGNORE);
                    fprintf(stderr, "Start Send tree\n");
                    t.MPI_send_tree();
                }
                adjust_all(t);
                use_bound=true;
                std::vector<size_t> nodes_to_search_idx;
                nodes_to_search_idx.reserve(nodes_to_search.size());
                for(const auto node:nodes_to_search) {
                    nodes_to_search_idx.push_back(node->dfs_index);
                }
                std::vector<size_t> defered_nodes;
                auto next_save_time=minutes_between_save?last_save_time+save_period:std::chrono::steady_clock::time_point::max();
                bool do_continue=true;
                auto search_stop_time=next_save_time;
                if (no_write_intermediate||search_end_time<next_save_time) {
                    search_stop_time=search_end_time;
                }
                optimize_tree_main_thread(nodes_to_search_idx, t,std::abs(radius),movalbe_src_log,allow_drift,log_moves?iteration:-1,defered_nodes,distribute,search_stop_time,do_continue,search_all_dir,isfirst_this_iter, callback
                                         );
                isfirst_this_iter=false;
                fprintf(stderr, "Defered %zu nodes\n",defered_nodes.size());
                nodes_to_search.reserve(defered_nodes.size());
                nodes_to_search.clear();
                for (auto idx : defered_nodes) {
                    if (t.get_node(idx)) {
                        nodes_to_search.push_back(t.get_node(idx));
                    }
                }
                auto curr_score=t.get_parsimony_score();
                if(curr_score>=new_score) {
                    nodes_to_search.clear();
                }
                new_score=curr_score;
                fprintf(stderr, "parsimony score after optimizing: %zu,with radius %d, second from start %ld \n\n",
                        new_score,std::abs(radius),std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_time).count());
                if(!no_write_intermediate) {
                    auto intermediate_writing=intermediate_template;
                    make_output_path(intermediate_writing);
                    auto save_start=std::chrono::steady_clock::now();
                    t.save_detailed_mutations(intermediate_writing);
                    rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
                    last_save_time=std::chrono::steady_clock::now();
                    fprintf(stderr, "Took %ldsecond to save intermediate protobuf\n",std::chrono::duration_cast<std::chrono::seconds>(last_save_time-save_start).count());
                }
                if(allow_drift) {
                    MAT::save_mutation_annotated_tree(t, intermediate_nwk_out+std::to_string(iteration)+".pb.gz");
                }
                if (std::chrono::steady_clock::now()>=search_end_time) {
                    break;
                }
                if (interrupted) {
                    break;
                }
                if (allow_drift) {
                    nodes_to_search.clear();
                }
                search_all_dir=true;
            }
            return new_score;
}
