#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "version.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <chrono>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <ctime>
//#include <malloc.h>
#include <fstream>
#include <iterator>
#include <limits>
#include <mpi.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <thread>
#include <unistd.h>
#include <sys/stat.h>
#include <unordered_set>
#include <boost/program_options.hpp>
#include <vector>
#include <iostream>
#include <sys/resource.h>
thread_local TlRng rng;
std::atomic_bool interrupted(false);
bool use_bound;
size_t get_memory() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}
int process_count;
int this_rank;
uint32_t num_threads;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
FILE* movalbe_src_log;
bool changing_radius=false;
void interrupt_handler(int) {
    fputs("interrupted\n", stderr);
    interrupted=true;
    fflush(movalbe_src_log);
}
void log_flush(int) {
    fflush(movalbe_src_log);
}
namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
static void make_output_path(std::string& path_template) {
    auto fd=mkstemps(const_cast<char*>(path_template.c_str()),3);
    close(fd);
}

void print_file_info(std::string info_msg,std::string error_msg,const std::string& filename) {
    struct stat stat_buf;
    errno=0;
    if(stat(filename.c_str(), &stat_buf)==0) {
        fprintf(stderr,"%s %s last modified: %s",info_msg.c_str(),filename.c_str(),ctime(&stat_buf.st_mtime));
    } else {
        perror(("Error accessing "+error_msg+" "+filename+" :").c_str());
        exit(EXIT_FAILURE);
    }
}
#define DRIFT_MASK 0x80000000
#define ALL_DIR_MASK 0x40000000
#define RADIUS_MASK 0x3fffffff
int main(int argc, char **argv) {
    int ignored;
    auto init_result=MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&ignored);
    if (init_result!=MPI_SUCCESS) {
        fprintf(stderr, "MPI init failed\n");
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    fprintf(stderr, "Running with %d processes\n",process_count);
    std::string output_path;
    std::string input_pb_path;
    std::string input_complete_pb_path;
    std::string input_nh_path;
    std::string input_vcf_path;
    std::string intermediate_pb_base_name="";
    std::string profitable_src_log;
    std::string transposed_vcf_path;
    float search_proportion=2;
    int rand_sel_seed=0;
    unsigned int max_optimize_hours;
    int radius;
    int drift_iterations=0;
    unsigned int minutes_between_save;
    int max_round;
    float min_improvement;
    bool no_write_intermediate;

    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
    ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) ")
    ("tree,t", po::value<std::string>(&input_nh_path)->default_value(""), "Input tree file")
    ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
    ("load-mutation-annotated-tree,i", po::value<std::string>(&input_pb_path)->default_value(""), "Load mutation-annotated tree object")
    ("save-mutation-annotated-tree,o", po::value<std::string>(&output_path)->required(), "Save output mutation-annotated tree object to the specified filename [REQUIRED]")
    ("radius,r", po::value<int32_t>(&radius)->default_value(-1),
     "Radius in which to restrict the SPR moves.")
    ("profitable-src-log,S", po::value<std::string>(&profitable_src_log)->default_value("/dev/null"),
     "The file to log from which node a profitable move can be found.")
    ("ambi-protobuf,a", po::value<std::string>(&input_complete_pb_path)->default_value(""),
     "Continue from intermediate protobuf")
    ("minutes-between-save,s",po::value<unsigned int>(&minutes_between_save)->default_value(0),"Minutes between saving intermediate protobuf")
    ("min-improvement,m",po::value<float>(&min_improvement)->default_value(0.0005),"Minimum improvement in the parsimony score as a fraction of the previous score in ordder to perform another iteration.")
    ("drift_iteration,d",po::value<int>(&drift_iterations)->default_value(0),"Iterations permiting equally parsimonious moves after parsimony score no longer improves")
    ("do-not-write-intermediate-files,n","Do not write intermediate files.")
    ("max-iterations,N", po::value<int>(&max_round)->default_value(1000), \
     "Maximum number of optimization iterations to perform.")
    ("max-hours,M",po::value(&max_optimize_hours)->default_value(0),"Maximium number of hours to run")
    ("transposed-vcf-path,V",po::value(&transposed_vcf_path)->default_value(""),"Auxiliary transposed VCF for ambiguous bases, used in combination with usher protobuf (-i)")
    ("version", "Print version number")
    ("node_proportion,z",po::value(&search_proportion)->default_value(2),"the proportion of nodes to search")
    ("node_sel,y",po::value(&rand_sel_seed),"Random seed for selecting nodes to search")
    ("help,h", "Print help messages");
    auto search_end_time=std::chrono::steady_clock::time_point::max();
    po::options_description all_options;
    all_options.add(desc);
    signal(SIGUSR2,interrupt_handler);
    signal(SIGUSR1,log_flush);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        if (this_rank==0) {
            if (vm.count("version")) {
                std::cout << "matOptimize (v" << PROJECT_VERSION << ")" << std::endl;
            } else {
                std::cerr << "matOptimize (v" << PROJECT_VERSION << ")" << std::endl;
                std::cerr << desc << std::endl;
            }
        }
        MPI_Finalize();
        // Return with error code 1 unless the user specifies help
        if(vm.count("help")) {
            return 0;
        } else
            return 1;
    }
    if (max_optimize_hours) {
        search_end_time=std::chrono::steady_clock::now()+std::chrono::hours(max_optimize_hours)-std::chrono::minutes(30);
        fprintf(stderr, "Set max opt time, will stop in %zu minutes\n",std::chrono::duration_cast<std::chrono::minutes>(search_end_time-std::chrono::steady_clock::now()).count());
    }
    Mutation_Annotated_Tree::Tree t;
    if (this_rank==0) {
        //std::string cwd=get_current_dir_name();
        no_write_intermediate=vm.count("do-not-write-intermediate-files");
        try {
            auto output_path_dir_name=boost::filesystem::system_complete(output_path).parent_path();
            if(!boost::filesystem::exists(output_path_dir_name)) {
                boost::filesystem::create_directories(output_path_dir_name);
            }
            errno=0;
        } catch(const boost::filesystem::filesystem_error& ex) {
            std::cerr<<"Cannot create parent directory of output path\n";
            std::cerr<<ex.what()<<std::endl;
            exit(EXIT_FAILURE);
        }
        auto fd=open(output_path.c_str(),O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR);
        if (fd==-1) {
            perror(("Cannot create output file"+output_path).c_str());
            exit(EXIT_FAILURE);
        } else {
            close(fd);
        }

        if (intermediate_pb_base_name==""&&(!no_write_intermediate)) {
            intermediate_pb_base_name=output_path+"intermediateXXXXXX.pb";
            auto fd=mkstemps(const_cast<char*>(intermediate_pb_base_name.c_str()), 3);
            close(fd);
        }
        std::string intermediate_writing;
        std::string intermediate_template;
        intermediate_template=intermediate_pb_base_name+"temp_eriting_XXXXXX.pb";
        fputs("Summary:\n",stderr);
        if (input_complete_pb_path!="") {
            print_file_info("Continue from", "continuation protobut,-a", input_complete_pb_path);
        } else if (input_pb_path!="") {
            if (input_vcf_path=="") {
                print_file_info("Load both starting tree and sample variant from", "input protobut,-i", input_pb_path);
            } else {
                print_file_info("Extract starting tree from", "input protobut,-i", input_pb_path);
                print_file_info("Load sample variant from", "sample vcf,-v", input_vcf_path);
            }
        } else if (input_nh_path!=""&&input_vcf_path!="") {
            print_file_info("Load starting tree from", "starting tree file,-t", input_nh_path);
            print_file_info("Load sample variant from", "sample vcf,-v", input_vcf_path);
        } else {
            fputs("Input file not completely specified. Please either \n"
                  "1. Specify an intermediate protobuf from last run with -a to continue optimization, or\n"
                  "2. Specify a usher-compatible protobuf with -i, and both starting tree and sample variants will be extracted from it, or\n"
                  "3. Specify a usher-compatible protobuf with -i for starting tree, and a VCF with -v for sample variants, or\n"
                  "4. Specify the starting tree in newick format with -t, and a VCF with -v for sample variants.\n",stderr);
            std::cerr << desc << std::endl;
            exit(EXIT_FAILURE);
        }
        fprintf(stderr,"Will output final protobuf to %s .\n",output_path.c_str());
        if (no_write_intermediate) {
            fprintf(stderr,"Will not write intermediate file. WARNNING: It will not be possible to continue optimization if this program is killed in the middle.\n");
        } else {
            fprintf(stderr,"Will output intermediate protobuf to %s. \n",intermediate_pb_base_name.c_str());
        }
        auto pid=getpid();
        bool log_moves=false;
        if (profitable_src_log!="/dev/null") {
            log_moves=true;
            fprintf(stderr,"Will write list of source nodes where a profitable move can be found to %s. \n",profitable_src_log.c_str());
            fprintf(stderr,"Run kill -s SIGUSR1 %d to flush the source node log\n",pid);
        }
        if (radius>=0) {
            fprintf(stderr,"Will consider SPR moves within a radius of %d. \n",radius);
        } else {
            fprintf(stderr,"Will double radius after each iteration\n");
            changing_radius=true;
        }
        fprintf(stderr,"Run kill -s SIGUSR2 %d to apply all the move found immediately, then output and exit.\n",pid);
        fprintf(stderr,"Using %d threads. \n",num_threads);
        auto start_time=std::chrono::steady_clock::now();
        auto save_period=std::chrono::minutes(minutes_between_save);
        fprintf(stderr, "Will drift for %d iterations \n",drift_iterations);

        //Loading tree
        Original_State_t origin_states;
        {
            tbb::task_scheduler_init init(process_count*num_threads);
            if (input_complete_pb_path!="") {
                t.load_detatiled_mutations(input_complete_pb_path);
            } else {
                if (input_vcf_path != "") {
                    fputs("Loading input tree\n",stderr);
                    if (input_nh_path != "") {
                        t = Mutation_Annotated_Tree::create_tree_from_newick(
                                input_nh_path);
                        fprintf(stderr, "Input tree have %zu nodes\n",t.all_nodes.size());
                    } else {
                        t = MAT::load_mutation_annotated_tree(input_pb_path);
                        t.uncondense_leaves();
                    }
                    fputs("Finished loading input tree, start reading VCF and assigning states \n",stderr);
                    load_vcf_nh_directly(t, input_vcf_path, origin_states);
                } else if(transposed_vcf_path!="") {
                    t=MAT::load_mutation_annotated_tree(input_pb_path);
#ifdef PROFILE_HEAP
                    raise(SIGUSR1);
#endif
                    //malloc_stats();
                    add_ambiguous_mutation(transposed_vcf_path.c_str(),t);
#ifdef PROFILE_HEAP
                    raise(SIGUSR1);
#endif
                    //malloc_stats();
                } else {
                    t = load_tree(input_pb_path, origin_states);
                }
            }

#ifndef NDEBUG
            //check_samples(t.root, origin_states, &t);
#endif
            t.populate_ignored_range();
            if(!no_write_intermediate&&input_complete_pb_path=="") {
                fputs("Checkpoint initial tree.\n",stderr);
                intermediate_writing=intermediate_template;
                make_output_path(intermediate_writing);
                t.save_detailed_mutations(intermediate_writing);
                rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
                fputs("Finished checkpointing initial tree.\n",stderr);
            }
        }
        auto last_save_time=std::chrono::steady_clock::now();
        size_t new_score;
        size_t score_before;
        int stalled = -1;
        score_before = t.get_parsimony_score();
        new_score = score_before;
        fprintf(stderr, "after state reassignment:%zu\n", score_before);
        fprintf(stderr, "Height:%zu\n", t.get_max_level());

        std::vector<MAT::Node *> nodes_to_search;
        std::vector<MAT::Node *> bfs_ordered_nodes;
        bfs_ordered_nodes = t.breadth_first_expansion();
        movalbe_src_log=fopen(profitable_src_log.c_str(),"w");
        if (!movalbe_src_log) {
            perror(("Error writing to log file "+profitable_src_log).c_str());
            movalbe_src_log=fopen("/dev/null", "w");
        }
        fprintf(movalbe_src_log, "source\tdestination\titeration\tscore.change\tdistance\tsubtree.size\n");
        bool isfirst=true;
        bool allow_drift=false;
        int iteration=1;
        tbb::task_scheduler_init init(num_threads);
        while(stalled<drift_iterations) {
            bfs_ordered_nodes = t.breadth_first_expansion();
            fputs("Start Finding nodes to move \n",stderr);
            bool search_all_nodes=false;
            bool search_all_dir=false;
            if (isfirst||allow_drift) {
                search_all_nodes=true;
                search_all_dir=true;
            } else if (radius<0&&radius>=-2*(int)t.max_level) {
                radius*=2;
                search_all_nodes=true;
                search_all_dir=true;
            }
            find_nodes_to_move(bfs_ordered_nodes, nodes_to_search,search_all_nodes,search_all_dir,radius,t);
            if (search_proportion<1) {
                std::vector<MAT::Node *> nodes_to_search_temp;
                nodes_to_search_temp.reserve(nodes_to_search.size()*search_proportion);
                std::mt19937_64 rng(rand_sel_seed);
                std::sample(nodes_to_search.begin(), nodes_to_search.end(), std::back_inserter(nodes_to_search_temp),size_t(std::round(nodes_to_search.size()*search_proportion)),rng);
                nodes_to_search.swap(nodes_to_search_temp);
            }
            isfirst=false;
            fprintf(stderr,"%zu nodes to search\n",nodes_to_search.size());
            if (nodes_to_search.empty()) {
                break;
            }
            bool isfirst_this_iter=true;
            //Actual optimization loop
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
                std::vector<MAT::Node*> defered_nodes;
                auto next_save_time=minutes_between_save?last_save_time+save_period:std::chrono::steady_clock::time_point::max();
                bool do_continue=true;
                auto search_stop_time=next_save_time;
                if (no_write_intermediate||search_end_time<next_save_time) {
                    search_stop_time=search_end_time;
                }
                optimize_tree_main_thread(nodes_to_search_idx, t,std::abs(radius),movalbe_src_log,allow_drift,log_moves?iteration:-1,defered_nodes,distribute,search_stop_time,do_continue,search_all_dir,isfirst_this_iter
#ifndef NDEBUG
                                          , origin_states
#endif
                                         );
                isfirst_this_iter=false;
                fprintf(stderr, "Defered %zu nodes\n",defered_nodes.size());
                nodes_to_search=std::move(defered_nodes);
                new_score=t.get_parsimony_score();
                fprintf(stderr, "parsimony score after optimizing: %zu,with radius %d, second from start %ld \n\n",
                        new_score,std::abs(radius),std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_time).count());
                if(!no_write_intermediate) {
                    intermediate_writing=intermediate_template;
                    make_output_path(intermediate_writing);
                    auto save_start=std::chrono::steady_clock::now();
                    t.save_detailed_mutations(intermediate_writing);
                    rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
                    last_save_time=std::chrono::steady_clock::now();
                    fprintf(stderr, "Took %ldsecond to save intermediate protobuf\n",std::chrono::duration_cast<std::chrono::seconds>(last_save_time-save_start).count());
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
            if (interrupted) {
                break;
            }
            float improvement=1-((float)new_score/(float)score_before);
            fprintf(stderr, "Last round improvement %f\n",improvement);
            if (improvement <min_improvement) {
                fprintf(stderr, "Less than minimium improvement,stalled for %d iterations\n",stalled);
                fprintf(stderr, "Will drift for %d iterations \n",drift_iterations);
                stalled++;
                allow_drift=true;
            } else {
                score_before = new_score;
                stalled = -1;
            }
            iteration++;
            if(iteration>=max_round) {
                fprintf(stderr, "Reached %d interations\n", iteration);
                break;
            }
            if (std::chrono::steady_clock::now()>=search_end_time) {
                fprintf(stderr, "Exceeded search time \n");
                break;
            }
        }
        int temp=0;
        MPI_Request req;
        MPI_Ibcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD,&req);
        fprintf(stderr, "Final Parsimony score %zu\n",t.get_parsimony_score());
        fclose(movalbe_src_log);
        save_final_tree(t, output_path);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    } else {
        tbb::task_scheduler_init init(num_threads);
        while(true) {
            MPI_Request request;
            fprintf(stderr, "Wait radius\n");
            MPI_Ibcast(&radius, 1, MPI_INT, 0, MPI_COMM_WORLD,&request);
            int done=0;
            while (done==0) {
                MPI_Test(&request, &done, MPI_STATUS_IGNORE);
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            bool do_drift=false;
            bool search_all_dir=false;
            if (radius&DRIFT_MASK) {
                do_drift=true;
            }
            if (radius&ALL_DIR_MASK) {
                search_all_dir=true;
                fprintf(stderr, "Search all directions\n");
            }
            radius&=RADIUS_MASK;
            fprintf(stderr, "recieved radius %d \n",radius);
            if(radius==0) {
                break;
            }
            t.MPI_receive_tree();
            adjust_all(t);
            use_bound=true;
            optimize_tree_worker_thread(t, radius,do_drift,search_all_dir);
            t.delete_nodes();
        }
    }
    for(auto& pos:mutated_positions) {
        delete pos.second;
    }
    if (this_rank==0) {
        t.delete_nodes();
    }
    fprintf(stderr, "Maximum memory usage from %d: %zu kb \n",this_rank,get_memory());
    MPI_Finalize();
}
