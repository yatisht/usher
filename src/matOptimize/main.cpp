#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "version.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
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
#include <limits>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <unistd.h>
#include <sys/stat.h>
#include <unordered_set>
#include <boost/program_options.hpp>
#include <vector>
#include <iostream>
#include <sys/resource.h>
thread_local TlRng rng;
void print_memory() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    fprintf(stderr, "Max resident %zu kb\n",usage.ru_maxrss);
}
uint32_t num_threads;
std::chrono::time_point<std::chrono::steady_clock> last_save_time;
bool no_write_intermediate;
size_t max_queued_moves;
std::chrono::steady_clock::duration save_period;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
FILE* movalbe_src_log;
bool interrupted;
bool changing_radius=false;
std::condition_variable progress_bar_cv;
bool timed_print_progress;
tbb::task_group_context search_context;
void interrupt_handler(int) {
    fputs("interrupted\n", stderr);
    search_context.cancel_group_execution();
    interrupted=true;
    fflush(movalbe_src_log);
}

void log_flush_handle(int) {
    fflush(movalbe_src_log);
    progress_bar_cv.notify_all();
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
int main(int argc, char **argv) {
    std::string output_path;
    std::string input_pb_path;
    std::string input_complete_pb_path;
    std::string input_nh_path;
    std::string input_vcf_path;
    std::string intermediate_pb_base_name="";
    std::string profitable_src_log;
    std::chrono::steady_clock::duration ori_save_period;
    std::string transposed_vcf_path;
    int drift_iter;
    unsigned int max_optimize_hours;
    int radius;
    unsigned int minutes_between_save;
    int max_round;
    float min_improvement;

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
    ("max-queued-moves,q",po::value<size_t>(&max_queued_moves)->default_value(INT_MAX),"Maximium number of profitable moves found before applying these moves")
    ("minutes-between-save,s",po::value<unsigned int>(&minutes_between_save)->default_value(0),"Minutes between saving intermediate protobuf")
    ("min-improvement,m",po::value<float>(&min_improvement)->default_value(0.0005),"Minimum improvement in the parsimony score as a fraction of the previous score in ordder to perform another iteration.")
    ("do-not-write-intermediate-files,n","Do not write intermediate files.")
    ("max-iterations,N", po::value<int>(&max_round)->default_value(1000), \
     "Maximum number of optimization iterations to perform.")
    ("exhaustive-mode,e","Search every non-root node as source node.")
    ("max-hours,M",po::value(&max_optimize_hours)->default_value(0),"Maximium number of hours to run")
    ("transposed-vcf-path,V",po::value(&transposed_vcf_path)->default_value(""),"Auxiliary transposed VCF for ambiguous bases, used in combination with usher protobuf (-i)")
    ("version", "Print version number")
    ("drift_iter,d",po::value(&drift_iter)->default_value(1),"Number of iteration to continue if no parsimony improvement")
    ("help,h", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);
    interrupted=false;
    signal(SIGUSR2,interrupt_handler);
#ifndef PROFILE_HEAP
    signal(SIGUSR1, log_flush_handle);
#endif
    po::variables_map vm;
    if (argc==1) {
        std::cerr << desc << std::endl;
        return EXIT_FAILURE;
    }
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        if (vm.count("version")) {
            std::cout << "matOptimize (v" << PROJECT_VERSION << ")" << std::endl;
        } else {
            std::cerr << "matOptimize (v" << PROJECT_VERSION << ")" << std::endl;
            std::cerr << desc << std::endl;
        }

        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }
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
    timed_print_progress=isatty(fileno(stderr));
    if (!timed_print_progress) {
        fprintf(stderr,"stderr is not a terminal, will not print progress\n");
        fprintf(stderr,"Run kill -s SIGUSR1 %d to print progress\n",pid);
    }
    if (radius>=0) {
        fprintf(stderr,"Will consider SPR moves within a radius of %d. \n",radius);
    } else {
        fprintf(stderr,"Will double radius after each iteration\n");
        changing_radius=true;
    }
    if (max_queued_moves) {
        fprintf(stderr, "Will stop search and apply all moves found when %zu moves are found\n",max_queued_moves);
    } else {
        max_queued_moves=std::numeric_limits<decltype(max_queued_moves)>::max();
    }
    if (minutes_between_save) {
        fprintf(stderr, "Will save intermediate result every %d minutes\n",minutes_between_save);
        ori_save_period=std::chrono::minutes(minutes_between_save);
    } else {
        ori_save_period=ori_save_period.max();
    }
    fprintf(stderr,"Run kill -s SIGUSR2 %d to apply all the move found immediately, then output and exit.\n",pid);
    fprintf(stderr,"Using %d threads. \n",num_threads);
    std::chrono::steady_clock::duration max_optimize_duration=std::chrono::hours(max_optimize_hours)-std::chrono::minutes(30);
    auto start_time=std::chrono::steady_clock::now();
    auto end_time=start_time+max_optimize_duration;
    tbb::task_scheduler_init init(num_threads);

    //Loading tree
    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t;

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
            print_memory();
            //malloc_stats();
            add_ambiguous_mutation(transposed_vcf_path.c_str(),t);
#ifdef PROFILE_HEAP
            raise(SIGUSR1);
#endif
            print_memory();
            //malloc_stats();
        } else {
            t = load_tree(input_pb_path, origin_states);
        }
        if(!no_write_intermediate) {
            fputs("Checkpoint initial tree.\n",stderr);
            intermediate_writing=intermediate_template;
            make_output_path(intermediate_writing);
            t.save_detailed_mutations(intermediate_writing);
            rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
            fputs("Finished checkpointing initial tree.\n",stderr);
        }
    }
    last_save_time=std::chrono::steady_clock::now();
    size_t new_score;
    size_t score_before;
    int stalled = 0;

#ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
#endif

    score_before = t.get_parsimony_score();
    new_score = score_before;
    fprintf(stderr, "after state reassignment:%zu\n", score_before);
    fprintf(stderr, "Height:%zu\n", t.get_max_level());

    tbb::concurrent_vector<MAT::Node *> nodes_to_search;
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
    size_t nodes_seached_last_iter=0;
    size_t nodes_seached_this_iter=0;
    while(stalled<drift_iter) {
        if (interrupted) {
            break;
        }
        bfs_ordered_nodes = t.breadth_first_expansion();
        fputs("Start Finding nodes to move \n",stderr);
        find_nodes_to_move(bfs_ordered_nodes, nodes_to_search,isfirst,radius,t);
        if (radius<0&&(isfirst||nodes_seached_this_iter>nodes_seached_last_iter)) {
            radius*=2;
            nodes_seached_last_iter=nodes_seached_this_iter;
        }
        isfirst=false;
        fprintf(stderr,"%zu nodes to search\n",nodes_to_search.size());
        if (allow_drift) {
            nodes_to_search=tbb::concurrent_vector<MAT::Node *>(bfs_ordered_nodes.begin(),bfs_ordered_nodes.end());
        }
        if (nodes_to_search.empty()) {
            break;
        }
        //Actual optimization loop
        bool searched_full=true;
        while (!nodes_to_search.empty()) {
            if (interrupted) {
                break;
            }
            bfs_ordered_nodes = t.breadth_first_expansion();
            if (max_optimize_hours) {
                save_period=std::min(ori_save_period,end_time-std::chrono::steady_clock::now());
            } else {
                save_period=ori_save_period;
            }
            auto res =
                optimize_tree(bfs_ordered_nodes, nodes_to_search, t,std::abs(radius),movalbe_src_log,allow_drift,log_moves?iteration:-1
#ifndef NDEBUG
                              , origin_states
#endif
                             );
            if (max_optimize_hours&&(std::chrono::steady_clock::now()-start_time>max_optimize_duration)) {
                interrupted=true;
                break;
            }
            new_score=res.first;
            if (searched_full) {
                nodes_seached_this_iter=res.second;
                searched_full=false;
            }
            fprintf(stderr, "parsimony score after optimizing: %zu,with radius %d, searched %zu arcs \n\n", new_score,std::abs(radius),res.second);
            auto save_start=std::chrono::steady_clock::now();
            if(std::chrono::steady_clock::now()-last_save_time>=save_period) {
                if(!no_write_intermediate) {
                    intermediate_writing=intermediate_template;
                    make_output_path(intermediate_writing);
                    t.save_detailed_mutations(intermediate_writing);
                    rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
                    last_save_time=std::chrono::steady_clock::now();
                    fprintf(stderr, "Took %ldsecond to save intermediate protobuf\n",std::chrono::duration_cast<std::chrono::seconds>(last_save_time-save_start).count());
                }
                last_save_time=std::chrono::steady_clock::now();
            }
        }
        float improvement=1-((float)new_score/(float)score_before);
        fprintf(stderr, "Last round improvement %f\n",improvement);
        if (improvement <min_improvement ) {
            if (nodes_seached_this_iter>nodes_seached_last_iter&&radius<0) {
                radius*=2;
            } else {
                fprintf(stderr, "Less than minimium improvement\n");
                stalled++;
                allow_drift=true;
            }
        } else {
            score_before = new_score;
            stalled = 0;
        }
        clean_tree(t);
        iteration++;
        if (max_optimize_hours&&(std::chrono::steady_clock::now()-start_time>max_optimize_duration)) {
            interrupted=true;
            break;
        }
        if(iteration>=max_round) {
            fprintf(stderr, "Reached %d interations\n", iteration);
            interrupted=true;
            break;
        }
    }
    fprintf(stderr, "Final Parsimony score %zu\n",t.get_parsimony_score());
    fclose(movalbe_src_log);
    save_final_tree(t, origin_states, output_path);
    for(auto& pos:mutated_positions) {
        delete pos.second;
    }
    t.delete_nodes();
}
