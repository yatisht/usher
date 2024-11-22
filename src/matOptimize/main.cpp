#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "version.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
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
#include <ios>
#include <iterator>
#include <limits>
#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <thread>
#include <unistd.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <boost/program_options.hpp>
#include <vector>
#include <iostream>
#include <sys/resource.h>
int count_back_mutation(const MAT::Tree& tree);
void get_pos_samples_old_tree(MAT::Tree& tree,std::vector<mutated_t>& output);
int follower_recieve_positions( std::vector<mutated_t>& to_recieve);
void MPI_min_back_reassign_states(MAT::Tree &tree,const std::vector<mutated_t> &mutations,
                                  int start_position);
void min_back_reassign_state_local(MAT::Tree& tree,const std::vector<mutated_t>& mutations);
thread_local TlRng rng;
std::atomic_bool interrupted(false);
bool use_bound;
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
static void load_reference(std::string fasta_fname){
    auto fh=fopen(fasta_fname.c_str(), "r");
    char* seq_name=nullptr;
    size_t seq_len=0;
    auto nchar=getline(&seq_name, &seq_len, fh);
    MAT::Mutation::chromosomes.emplace_back(seq_name+1,seq_name+nchar-1);
    MAT::Mutation::chromosome_map.emplace(MAT::Mutation::chromosomes[0],0);
    free(seq_name);
    auto read=fgetc(fh);
    MAT::Mutation::refs.push_back(0);
    while (read!=EOF) {
        if (read!='\n') {
            auto parsed_nuc=MAT::get_nuc_id(read);
            if (parsed_nuc==0xf) {
                parsed_nuc=0;
            }
            MAT::Mutation::refs.push_back(parsed_nuc);
        }
        read=fgetc(fh);
    }
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
static void remove_sibling(size_t node_id,MAT::Tree& t,std::unordered_map<size_t, bool>& filtered){
    auto this_node=t.get_node(node_id);
    if(this_node->parent){
        auto par_node=this_node->parent;
        auto iter=filtered.find(par_node->node_id);
        if(iter!=filtered.end()){
            iter->second=true;
        }
        for (auto child : par_node->children) {
            if(child!=this_node){
                auto iter=filtered.find(child->node_id);
                if(iter!=filtered.end()){
                    iter->second=true;
                }
            }
        }
    }
}
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
    std::string intermediate_nwk_out="";
    bool reduce_back_mutations=true;
    std::string profitable_src_log;
    std::string ref_file;
    std::string transposed_vcf_path;
    std::string black_list_node_file;
    float search_proportion=2;
    int rand_sel_seed=0;
    unsigned int max_optimize_hours;
    int radius;
    int drift_iterations=0;
    unsigned int minutes_between_save;
    int max_round;
    float min_improvement;
    bool no_write_intermediate;
    std::string diff_file_path;
    std::string branch_support_newick_out;

    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
    ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) ")
    ("tree,t", po::value<std::string>(&input_nh_path)->default_value(""), "Input tree file")
    ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
    ("load-mutation-annotated-tree,i", po::value<std::string>(&input_pb_path)->default_value(""), "Load mutation-annotated tree object")
    ("save-mutation-annotated-tree,o", po::value<std::string>(&output_path), "Save output mutation-annotated tree object to the specified filename [REQUIRED]")
    ("epps_on_branch_len,E", po::value<std::string>(&branch_support_newick_out), "Output a newick with number of equally parsimonious placements on the branch length field ")
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
    ("diff_file_path,D",po::value(&diff_file_path)->default_value(""),"Diff file from MAPLE, used with newick tree (-t)")
    ("reference,R",po::value(&ref_file)->default_value(""),"Reference file, use with diff file (-D)")
    ("version", "Print version number")
    ("node_proportion,z",po::value(&search_proportion)->default_value(2),"the proportion of nodes to search")
    ("node_sel,y",po::value(&rand_sel_seed),"Random seed for selecting nodes to search")
    ("drift_nwk_file,b",po::value(&intermediate_nwk_out)->default_value(""),"Newick filename stem for drifting")
    ("black_list_node_file",po::value(&black_list_node_file)->default_value(""),"Nodes that won't be moved")
    ("no_reduce_back_mutations,c","skip FS that reduce back mutations in the end")
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
    if (vm.count("no_reduce_back_mutations")) {
        reduce_back_mutations=false;
    } else {
        reduce_back_mutations=true;
    }
    if (drift_iterations) {
        min_improvement=0.000000001;
    }
    if(output_path==""&&branch_support_newick_out==""){
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
        fputs("Please either supply output path for protobuf or output EPPs annnotated newick file path",stderr);
        exit(EXIT_FAILURE);
    }
    if (max_optimize_hours) {
        search_end_time=std::chrono::steady_clock::now()+std::chrono::hours(max_optimize_hours)-std::chrono::minutes(30);
        fprintf(stderr, "Set max opt time, will stop in %zu minutes\n",std::chrono::duration_cast<std::chrono::minutes>(search_end_time-std::chrono::steady_clock::now()).count());
    }
    Mutation_Annotated_Tree::Tree t;
    if (this_rank==0) {
        //std::string cwd=get_current_dir_name();
        if(output_path!=""){
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

        if (intermediate_pb_base_name==""&&(!no_write_intermediate)&&output_path!="") {
            intermediate_pb_base_name=output_path+"intermediateXXXXXX.pb";
            auto fd=mkstemps(const_cast<char*>(intermediate_pb_base_name.c_str()), 3);
            close(fd);
        }
        }
        std::string intermediate_writing;
        std::string intermediate_template;
        if(output_path!=""){
            intermediate_template=intermediate_pb_base_name+"temp_eriting_XXXXXX.pb";
        }
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
        } else if (input_nh_path!=""&&diff_file_path!="") {
            print_file_info("Load starting tree from", "starting tree file,-t", input_nh_path);
            print_file_info("Load sample variant from", "sample maple/diff,-D", diff_file_path);
        }
        else {
            fputs("Input file not completely specified. Please either \n"
                  "1. Specify an intermediate protobuf from last run with -a to continue optimization, or\n"
                  "2. Specify a usher-compatible protobuf with -i, and both starting tree and sample variants will be extracted from it, or\n"
                  "3. Specify a usher-compatible protobuf with -i for starting tree, and a VCF with -v for sample variants, or\n"
                  "4. Specify the starting tree in newick format with -t, and either a VCF with -v or a Maple/diff with -D for sample variants.\n",stderr);
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
                        //fprintf(stderr, "Input tree have %zu nodes\n",t.all_nodes.size());
                    } else {
                        if(!MAT::load_mutation_annotated_tree(input_pb_path,t)) {
                            exit(EXIT_FAILURE);
                        }
                        t.uncondense_leaves();
                    }
                    fputs("Finished loading input tree, start reading VCF and assigning states \n",stderr);
                    load_vcf_nh_directly(t, input_vcf_path, origin_states);
                } else if(transposed_vcf_path!="") {
                    if(!MAT::load_mutation_annotated_tree(input_pb_path,t)) {
                        exit(EXIT_FAILURE);
                    }
#ifdef PROFILE_HEAP
                    raise(SIGUSR1);
#endif
                    //malloc_stats();
                    add_ambiguous_mutation(transposed_vcf_path.c_str(),t,false);
#ifdef PROFILE_HEAP
                    raise(SIGUSR1);
#endif
                    //malloc_stats();
                } else if (diff_file_path!="") {
                    if (input_nh_path=="") {
                        fprintf(stderr, "expect newick file\n");
                        exit(EXIT_FAILURE);
                    }
                    if(ref_file==""){
                        fprintf(stderr, "expect reference fasta file\n");
                        exit(EXIT_FAILURE);

                    }
                    load_reference(ref_file);
                    t = Mutation_Annotated_Tree::create_tree_from_newick(
                                input_nh_path);
                    add_ambiguous_mutation(diff_file_path.c_str(),t,true);
                    
                } else {
                    t = load_tree(input_pb_path, origin_states);
                }
            }

#ifndef NDEBUG
            //check_samples(t.root, origin_states, &t);
#endif
            t.populate_ignored_range();
            if(!no_write_intermediate&&input_complete_pb_path==""&&output_path!="") {
                fputs("Checkpoint initial tree.\n",stderr);
                intermediate_writing=intermediate_template;
                make_output_path(intermediate_writing);
                t.save_detailed_mutations(intermediate_writing);
                rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
                fputs("Finished checkpointing initial tree.\n",stderr);
            }
        }
        if (radius==0) {
            save_final_tree(t, output_path);
            return 0;
        }
        std::unordered_set<size_t> node_id_to_ignore;
        if (black_list_node_file!="") {
            std::fstream f(black_list_node_file,std::ios::in);
            std::string temp;
            if(!f){
                auto err_string="unable to open ignored nodes file"+black_list_node_file;
                perror(err_string.c_str());
            }
            while (f) {
                std::getline(f,temp);
                if (temp=="") {
                    continue;
                }
                auto node_id=t.get_node(temp);
                if (node_id) {
                    node_id_to_ignore.insert(node_id->node_id);
                }else {
                    fprintf(stderr, "ignored node %s not in the tree\n", temp.c_str());
                }
            }
        }
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
        if(branch_support_newick_out!=""){
            if(radius<0){
                radius=2*t.get_max_level();
            }
            t.breadth_first_expansion();
            auto all_nodes=t.depth_first_expansion();
            adjust_all(t);
            use_bound=true;
            std::atomic_size_t searched(0);
            auto epp_fh=fopen("epps_dump", "w");;
            tbb::parallel_for(tbb::blocked_range<size_t>(0,all_nodes.size()),[epp_fh,&t,&searched,radius,&all_nodes](tbb::blocked_range<size_t> r){
                for(size_t idx=r.begin();idx<r.end();idx++){
                    output_t out;
                    out.moves=new std::vector<Profitable_Moves_ptr_t>;
                    auto node=all_nodes[idx];
                    Reachable reachable{true,true};
                    find_moves_bounded(node, out, radius, true, reachable);
                    std::unordered_map<size_t, bool> filtered;
                    if(out.score_change==-1){
                        filtered.emplace(node->node_id,false);
                    }
                    for (const auto& move : *out.moves) {
                        filtered.emplace(move->dst->node_id,false);
                    }
                    if(out.score_change==-1){
                        remove_sibling(node->node_id, t, filtered);
                    }
                    for(auto& node_id_filtered:filtered){
                        if(node_id_filtered.second){
                            continue;
                        }
                        remove_sibling(node_id_filtered.first, t, filtered);
                    }
                    std::vector<size_t> filtered_nodes;
                    filtered_nodes.reserve(filtered.size());
                    for (auto& node_id_filtered : filtered) {
                        if(!node_id_filtered.second){
                            filtered_nodes.emplace_back(node_id_filtered.first);
                        }
                    }
                    node->branch_length=filtered_nodes.size();
                    if (node->branch_length<1||filtered_nodes.size()>(out.moves->size()+1)) {
                        raise(SIGTRAP);
                    }
                    if(filtered_nodes.size()>1){
                        std::string out_str=t.get_node_name_for_log_output(node->node_id)+":";
                        for (const auto& node_id : filtered_nodes) {
                            if(node_id==node->node_id){
                                continue;
                            }
                            out_str+=(t.get_node_name_for_log_output(node_id)+",");
                        }
                        out_str.pop_back();
                        out_str.push_back('\n');
                        fputs(out_str.c_str(), epp_fh);
                    }
                    delete out.moves;
                }
                searched+=r.size();
                printf("searched %zu out of %zu\r",searched.load(std::memory_order_relaxed),all_nodes.size());

            });
            fclose(epp_fh);
            std::fstream out_f(branch_support_newick_out,std::ios::out);
            out_f<<t.get_newick_string(true,true,true,true);
            return EXIT_SUCCESS;
        }
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
            if(!node_id_to_ignore.empty()){
                nodes_to_search.erase(std::remove_if(nodes_to_search.begin(), nodes_to_search.end(), 
                    [&node_id_to_ignore](MAT::Node* node){
                        return node_id_to_ignore.find(node->node_id)!=node_id_to_ignore.end();
                    }),nodes_to_search.end());
            }
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
            //Actual optimization loop
            new_score=optimize_inner_loop(nodes_to_search,t,radius,
#ifdef CHECK_STATE_REASSIGN
            origin_states,
#endif
            allow_drift,search_all_dir,minutes_between_save,no_write_intermediate,search_end_time,start_time,log_moves,iteration,intermediate_template,intermediate_pb_base_name,intermediate_nwk_out);
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
        if (reduce_back_mutations) {
            fprintf(stderr, "Parsimony score before %zu\n",t.get_parsimony_score());
            fprintf(stderr, "Back mutation count before %d\n",count_back_mutation(t));
            std::vector<mutated_t> output(MAT::Mutation::refs.size());
            get_pos_samples_old_tree(t, output);
            if (process_count==1) {
                min_back_reassign_state_local(t,output);
            } else {
                MPI_min_back_reassign_states(t, output, 0);
            }
            fprintf(stderr, "Back mutation count after %d\n",count_back_mutation(t));
        }
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
        if (reduce_back_mutations) {
            std::vector<mutated_t> to_recieve;
            int pos= follower_recieve_positions(to_recieve);
            MPI_min_back_reassign_states(t, to_recieve, pos);
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
