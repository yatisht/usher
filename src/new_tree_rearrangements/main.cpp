#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <tbb/task.h>
#include <cstdio>
#include <string>

#include <tbb/task_scheduler_init.h>
#include <unordered_set>
#include <boost/program_options.hpp> 
#include <vector>
#include <iostream>
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
FILE* movalbe_src_log;
int early_stop_saving;
bool interrupted;
tbb::task_group_context search_context;
void interrupt_handler(int){
    fputs("interrupted\n", stderr);
    search_context.cancel_group_execution();
    interrupted=true;
    fflush(movalbe_src_log);
}

void log_flush_handle(int){
    fflush(movalbe_src_log);
}
namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
int main(int argc, char **argv) {
    std::string output_path;
    std::string input_pb_path;
    std::string input_complete_pb_path;
    std::string input_nh_path;
    std::string input_vcf_path;
    std::string intermediate_pb_base_name;
    std::string profitable_src_log;
    int radius;
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    early_stop_saving=0;
    desc.add_options()
        ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
        ("tree,t", po::value<std::string>(&input_nh_path)->default_value(""), "Input tree file")
        ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("load-mutation-annotated-tree,i", po::value<std::string>(&input_pb_path)->default_value(""), "Load mutation-annotated tree object")
        ("save-mutation-annotated-tree,o", po::value<std::string>(&output_path)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
        ("save-intermediate-mutation-annotated-tree,m", po::value<std::string>(&intermediate_pb_base_name)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
        ("radius,r", po::value<int32_t>(&radius)->default_value(10),
         "Radius in which to restrict the SPR moves.")
        ("profitable_src_log,s", po::value<std::string>(&profitable_src_log)->default_value("/dev/null"),
         "The file to log from which node a profitable move can be found.")
        ("ambi_protobuf,a", po::value<std::string>(&input_complete_pb_path)->default_value(""),
         "Continue from intermediate protobuf");
        
    po::options_description all_options;
    all_options.add(desc);
    interrupted=false;
    signal(SIGINT,interrupt_handler);
    signal(SIGUSR1, log_flush_handle);
    po::variables_map vm;
    if (argc==1) {
        std::cerr << desc << std::endl;
        return EXIT_FAILURE;
    }
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    tbb::task_scheduler_init init(num_threads);

    //Loading tree
    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t;
    
    int iteration=0;
    if (input_complete_pb_path!="") {
        t.load_detatiled_mutations(input_complete_pb_path);
    }else{
        if (input_pb_path!="") {
            t=load_tree(input_pb_path, origin_states);
        }else {
            t=load_vcf_nh_directly(input_nh_path, input_vcf_path, origin_states);
        }
        puts("Checkpoint initial tree.\n");
        t.save_detailed_mutations(std::string(intermediate_pb_base_name).append(std::to_string(iteration++)).append(".pb"));
    }
    
    size_t new_score;
    size_t score_before;
    int stalled = 0;

    #ifndef NDEBUG
    Original_State_t origin_state_to_check(origin_states);
    check_samples(t.root, origin_state_to_check, &t);
    #endif

    score_before = t.get_parsimony_score();
    new_score = score_before;
    fprintf(stderr, "after state reassignment:%zu\n", score_before);

/*    t.breadth_first_expansion();
    t.depth_first_expansion();
    populate_mutated_pos(origin_state_to_check);
    clean_tree(t);
    fprintf(stderr, "nodes %zu\n",t.all_nodes.size());
    auto src=t.get_node("37308");
    auto dst=t.get_node("37313");
    output_t out;
    individual_move(src, dst, get_LCA(src,dst), out);
    //find_profitable_moves(src, out, 8);
    //Find nodes to search
    */
    tbb::concurrent_vector<MAT::Node *> nodes_to_search;
    std::vector<MAT::Node *> bfs_ordered_nodes;
    bfs_ordered_nodes = t.breadth_first_expansion();
    /*apply_moves(out.moves, t, bfs_ordered_nodes, nodes_to_search);
    fprintf(stderr, "%zu,nodes %zu\n",t.get_parsimony_score(),t.all_nodes.size());
    clean_tree(t);
    fprintf(stderr, "%zu,nodes %zu\n",t.get_parsimony_score(),t.all_nodes.size());*/
    size_t inner_loop_score_before = score_before;
    movalbe_src_log=fopen(profitable_src_log.c_str(),"w");
    if (!movalbe_src_log) {
        perror(("Error writing to log file "+profitable_src_log).c_str());
        movalbe_src_log=fopen("/dev/null", "w");
    }
    bool isfirst=true;
    while(stalled<=1){
    if (interrupted) {
        break;
    }
    bfs_ordered_nodes = t.breadth_first_expansion();
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search,isfirst,radius);
    isfirst=false;
    printf("%zu nodes to search\n",nodes_to_search.size());
    if (nodes_to_search.empty()) {
        break;
    }
    //Actual optimization loop
    while (!nodes_to_search.empty()) {
        if (interrupted) {
            break;
        }
        bfs_ordered_nodes = t.breadth_first_expansion();
        new_score =
            optimize_tree(bfs_ordered_nodes, nodes_to_search, t,radius,movalbe_src_log
            #ifndef NDEBUG
            , origin_states
            #endif
            );
        fprintf(stderr, "after optimizing:%zu\n", new_score);
        t.save_detailed_mutations(std::string(intermediate_pb_base_name).append(std::to_string(iteration++)).append(".pb"));
        if (new_score >= inner_loop_score_before) {
            stalled++;
        } else {
            inner_loop_score_before = new_score;
            stalled = 0;
        }

    }
        clean_tree(t);
    }
    fclose(movalbe_src_log);
    save_final_tree(t, origin_states, output_path);
    for(auto& pos:mutated_positions){
        delete pos.second;
    }
    t.delete_nodes();
    fprintf(stderr, "early stop savings:%d\n",early_stop_saving);
}