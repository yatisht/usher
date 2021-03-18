#include "../tree_rearrangement.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <atomic>
#include <boost/program_options.hpp> 
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <string>
uint32_t num_cores;
extern std::atomic_bool interrupted;

static void catch_interrupt(int signal){
    fprintf(stderr, "interrupt recieved, finishing off \n");
    interrupted.store(true);
}

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
int main(int argc,char** argv){
    std::string din_filename;
    std::string dout_filename;
    int radius;
    int min_batch;
    int conflict_pct;
    
    num_cores = tbb::task_scheduler_init::default_num_threads();
    po::options_description desc{"Options"};
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
    ("load-mutation-annotated-tree,i", po::value<std::string>(&din_filename)->default_value(""), "Load mutation-annotated tree object")
    ("save-mutation-annotated-tree,o", po::value<std::string>(&dout_filename)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
    ("radius,r", po::value<int>(&radius)->default_value(2), "radius to search for new place to move a node")
    ("min_batch,m", po::value<int>(&min_batch)->default_value(50), "minimum batch size before restart the search")
    ("conflict_pct,c", po::value<int>(&conflict_pct)->default_value(20), "The percentage of conflicting moves before restarting the search");

    po::options_description all_options;
    all_options.add(desc);

    po::variables_map vm;
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

    if (din_filename=="") {
        fprintf(stderr, "Input MAT need to be specified.");
        exit(EXIT_FAILURE);
    }
    
    
    if (dout_filename=="") {
        fprintf(stderr, "Output MAT base name need to be specified.");
        exit(EXIT_FAILURE);
    }
    
    auto tmp_T = MAT::load_mutation_annotated_tree(din_filename);
    if(signal(SIGINT, catch_interrupt)==SIG_ERR){
        fputs( "failed to register interrupt signal handler",stderr);
    }
    Tree_Rearrangement::refine_trees(tmp_T,radius,min_batch,conflict_pct,dout_filename);
}