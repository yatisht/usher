#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <cstdio>
#include <string>

#include <tbb/task_scheduler_init.h>
#include <unordered_set>
#include <boost/program_options.hpp> 
#include <vector>
#include <iostream>
namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

int main(int argc, char **argv) {
    std::string output_path;
    std::string input_pb_path;
    std::string input_nh_path;
    std::string input_vcf_path;
    std::string intermediate_pb_base_name;
    int radius=10;
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    desc.add_options()
        ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
        ("tree,t", po::value<std::string>(&input_nh_path)->default_value(""), "Input tree file")
        ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("load-mutation-annotated-tree,i", po::value<std::string>(&input_pb_path)->default_value(""), "Load mutation-annotated tree object")
        ("save-mutation-annotated-tree,o", po::value<std::string>(&output_path)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
        ("save-intermediate-mutation-annotated-tree,m", po::value<std::string>(&intermediate_pb_base_name)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
        ("radius,r", po::value<uint32_t>()->default_value(10),
         "Radius in which to restrict the SPR moves.");

        
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

    tbb::task_scheduler_init init(num_threads);

    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t=(input_pb_path!="")?load_tree(input_pb_path, origin_states):load_vcf_nh_directly(input_nh_path, input_vcf_path, origin_states);
    int iteration=0;
    puts("Checkpoint initial tree.\n");
    t.save_detailed_mutations(std::string(intermediate_pb_base_name).append(std::to_string(iteration++)).append(".pb"));
    #ifndef NDEBUG
    Original_State_t origin_state_to_check(origin_states);
    check_samples(t.root, origin_state_to_check, &t);
    #endif
    size_t score_before = t.get_parsimony_score();
    size_t new_score = score_before;
    fprintf(stderr, "after state reassignment:%zu\n", score_before);
    int stalled = 0;
    tbb::concurrent_vector<MAT::Node *> nodes_to_search;
    std::vector<MAT::Node *> bfs_ordered_nodes;
    bfs_ordered_nodes = t.breadth_first_expansion();
    for(auto node:bfs_ordered_nodes){
        node->tree=&t;
    }
    size_t inner_loop_score_before = score_before;
    /*std::vector<MAT::Node *> old_nodes = t.breadth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,old_nodes.size()),[&old_nodes](const tbb::blocked_range<size_t>&r){
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            auto node=old_nodes[idx];
            assert(node->is_root() || node->is_leaf() || node->children.size() > 1);
            if (node->parent){
            for(const auto mut:node->mutations){
            auto& par_mutations=node->parent->mutations;
            auto iter=par_mutations.find(mut.get_position());
            if(iter!=par_mutations.end()){
                assert(iter->get_par_one_hot()!=mut.get_mut_one_hot()||(!mut.is_valid()));
            }
        }}
        }
    });
    MAT::Node* src=t.get_node("MT834209.1|USA/WA-S1536/2020|20-05-25");
    MAT::Node* dst=t.get_node("6379");
    MAT::Node* LCA=get_LCA(src, dst);
    individual_move(src,dst,LCA);*/
    FILE* log=fopen("profitable_src","w");
    perror("");
    bool isfirst=true;
    while(stalled<=1){
    bfs_ordered_nodes = t.breadth_first_expansion();
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search,isfirst,radius);
    isfirst=false;
    printf("%zu nodes to search\n",nodes_to_search.size());
    if (nodes_to_search.empty()) {
        break;
    }
    while (!nodes_to_search.empty()) {
        bfs_ordered_nodes = t.breadth_first_expansion();
        new_score =
            optimize_tree(bfs_ordered_nodes, nodes_to_search, t,radius,log
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
        std::unordered_set<std::string> changed_nodes;
        clean_tree(t, changed_nodes);
    }
    fclose(log);
    save_final_tree(t, origin_states, output_path);
    for(auto& pos:mutated_positions){
        delete pos.second;
    }
    t.delete_nodes();
}