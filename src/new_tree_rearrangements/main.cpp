#include "Fitch_Sankoff.hpp"
#include "check_samples.hpp"
#include "priority_conflict_resolver.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <cstdio>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <unordered_set>
#include <boost/program_options.hpp> 
#include <vector>
#include <iostream>
namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
static void save_final_tree(MAT::Tree &t, Original_State_t origin_states,
                            const std::string &output_path) {
    std::vector<MAT::Node *> dfs = t.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
                      [&dfs](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              dfs[i]->mutations.remove_invalid();
                          }
                      });
    fix_condensed_nodes(&t);
    fprintf(stderr, "%d condensed_nodes",t.condensed_leaves.size());
    check_samples(t.root, origin_states, &t);
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, output_path);
}
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst){
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src=src->parent;
        }
        else if (src->dfs_index<dst->dfs_index) {
            dst=dst->parent;
        }
    }
    return src;
}

static size_t
optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
              tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
              MAT::Tree &t, Original_State_t origin_states,int radius) {
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", bfs_ordered_nodes.size());
    fprintf(stderr, "Internal node size %zu\n", t.curr_internal_node);

    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    Conflict_Resolver resolver(bfs_ordered_nodes.size());
    output_t out;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                      [&nodes_to_search, &resolver,
                       &deferred_nodes,radius](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              output_t out;
                              find_profitable_moves(nodes_to_search[i], out, radius);
                              if (!out.moves.empty()) {
                                  deferred_nodes.push_back(
                                      out.moves[0]->get_src());
                                  resolver(out.moves);
                              }
                          }
                      });
    std::vector<Profitable_Moves_ptr_t> all_moves;
    resolver.schedule_moves(all_moves);
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
    ,origin_states
#endif
    );
    check_samples(t.root, origin_states, &t);
    nodes_to_search = std::move(deferred_nodes);
    return t.get_parsimony_score();
}

int main(int argc, char **argv) {
    std::string output_path;
    std::string input_pb_path;
    std::string input_nh_path;
    std::string input_vcf_path;
    std::string intermediate_pb_base_name;
    int radius=10;
    po::options_description desc{"Options"};

    desc.add_options()
        ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
        ("tree,t", po::value<std::string>(&input_nh_path)->default_value(""), "Input tree file")
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

    Original_State_t origin_states;
    Mutation_Annotated_Tree::Tree t=(input_pb_path!="")?load_tree(input_pb_path, origin_states):load_vcf_nh_directly(input_nh_path, input_vcf_path, origin_states);
    int iteration=0;
    save_intermediate_tree(t, std::string(intermediate_pb_base_name).append(std::to_string(iteration++)).append(".pb"));
    Original_State_t origin_state_to_check(origin_states);
    check_samples(t.root, origin_state_to_check, &t);
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

    //MAT::Node* src=t.get_node("6360");
    //MAT::Node* dst=t.get_node("6347");
    //MAT::Node* LCA=get_LCA(src, dst);
    //individual_move(src,dst,LCA);
    while(stalled<=1){
    bfs_ordered_nodes = t.breadth_first_expansion();
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search);
    while (!nodes_to_search.empty()) {
        new_score =
            optimize_tree(bfs_ordered_nodes, nodes_to_search, t, origin_states,radius);
        fprintf(stderr, "after optimizing:%zu\n", new_score);
        save_intermediate_tree(t,std::string(intermediate_pb_base_name).append(std::to_string(iteration++)).append(".pb"));
        if (new_score >= inner_loop_score_before) {
            stalled++;
        } else {
            inner_loop_score_before = new_score;
            stalled = 0;
        }
    }}
    save_final_tree(t, origin_states, output_path);
}