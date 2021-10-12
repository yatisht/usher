#include "mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <chrono>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <unistd.h>
uint32_t num_threads;
std::chrono::time_point<std::chrono::steady_clock> last_save_time;
bool no_write_intermediate;
size_t max_queued_moves;
std::chrono::steady_clock::duration save_period;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
FILE* movalbe_src_log;
bool interrupted;
std::condition_variable progress_bar_cv;
bool timed_print_progress;
tbb::concurrent_unordered_map<MAT::Mutation, tbb::concurrent_unordered_map<std::string, nuc_one_hot>*,Mutation_Pos_Only_Hash,
    Mutation_Pos_Only_Comparator>
    mutated_positions;
#undef NDEBUG
#include <cassert>
int main(int argc,char** argv) {
    Original_State_t origin_states;
    tbb::task_scheduler_init init(20);
    Mutation_Annotated_Tree::Tree ori_tree=Mutation_Annotated_Tree::create_tree_from_newick(argv[1]);
    load_vcf_nh_directly(ori_tree, argv[2], origin_states);
    unlink("intermediate_mutations_test_out.pb");
    ori_tree.populate_ignored_range();
    auto save_start=std::chrono::steady_clock::now();
    ori_tree.save_detailed_mutations("intermediate_mutations_test_out.pb");
    auto save_end=std::chrono::steady_clock::now();
    fprintf(stderr, "Save took %ld seconds",std::chrono::duration_cast<std::chrono::seconds>(save_end-save_start).count());
    auto load_start=std::chrono::steady_clock::now();
    Mutation_Annotated_Tree::Tree loaded_tree;
    loaded_tree.load_detatiled_mutations("intermediate_mutations_test_out.pb");
    auto load_end=std::chrono::steady_clock::now();
    fprintf(stderr, "load took %ld seconds",std::chrono::duration_cast<std::chrono::seconds>(load_end-load_start).count());
    auto ori_tree_dfs=ori_tree.depth_first_expansion();
    auto loaded_tree_dfs=loaded_tree.depth_first_expansion();
    assert(ori_tree_dfs.size()==loaded_tree_dfs.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,loaded_tree_dfs.size()),
    [&](tbb::blocked_range<size_t>& in) {
        for (size_t node_idx=in.begin(); node_idx<in.end(); node_idx++) {
            const auto ori_tree_node=ori_tree_dfs[node_idx];
            const auto loaded_tree_node=loaded_tree_dfs[node_idx];
            assert(ori_tree_node->identifier==loaded_tree_node->identifier);
            assert(ori_tree_node->children.size()==loaded_tree_node->children.size());
            for (size_t child_idx=0; child_idx<ori_tree_node->children.size(); child_idx++) {
                assert(ori_tree_node->children[child_idx]->identifier==loaded_tree_node->children[child_idx]->identifier);
            }
            assert(ori_tree_node->ignore.size()==loaded_tree_node->ignore.size());
            const auto& ori_node_mutation=ori_tree_node->mutations;
            const auto& loaded_node_mutation=loaded_tree_node->mutations;
            assert(ori_node_mutation.size()==loaded_node_mutation.size());
            for (size_t mut_idx=0; mut_idx<ori_node_mutation.size(); mut_idx++) {
                assert(ori_node_mutation[mut_idx]==loaded_node_mutation[mut_idx]);
            }
        }
    });
    ori_tree.delete_nodes();
    loaded_tree.delete_nodes();
}