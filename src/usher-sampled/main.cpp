#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <chrono>
#include <climits>
#include <complex>
#include <csignal>
#include <cstdio>
#include <random>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <utility>
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
void fix_condensed_nodes(MAT::Tree *tree) {
    std::vector<MAT::Node *> nodes_to_fix;
    for (auto iter : tree->all_nodes) {
        if (tree->condensed_nodes.count(iter.first) &&
                (!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for (auto node : nodes_to_fix) {
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier,
                          std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier, node);
    }
}
static int set_descendant_count(MAT::Node* root){
    size_t child_count=0;
    for (auto child : root->children) {
        child_count+=set_descendant_count(child);
    }
    root->bfs_index=child_count;
    return child_count;
}

int main(int argc, char **argv) {
    std::string vcf_filename;
    std::string protobuf_in;
    std::string protobuf_out;
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    int first_n_sample;
    int sampling_radius;
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
    fprintf(stderr, "Sampling at radius %d \n",sampling_radius);
    auto start_time=std::chrono::steady_clock::now();
    tbb::task_scheduler_init init(num_threads);
    MAT::Tree tree=MAT::load_mutation_annotated_tree(protobuf_in);
    tree.uncondense_leaves();
    std::vector<Sample_Muts> samples_to_place;
    Sample_Input(vcf_filename.c_str(),samples_to_place,tree);
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
        node->children.reserve(4*node->children.size());
        #endif
    }
    assign_descendant_muts(tree);
    assign_levels(tree.root);
    //std::minstd_rand rng(5);
    //std::shuffle(samples_to_place.begin(),samples_to_place.end(),rng);
    /*for (auto && to_place : samples_to_place) {
        if (to_place.sample_name=="s2138638s") {
            place_sample(std::move(to_place),sampled_tree_root
                  , tree,
                  sampling_radius
#ifndef NDEBUG
                  ,
                  ori_state
#endif
);
        }
    }
    for (auto && to_place : samples_to_place) {
        if (to_place.sample_name=="s1429085s") {
            place_sample(std::move(to_place),sampled_tree_root
                  , tree,
                  sampling_radius
#ifndef NDEBUG
                  ,
                  ori_state
#endif
);
        }
    }
    for (auto && to_place : samples_to_place) {
        if (to_place.sample_name=="s483795s") {
            place_sample(std::move(to_place),sampled_tree_root
                  , tree,
                  sampling_radius
#ifndef NDEBUG
                  ,
                  ori_state
#endif
);
        }
    }*/
    set_descendant_count(tree.root);
    place_sample(samples_to_place,tree,4
    #ifndef NDEBUG
        ,ori_state
    #endif
    );
    dfs=tree.depth_first_expansion();
    clean_up_leaf(dfs);
    fix_condensed_nodes(&tree);
    MAT::save_mutation_annotated_tree(tree, protobuf_out);
    auto duration=std::chrono::steady_clock::now()-start_time;
    fprintf(stderr, "Took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
}