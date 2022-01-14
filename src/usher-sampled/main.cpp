#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <chrono>
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
int main(int argc, char **argv) {
    std::string vcf_filename;
    std::string protobuf_in;
    std::string protobuf_out;
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
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
        ("Sampling radius,r",po::value(&sampling_radius)->default_value(10))
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
    #ifndef NDEBUG
    check_samples(tree.root, ori_state, &tree);
    fprintf(stderr, "\n------\n%zu samples\n",ori_state.size());
    #endif
    Sampled_Tree_Node *sampled_tree_root=sample_tree(tree, sampling_radius);
    
    #ifndef NDEBUG
    std::vector<Sampled_Tree_Node *> output;
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(tree,output,sampling_radius);
    #endif
    std::minstd_rand rng(0);
    std::shuffle(samples_to_place.begin(),samples_to_place.end(),rng);
    /*for (auto && to_place : samples_to_place) {
        if (to_place.sample_name=="s1785550s") {
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
    for (auto&& to_place : samples_to_place) {
        fprintf(stderr, "placing sample %s\n",to_place.sample_name.c_str());
        place_sample(std::move(to_place),sampled_tree_root
                  , tree,
                  sampling_radius
#ifndef NDEBUG
                  ,
                  ori_state
#endif
);
    }
    dfs=tree.depth_first_expansion();
    clean_up_leaf(dfs);
    MAT::save_mutation_annotated_tree(tree, protobuf_out);
    auto duration=std::chrono::steady_clock::now()-start_time;
    fprintf(stderr, "Took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
}