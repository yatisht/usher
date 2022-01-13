#include "src/matOptimize/check_samples.hpp"
#include "usher.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <cstdio>
#include <tbb/task_scheduler_init.h>
#include <utility>
namespace po = boost::program_options;

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

    #ifndef NDEBUG
    std::vector<Sampled_Tree_Node *> output;
    Sampled_Tree_Node *sampled_tree_root=sample_tree(tree, sampling_radius);
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(tree,output,sampling_radius);
    #endif
    //std::random_shuffle(samples_to_place.begin(),samples_to_place.end());
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
    MAT::save_mutation_annotated_tree(tree, protobuf_out);

}