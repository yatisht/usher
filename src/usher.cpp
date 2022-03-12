#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <memory>
#include <limits>
#include "boost/filesystem.hpp"
#include "usher_graph.hpp"
#include "parsimony.pb.h"
#include "version.hpp"
#include "usher_common.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

int main(int argc, char** argv) {
    //Variables to load command-line options using Boost program_options
    std::string tree_filename;
    std::string din_filename;
    std::string dout_filename;
    std::string outdir;
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    uint32_t max_trees;
    uint32_t max_uncertainty;
    uint32_t max_parsimony;;
    bool sort_before_placement_1 = false;
    bool sort_before_placement_2 = false;
    bool sort_before_placement_3 = false;
    bool reverse_sort = false;
    bool collapse_tree=false;
    bool collapse_output_tree=false;
    bool print_uncondensed_tree = false;
    bool print_parsimony_scores = false;
    bool retain_original_branch_len = false;
    bool no_add = false;
    bool detailed_clades = false;
    size_t print_subtrees_size=0;
    size_t print_subtrees_single=0;
    po::options_description desc{"Options"};

    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
    ("vcf,v", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
    ("tree,t", po::value<std::string>(&tree_filename)->default_value(""), "Input tree file")
    ("outdir,d", po::value<std::string>(&outdir)->default_value("."), "Output directory to dump output and log files [DEFAULT uses current directory]")
    ("load-mutation-annotated-tree,i", po::value<std::string>(&din_filename)->default_value(""), "Load mutation-annotated tree object")
    ("save-mutation-annotated-tree,o", po::value<std::string>(&dout_filename)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
    ("sort-before-placement-1,s", po::bool_switch(&sort_before_placement_1), \
     "Sort new samples based on computed parsimony score and then number of optimal placements before the actual placement [EXPERIMENTAL].")
    ("sort-before-placement-2,S", po::bool_switch(&sort_before_placement_2), \
     "Sort new samples based on the number of optimal placements and then the parsimony score before the actual placement [EXPERIMENTAL].")
    ("sort-before-placement-3,A", po::bool_switch(&sort_before_placement_3), \
     "Sort new samples based on the number of ambiguous bases [EXPERIMENTAL].")
    ("reverse-sort,r", po::bool_switch(&reverse_sort), \
     "Reverse the sorting order of sorting options (sort-before-placement-1 or sort-before-placement-2) [EXPERIMENTAL]")
    ("collapse-tree,c", po::bool_switch(&collapse_tree), \
     "Collapse internal nodes of the input tree with no mutations and condense identical sequences in polytomies into a single node and the save the tree to file condensed-tree.nh in outdir")
    ("collapse-output-tree,C", po::bool_switch(&collapse_output_tree), \
     "Collapse internal nodes of the output tree with no mutations before the saving the tree to file final-tree.nh in outdir")
    ("max-uncertainty-per-sample,e", po::value<uint32_t>(&max_uncertainty)->default_value(1e6), \
     "Maximum number of equally parsimonious placements allowed per sample beyond which the sample is ignored")
    ("max-parsimony-per-sample,E", po::value<uint32_t>(&max_parsimony)->default_value(1e6), \
     "Maximum parsimony score of the most parsimonious placement(s) allowed per sample beyond which the sample is ignored")
    ("write-uncondensed-final-tree,u", po::bool_switch(&print_uncondensed_tree), "Write the final tree in uncondensed format and save to file uncondensed-final-tree.nh in outdir")
    ("write-subtrees-size,k", po::value<size_t>(&print_subtrees_size)->default_value(0), \
     "Write minimum set of subtrees covering the newly added samples of size equal to this value")
    ("write-single-subtree,K", po::value<size_t>(&print_subtrees_single)->default_value(0), \
     "Similar to write-subtrees-size but produces a single subtree with all newly added samples along with random samples up to the value specified by this argument")
    ("write-parsimony-scores-per-node,p", po::bool_switch(&print_parsimony_scores), \
     "Write the parsimony scores for adding new samples at each existing node in the tree without modifying the tree in a file names parsimony-scores.tsv in outdir")
    ("multiple-placements,M", po::value<uint32_t>(&max_trees)->default_value(1), \
     "Create a new tree up to this limit for each possibility of parsimony-optimal placement")
    ("retain-input-branch-lengths,l", po::bool_switch(&retain_original_branch_len), \
     "Retain the branch lengths from the input tree in out newick files instead of using number of mutations for the branch lengths.")
    ("no-add,n", po::bool_switch(&no_add), \
     "Do not add new samples to the tree")
    ("detailed-clades,D", po::bool_switch(&detailed_clades), \
     "In clades.txt, write a histogram of annotated clades and counts across all equally parsimonious placements")
    ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
    ("version", "Print version number")
    ("help,h", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        if (vm.count("version")) {
            std::cout << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
        } else {
            std::cerr << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
            std::cerr << desc << std::endl;
        }
        // Return with error code 1 unless the user specifies help or version
        if(vm.count("help") || vm.count("version"))
            return 0;
        else
            return 1;
    }

    if (vm.count("version")) {
        std::cout << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
    }

    // timer object to be used to measure runtimes of individual stages
    Timer timer;

    fprintf(stderr, "Initializing %u worker threads.\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);

    MAT::Tree tmp_T;
    MAT::Tree* T = NULL;
    // Variables below used to store the different fields of the input VCF file
    std::vector<Missing_Sample> missing_samples;

    // Vectore to store the names of samples which have a high number of
    // parsimony-optimal placements
    std::vector<std::string> low_confidence_samples;

    // If tree filename is specified, UShER needs to first load the tree and
    // a create a new mutation-annotated tree object (stored in optimal_trees)
    // using the sample variants in the input VCF file. If the VCF contains
    // samples missing in the input tree, they get added to missing_samples
    if (tree_filename != "") {
        fprintf(stderr, "Loading input tree.\n");
        timer.Start();
        // Create a new tree from the input newick file
        tmp_T = MAT::create_tree_from_newick(tree_filename);
        if (tmp_T.root == NULL) {
            fprintf(stderr, "ERROR: Empty tree.\n");
            exit(1);
        }
        T = &tmp_T;

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        MAT::read_vcf(T, vcf_filename, missing_samples, true);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    // Check if input mutation-annotated tree filename is specified
    else if (din_filename != "") {
        timer.Start();

        fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", din_filename.c_str());

        // Load mutation-annotated tree
        tmp_T = MAT::load_mutation_annotated_tree(din_filename);

        if (tmp_T.root == NULL) {
            fprintf(stderr, "ERROR: Empty tree.\n");
            exit(1);
        }
        T = &tmp_T;
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        timer.Start();

        MAT::read_vcf(T, vcf_filename, missing_samples, false);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    } else {
        fprintf(stderr, "Error! No input tree or assignment file provided!\n");
        exit(1);
    }

    int return_val = usher_common(dout_filename, outdir, max_trees, max_uncertainty, max_parsimony,
                                  sort_before_placement_1, sort_before_placement_2, sort_before_placement_3, reverse_sort, collapse_tree,
                                  collapse_output_tree, print_uncondensed_tree, print_parsimony_scores, retain_original_branch_len, no_add,
                                  detailed_clades, print_subtrees_size, print_subtrees_single, missing_samples, low_confidence_samples, T);

    return return_val;
}
