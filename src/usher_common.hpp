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

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

int usher_common(std::string dout_filename, std::string outdir, uint32_t num_threads, uint32_t max_trees,
                 uint32_t max_uncertainty, uint32_t max_parsimony, bool sort_before_placement_1, bool sort_before_placement_2, bool sort_before_placement_3,
                 bool reverse_sort, bool collapse_tree, bool collapse_output_tree, bool print_uncondensed_tree, bool print_parsimony_scores,
                 bool retain_original_branch_len, bool no_add, bool detailed_clades, size_t print_subtrees_size, size_t print_subtrees_single,
                 std::vector<Missing_Sample>& missing_samples, std::vector<std::string>& low_confidence_samples, MAT::Tree* loaded_MAT);
