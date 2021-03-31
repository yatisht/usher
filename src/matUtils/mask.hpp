#include "common.hpp"

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void create_pseudosamples(MAT::Tree* T, size_t min_terminals);
void simplify_tree(MAT::Tree* T);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
