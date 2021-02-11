#include "common.hpp"

po::variables_map parse_prune_command(po::parsed_options parsed);
void prune_main(po::parsed_options parsed);
MAT::Tree prune_leaves (const MAT::Tree& T, std::string sample_filename);
MAT::Tree prune_all_but_leaves (const MAT::Tree& T, std::string sample_filename);
