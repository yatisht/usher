#include "common.hpp"

po::variables_map parse_prune_command(po::parsed_options parsed);
void prune_main(po::parsed_options parsed);
MAT::Tree prune_leaves (const MAT::Tree& T, std::vector<std::string> sample_names);
MAT::Tree get_sample_subtree (const MAT::Tree& T, std::vector<std::string> sample_names);
MAT::Tree get_sample_prune (const MAT::Tree& T, std::vector<std::string> sample_names);
