#include "common.hpp"

po::variables_map parse_describe_command(po::parsed_options parsed);
void describe_main(po::parsed_options parsed);
void mutation_paths(const MAT::Tree& T, std::string sample_filename);
