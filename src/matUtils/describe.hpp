#include "common.hpp"

po::variables_map parse_describe_command(po::parsed_options parsed);
void describe_main(po::parsed_options parsed);
std::vector<std::string> mutation_paths(const MAT::Tree& T, std::vector<std::string> samples);
std::vector<std::string> clade_paths(MAT::Tree T, std::vector<std::string> clades);