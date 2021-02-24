#include "common.hpp"

po::variables_map parse_filter_command(po::parsed_options parsed);
void filter_main(po::parsed_options parsed);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
