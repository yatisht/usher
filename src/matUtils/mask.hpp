#include "common.hpp"

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
