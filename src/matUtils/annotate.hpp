#include "common.hpp"

po::variables_map parse_annotate_command(po::parsed_options parsed); 
void annotate_main(po::parsed_options parsed);
void assignLineages (MAT::Tree& T, const std::string& lineage_filename, float min_freq, float set_overlap, bool clear_current = false);
