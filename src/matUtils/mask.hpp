#include "common.hpp"

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void simplify_tree(MAT::Tree* T);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
void restrictMutationsLocally (std::string mutations_filename, MAT::Tree* T, bool global = false);
void renameSamples(std::string rename_filename, MAT::Tree& T);
void moveNodes (std::string node_filename, MAT::Tree* T);