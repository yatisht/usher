#include "common.hpp"
po::variables_map parse_summary_command(po::parsed_options parsed);
bool consistent(MAT::Tree T, MAT::Tree B);
bool chelper(MAT::Node* a, MAT::Node* b);
void merge_main(po::parsed_options parsed);
