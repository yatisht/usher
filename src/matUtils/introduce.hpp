#include "common.hpp"
#include "select.hpp"

po::variables_map parse_introduce_command(po::parsed_options parsed);
std::vector<std::string> find_introductions(MAT::Tree* T, std::vector<std::string> samples); 
void introduce_main(po::parsed_options parsed);