#include "common.hpp"
#include "../usher_graph.hpp"

typedef tbb::concurrent_unordered_map<std::string, std::string> concurMap;

po::variables_map parse_summary_command(po::parsed_options parsed);
bool consistent(MAT::Tree T, MAT::Tree B, concurMap& consistNodes);
bool chelper(MAT::Node* a, MAT::Node* b, concurMap& consistNodes);
void merge_main(po::parsed_options parsed);
MAT::Tree subtree(MAT::Tree tree, std::vector<std::string> common);
