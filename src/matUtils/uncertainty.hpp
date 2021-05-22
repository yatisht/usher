#include "common.hpp"

std::vector<MAT::Node*> get_common_nodes (std::vector<std::vector<MAT::Node*>> nodepaths);
std::vector<float> get_all_distances(MAT::Node* target, std::vector<std::vector<MAT::Node*>> paths);
size_t get_neighborhood_size(std::vector<MAT::Node*> nodes, MAT::Tree* T);
void findEPPs (MAT::Tree* T, MAT::Node* node, bool get_nsize, size_t* nbest, size_t* nsize);
void findEPPs_wrapper (MAT::Tree Tobj, std::string sample_file, std::string fepps, std::string fneigh);
std::vector<std::string> get_samples_epps (MAT::Tree* T, size_t max_epps, std::vector<std::string> to_check);
po::variables_map parse_uncertainty_command(po::parsed_options parsed);
void uncertainty_main(po::parsed_options parsed);