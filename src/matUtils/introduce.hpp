#include "common.hpp"
// #include "select.hpp"

po::variables_map parse_introduce_command(po::parsed_options parsed);
std::map<std::string, std::vector<std::string>> read_two_column (std::string sample_filename);
void record_clade_regions(MAT::Tree* T, std::map<std::string, std::map<std::string, float>> region_assignments, std::string filename);
size_t get_monophyletic_cladesize(MAT::Tree* T, std::map<std::string, float> assignments, MAT::Node* subroot);
float get_association_index(MAT::Tree* T, std::map<std::string, float> assignments, MAT::Node* subroot);
std::vector<std::string> find_introductions(MAT::Tree* T, std::vector<std::string> samples, std::vector<std::string> regions, bool add_info, std::string clade_output); 
std::map<std::string, float> get_assignments(MAT::Tree* T, std::unordered_set<std::string> sample_set);
void introduce_main(po::parsed_options parsed);