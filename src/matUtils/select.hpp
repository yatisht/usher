#include "common.hpp"

std::vector<std::string> read_sample_names (std::string sample_filename);
std::vector<std::string> get_clade_samples (const MAT::Tree& T, std::string clade_name);
std::vector<std::string> get_mutation_samples (const MAT::Tree& T, std::string mutation_id);
std::vector<std::string> get_samples_epps (const MAT::Tree& T, size_t max_epps);