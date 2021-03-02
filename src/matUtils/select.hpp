#include "common.hpp"

std::vector<std::string> read_sample_names (std::string sample_filename);
std::vector<std::string> get_clade_samples (MAT::Tree T, std::string clade_name);
std::vector<std::string> get_mutation_samples (MAT::Tree T, std::string mutation_id);
std::vector<std::string> get_parsimony_samples (MAT::Tree T, float max_parsimony);
std::vector<std::string> get_clade_representatives(MAT::Tree T);
std::vector<std::string> sample_intersect (std::vector<std::string> samples, std::vector<std::string> nsamples);