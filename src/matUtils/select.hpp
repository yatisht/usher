#include "common.hpp"
#include <regex>

std::vector<std::string> read_sample_names (std::string sample_filename);
std::vector<std::string> get_clade_samples (MAT::Tree* T, std::string clade_name);
std::vector<std::string> get_mutation_samples (MAT::Tree* T, std::string mutation_id);
std::vector<std::string> get_parsimony_samples (MAT::Tree* T, std::vector<std::string> samples_to_check, int max_parsimony);
std::vector<std::string> get_clade_representatives(MAT::Tree* T, size_t samples_per_clade);
std::vector<std::string> sample_intersect (std::unordered_set<std::string> samples, std::vector<std::string> nsamples);
std::vector<std::string> get_nearby (MAT::Tree* T, std::string sample_id, int number_to_get);
std::vector<std::string> get_short_steppers(MAT::Tree* T, std::vector<std::string> samples_to_check, int max_mutations);
std::vector<std::string> get_short_paths(MAT::Tree* T, std::vector<std::string> samples_to_check, int max_path);
std::unordered_map<std::string,std::unordered_map<std::string,std::string>> read_metafile(std::string metainf, std::set<std::string> samples_to_use);
std::vector<std::string> get_sample_match(MAT::Tree* T, std::vector<std::string> samples_to_check, std::string substring);
std::vector<std::string> fill_random_samples(MAT::Tree* T, std::vector<std::string> current_samples, size_t target_size, bool lca_limit = false);
std::vector<std::string> get_mrca_samples(MAT::Tree* T, std::vector<std::string> current_samples);
std::pair<std::vector<std::string>, size_t> get_closest_samples(MAT::Tree* T, std::string nid);
