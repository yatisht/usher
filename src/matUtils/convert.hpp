#include "common.hpp"

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes);

void make_json (MAT::Tree T, std::string json_filename);
MAT::Tree load_mat_from_json(std::string json_filename);