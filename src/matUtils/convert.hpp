#include "common.hpp"

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes);
void write_json_from_mat(MAT::Tree* T, std::string output_filename, std::map<std::string,std::map<std::string,std::string>> catmeta);
MAT::Tree load_mat_from_json(std::string json_filename);