#include "common.hpp"
#include "translate.hpp"

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes, std::vector<std::string> samples_vec = {});
void write_json_from_mat(MAT::Tree* T, std::string output_filename, std::vector<std::map<std::string,std::map<std::string,std::string>>>* catmeta);
MAT::Tree load_mat_from_json(std::string json_filename);
void get_minimum_subtrees(MAT::Tree* T, std::vector<std::string> samples, size_t target_size, std::string output_dir, std::vector<std::map<std::string,std::map<std::string,std::string>>>* catmeta, std::string json_n, std::string newick_n, bool retain_original_branch_len = false);
std::vector<std::string> get_nearby (MAT::Tree* T, std::string sample_id, int number_to_get);
void save_taxodium_tree (MAT::Tree &tree, std::string out_filename, std::vector<std::string> meta_filenames, std::string gtf_filename, std::string fasta_filenames, std::string title, std::string description);
std::unordered_map<std::string, std::vector<std::string>> read_metafiles_tax(std::vector<std::string> filenames, Taxodium::AllData &all_data);
