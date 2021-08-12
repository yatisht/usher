#include "common.hpp"
#include "translate.hpp"
#include <google/protobuf/text_format.h> //TODO: remove (for debugging)

void save_taxodium_tree (MAT::Tree &tree, std::string out_filename, std::vector<std::string> meta_filenames, std::string gtf_filename, std::string fasta_filename);
std::unordered_map<std::string, std::vector<std::string>> read_metafiles_tax(std::vector<std::string> filenames, Taxodium::AllData &all_data);

