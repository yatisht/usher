#include "common.hpp"

po::variables_map parse_annotate_command(po::parsed_options parsed);
void annotate_main(po::parsed_options parsed);
void assignLineages (MAT::Tree& T, const std::string& lineage_filename, bool clear_current = false);
void assignLineages (MAT::Tree& T, const std::string& clade_filename, const std::string& clade_mutations_filename, const std::string& clade_paths_filename, float min_freq, float mask_freq, float set_overlap, float clip_sample_frequency, bool clear_current = false, const std::string& mutations_filename = "", const std::string& details_filename = "");
void assignLineagesFromPaths (MAT::Tree& T, const std::string& clade_paths_filename, std::unordered_set<std::string>& clades_already_assigned);
