#include "common.hpp"

po::variables_map parse_summary_command(po::parsed_options parsed);
void write_sample_table(MAT::Tree* T, std::string filename);
void write_clade_table(MAT::Tree* T, std::string filename);
void write_mutation_table(MAT::Tree* T, std::string filename);
void write_translate_table(MAT::Tree* T, std::string output_filename, std::string gtf_filename, std::string fasta_filename);
void write_aberrant_table(MAT::Tree* T, std::string filename);
std::map<std::set<std::string>,size_t> count_haplotypes(MAT::Tree* T);
void write_haplotype_table(MAT::Tree* T, std::string filename);
void write_sample_clades_table (MAT::Tree* T, std::string sample_clades);
void translate_main(MAT::Tree *T, std::string output_filename, std::string gff_filename, std::string fasta_filename );
void summary_main(po::parsed_options parsed);
