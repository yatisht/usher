#include "common.hpp"

po::variables_map parse_convert_command(po::parsed_options parsed);
void convert_main(po::parsed_options parsed);
void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes);
