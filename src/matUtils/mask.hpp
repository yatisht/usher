#include "common.hpp"
#include <sstream>
#include <string>
#include <future>

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void simplify_tree(MAT::Tree* T);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
void restrictMutationsLocally (std::string mutations_filename, MAT::Tree* T, bool global = false);
void renameSamples(std::string rename_filename, MAT::Tree& T);
void moveNodes (std::string node_filename, MAT::Tree* T);
void localMask (uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file, std::string filename, uint32_t num_threads);
std::map<std::string, std::map<int, int>> readDiff (const std::string& diff_file);
//void traverseToMRCA(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data);
//void traverseFromLeafParentToMRCA(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data);

void nodeComp(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data);
bool prev_check(std::pair<int, int>& prev, std::pair<int, int> line);
void combine_missing(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::list<std::pair<int, int>>& missing_data, std::map<std::string, std::map<int, int>>& diff_data); 
void getDistance(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::map<std::string, std::set<std::string>>& comparisons);
