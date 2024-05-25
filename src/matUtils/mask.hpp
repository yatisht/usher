#include "common.hpp"
#include "maskselect.hpp"
#include <sstream>
#include <string>

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void simplify_tree(MAT::Tree* T);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
void restrictMutationsLocally (std::string mutations_filename, MAT::Tree* T, bool global = false);
void renameSamples(std::string rename_filename, MAT::Tree& T);
void moveNodes (std::string node_filename, MAT::Tree* T);
void localMask (uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file, std::string filename, uint32_t num_threads);
std::map<std::string, std::map<int, int>> readDiff (const std::string& diff_file);
//void mask_closest_samples_dfs(MAT::Node *node, MAT::Node *target, size_t path_length, size_t max_path_length, std::vector<std::pair<MAT::Node *, size_t>> &leaves, bool fixed_k);
void nodeComp(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::map<std::string, std::map<int, int>>& diff_data);
void getDistance(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data);

//std::pair<std::vector<std::string>, size_t> mask_get_closest_samples(MAT::Tree* T, std::string nid, bool fixed_k, size_t k);
//void subtreeMask (Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data);
void dfs(MAT::Node* l, int bl,Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data, uint32_t snp_distance, std::set<std::string>& visited, std::set<std::string>& leaves);
void dfsUtil(MAT::Node* l, int bl, Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data, std::set<std::string>& visited, uint32_t snp_distance, std::set<std::string>& leaves);
void processNode(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::map<std::string, std::map<int, int>>& diff_data);