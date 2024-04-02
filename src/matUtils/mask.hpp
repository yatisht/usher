#include "common.hpp"
#include <sstream>
#include <string>

po::variables_map parse_mask_command(po::parsed_options parsed);
void mask_main(po::parsed_options parsed);
void simplify_tree(MAT::Tree* T);
void restrictSamples (std::string samples_filename, MAT::Tree& T);
void restrictMutationsLocally (std::string mutations_filename, MAT::Tree* T, bool global = false);
void renameSamples(std::string rename_filename, MAT::Tree& T);
void moveNodes (std::string node_filename, MAT::Tree* T);
void localMask (uint32_t snp_distance, uint32_t max_snp_distance ,MAT::Tree& T, std::string diff_file);
std::map<std::string, std::map<int, int>> readDiff (const std::string& diff_file);
//void subtreeMask (Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data);
void dfs(Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data);
void dfsUtil(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* anc_node, std::map<std::string, std::map<int, int>>& diff_data, std::vector<std::tuple<int, int>>& missing, std::set<std::string>& visited);
void processNode(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* anc_node, std::vector<std::tuple<int, int>>& missing, std::map<std::string, std::map<int, int>>& diff_data);