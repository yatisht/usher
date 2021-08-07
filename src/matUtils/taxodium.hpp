#include "common.hpp"
#include "taxodium.pb.h"

void save_taxodium_tree (MAT::Tree &tree, std::string filename, std::vector<std::unordered_map<std::string,std::unordered_map<std::string,std::string>>> &catmeta);
