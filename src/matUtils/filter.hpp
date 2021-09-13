#include "common.hpp"

MAT::Tree filter_master(const MAT::Tree& T, std::vector<std::string> sample_names, bool prune, bool keep_clade_annotations=false);
MAT::Tree prune_leaves (const MAT::Tree& T, std::vector<std::string> sample_names);
MAT::Tree get_sample_subtree (const MAT::Tree& T, std::vector<std::string> sample_names, bool keep_clade_annotations=false);
MAT::Tree get_sample_prune (const MAT::Tree& T, std::vector<std::string> sample_names, bool keep_clade_annotations=false);
void resolve_polytomy(MAT::Tree &T, std::vector<MAT::Node*> polytomy_nodes);
MAT::Tree resolve_all_polytomies(MAT::Tree T);
void reroot_tree(MAT::Tree* T, std::string rnid);