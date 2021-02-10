#include "tree_rearrangement_internal.hpp"
Move* check_mut_conflict(MAT::Node* node, const std::vector<Fitch_Sankoff_Result*>& states);
bool check_loop_conflict(MAT::Node* src, MAT::Node* dst);
void register_mut_conflict(MAT::Node* node, const std::vector<Fitch_Sankoff_Result*>& states);
void register_loop_conflict(MAT::Node* src, MAT::Node* dst);

std::unordered_set<MAT::Node*> crossing_paths;
std::unordered_map<MAT::Node*, std::unordered_map<int,Move*>> subtree_moved_mutations;

