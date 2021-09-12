#include "common.hpp"
typedef tbb::concurrent_hash_map<std::string, std::string> concurMap;
po::variables_map parse_summary_command(po::parsed_options parsed);
bool consistent(MAT::Tree T, MAT::Tree B);
bool chelper(MAT::Node* a, MAT::Node* b);
void merge_main(po::parsed_options parsed);
extern concurMap consistNodes;
MAT::Tree subtree(MAT::Tree tree, std::vector<std::string> common);
struct placement_input {
    std::string missing_sample;
    MAT::Tree* T;
    MAT::Node* node;
    std::vector<MAT::Mutation>* missing_sample_mutations;

    int* best_set_difference;
    int* set_difference;
    size_t* best_node_num_leaves;
    size_t j;
    size_t* best_j;
    size_t distance;
    size_t* best_distance;
    size_t* num_best;
    MAT::Node** best_node;

    std::vector<bool>* node_has_unique;
    std::vector<size_t>* best_j_vec;

    bool* has_unique;

    std::vector<MAT::Mutation>* excess_mutations;
    std::vector<MAT::Mutation>* imputed_mutations;

    placement_input () {
        distance = 0;
        best_distance = &distance;
    }
};

void placement(placement_input& input, bool compute_parsimony_scores, bool compute_vecs);

