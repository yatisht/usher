#include "mutation_annotated_tree.hpp"
#include <cstddef>
#include <tbb/concurrent_vector.h>
struct backward_pass_range {
    union {
        size_t first_child_bfs_idx;
        const std::string* identifier;
    };
    size_t child_size;
};

struct forward_pass_range {
    size_t parent_bfs_idx;
    size_t child_size;
    size_t left_child_idx;
    size_t right_child_idx;
};
typedef std::vector<std::pair<long, nuc_one_hot>> mutated_t;
struct mutated_t_comparator {
    bool operator()(const std::pair<long, nuc_one_hot>& lhs,const std::pair<long, nuc_one_hot>& rhs) const {
        return lhs.first>rhs.first;
    }
};
void Fitch_Sankoff_prep(const std::vector<Mutation_Annotated_Tree::Node*>& bfs_ordered_nodes, std::vector<backward_pass_range>& child_idx_range,std::vector<forward_pass_range>& parent_idx);
void Fitch_Sankoff_Whole_Tree(const std::vector<backward_pass_range>& child_idx_range,const std::vector<forward_pass_range>& parent_idx,const Mutation_Annotated_Tree::Mutation & base,const mutated_t& mutated,std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>& output,Mutation_Annotated_Tree::Tree* try_similar=nullptr);
#if defined CHECK_STATE_REASSIGN||defined DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
void FS_backward_pass(const std::vector<Mutation_Annotated_Tree::Node*> bfs_ordered_nodes, std::vector<uint8_t>& boundary1_major_allele,const std::unordered_map<std::string, nuc_one_hot>& mutated,nuc_one_hot ref_nuc);
int FS_forward_assign_states_only(const std::vector<Mutation_Annotated_Tree::Node*>& bfs_ordered_nodes,const std::vector<uint8_t>& boundary1_major_allele,const nuc_one_hot parent_state,std::vector<uint8_t>& states_out,std::vector<std::vector<Mutation_Annotated_Tree::Node*>>& children_mutation_count);
#endif
void set_state_from_cnt(const std::array<int,4>& data, uint8_t& boundary1_major_allele_out);