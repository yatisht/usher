#include "Profitable_Moves_Enumerators.hpp"
void get_intermediate_nodes_mutations(
    const MAT::Node *this_node,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change
// const New_Tie_Collection_t &parent_ties_in,
// New_Tie_Collection_t &decrease_if_match_parent_state_out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg>& debug,
    const std::vector<MAT::Node *>& node_stack
#endif
);
void get_parsimony_score_change_from_add(
    MAT::Node *parent_node,
    const Mutation_Count_Change_Collection &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations, bool is_terminal,
    // New_Tie_Collection_t &decrease_if_match_parent_state_out,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg>& debug,
    const std::vector<MAT::Node *>& node_stack
#endif
) ;
void get_parent_altered_remove(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    //,New_Tie_Collection_t &decrease_if_match_parent_state,
    const MAT::Node *src, int &parent_parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg>& debug,
    const std::vector<MAT::Node *>& node_stack
#endif
) ;
void get_LCA_mutation(
    const MAT::Node *LCA, MAT::Node* is_src_terminal, MAT::Node* is_dst_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    // const New_Tie_Collection_t &unresolved_ties_from_src,
    // const New_Tie_Collection_t &unresolved_ties_from_dst,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    // New_Tie_Collection_t &new_ties,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<LCA_merged_states> &dbg_states_in,
    std::vector<MAT::Node *> node_stack, const MAT::Node *src_branch,
    const MAT::Node *dst_branch
#endif
);
