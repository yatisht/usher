#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#define NEW_MARK -2
struct state_change {
    int position;
    uint8_t chr_idx;
    uint8_t old_state;
    uint8_t new_state;
    state_change() : position(NEW_MARK) {}
    state_change(const MAT::Mutation &in, uint8_t old_state)
        : old_state(old_state) {
        position = in.get_position();
        chr_idx = in.get_chromIdx();
        new_state = in.get_mut_one_hot();
    }
};
typedef std::vector<state_change> State_Change_Collection;
struct Altered_Node_t {
    MAT::Node *altered_node;
    State_Change_Collection changed_states;
    Altered_Node_t() {}
    Altered_Node_t(MAT::Node *node) : altered_node(node) {}
};
#define SIBLING_UNIQUE_SHAMT 0
#define HAVE_SHARED_SHAMT 1
#define SIBLING_INCONSISTENT_SHAMT 2
#define NEW_NODE_INCONSISTENT_SHAMT 3
char merge_new_node_mutations(
    const MAT::Mutations_Collection &new_node_mutations,
    const MAT::Mutations_Collection &sibling_node_mutations,
    MAT::Mutations_Collection &shared_node_mutations_out,
    MAT::Mutations_Collection &sibling_node_mutations_out,
    MAT::Mutations_Collection &new_node_mutations_out,
    MAT::Mutations_Collection &to_merge_if_children);
MAT::Node *replace_with_internal_node(MAT::Node *to_replace,
                                             MAT::Tree &tree);