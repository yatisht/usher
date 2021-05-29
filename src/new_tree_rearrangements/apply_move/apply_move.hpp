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
void move_node(MAT::Node *src, MAT::Node *dst,
               std::vector<MAT::Node *> &altered_node, MAT::Tree &tree,
               std::unordered_set<size_t> &deleted,
               std::vector<MAT::Node *> &nodes_to_clean
);

#ifdef CHECK_STATE_REASSIGN
MAT::Tree reassign_state_full(MAT::Tree &tree_in);
void compare_mutation_tree(MAT::Tree &t,MAT::Tree &new_tree);
#endif

void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
    ,
    MAT::Tree &new_tree
#endif
) ;

void forward_pass(std::vector<Altered_Node_t> &in
#ifdef CHECK_STATE_REASSIGN
                  ,
                  MAT::Tree &new_tree
#endif
);
void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
    ,
    MAT::Tree &new_tree
#endif
);
bool
merge_mutation_single_child(MAT::Node *node,
                            const MAT::Mutations_Collection &merge_with);
                            void clean_up_src_states(MAT::Node *src,
                                std::vector<Altered_Node_t> &out);