#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
//mark a state change collection is produced by merging in forward pass
#define NEW_MARK -2
#define NEW_SRC_MARK -3
//Represent a change in parent state
struct state_change {
    int position;
    uint8_t chr_idx;
    uint8_t old_state;
    uint8_t new_state;
    state_change() : position(NEW_MARK) {}
    state_change(int pos) : position(pos) {}
    state_change(const MAT::Mutation &in, uint8_t old_state)
        : old_state(old_state) {
        position = in.get_position();
        chr_idx = in.get_chromIdx();
        new_state = in.get_mut_one_hot();
    }
};
typedef std::vector<state_change> State_Change_Collection;
//used in backward pass to record nodes need to do backward pass, and their mut_nuc change
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
/**
 * @brief Get the mutation vector of the src node, dst node, and node from branch splitting
    If not LCA                             or                       LCA
        |                                                   /               \
        shared_node                               new_node(moved src)       sibling_node(src_branch_node)
    /                   \
new_node(moved src)     sibling node (dst)
 * @param new_node_mutations Mutation vector of src node to preserve its state after moving as sibling of sibling node
 * @param sibling_node_mutations the dst node if not moving to LCA, src_branch_node if moving to LCA
 * @param shared_node_mutations_out mutation vector of the node from spliting
 * @param sibling_node_mutations_out of sibling node
 * @param new_node_mutations_out of moved src node
 * @param to_merge_if_children if not to split branch, from no shared mutation, mutations need to be merged with new_node_mutations_out to perserve its state
 * @return bit flag with value SIBLING_UNIQUE: sibling node have mutation not found in new node
        HAVE_SHARED: sibling node and new node share some mutations
        SIBLING_INCONSISTENT (not used): sibling node need to be cleaned, some loci can change to follow paret state
        NEW_NODE_INCONSISTENT (not used): new node need to be cleaned, some loci can change to follow paret state
 */
char merge_new_node_mutations(
    const MAT::Mutations_Collection &new_node_mutations,
    const MAT::Mutations_Collection &sibling_node_mutations,
    MAT::Mutations_Collection &shared_node_mutations_out,
    MAT::Mutations_Collection &sibling_node_mutations_out,
    MAT::Mutations_Collection &new_node_mutations_out);
//insert a node between to_replace and its parent, then return it
MAT::Node *replace_with_internal_node(MAT::Node *to_replace,
                                             MAT::Tree &tree);
/**
 * @brief Move src to dst, while perserving the state of all nodes
 * @param src
 * @param dst
 * @param altered_node nodes whose children set have changed, need to do backward pass on them
 * @param tree
 * @param deleted internal node lost all their children, or nodes with only one children left, so delete
 * @param nodes_to_clean src node or src_branch_node, or dst node that have their parent state changed due to branch splitting
 */
void move_node(MAT::Node *src, MAT::Node *dst,
               std::vector<MAT::Node *> &altered_node, MAT::Tree &tree,
               std::unordered_set<size_t> &deleted,
               std::vector<MAT::Node *> &nodes_to_clean
);

#ifdef CHECK_STATE_REASSIGN
MAT::Tree reassign_state_full(MAT::Tree &tree_in);
void compare_mutation_tree(MAT::Tree &t,MAT::Tree &new_tree);
#endif
//update major allele set and boundary allele set of altered_nodes_in, output change in mut_nuc in nodes_with_changed_states_out(such changes are done with the assumption that the state of parent node is unchanged, but this assumption will be fixed in forward pass)
void reassign_backward_pass(
    const std::vector<MAT::Node *> &altered_nodes_in,
    std::vector<Altered_Node_t> &nodes_with_changed_states_out
#ifdef CHECK_STATE_REASSIGN
    ,
    MAT::Tree &new_tree
#endif
) ;
//reassign state of nodes in in to match with parent state if possible
void forward_pass(std::vector<Altered_Node_t> &in
#ifdef CHECK_STATE_REASSIGN
                  ,
                  MAT::Tree &new_tree
#endif
);
//merge mutation vector of nodes with their parent node deleted because it only have 1 child left
bool merge_mutation_single_child(MAT::Node *node,const MAT::Mutations_Collection &merge_with);
void clean_up_src_states(MAT::Node *src,std::vector<Altered_Node_t> &out);
