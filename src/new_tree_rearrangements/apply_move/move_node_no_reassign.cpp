#include "src/new_tree_rearrangements/check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include "apply_move.hpp"
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef CHECK_PRIMARY_MOVE
static void get_mutation_set_from_root(MAT::Node *node,
                                       MAT::Mutations_Collection &mutations) {
    MAT::Mutations_Collection temp;
    while (node) {
        temp.merge(node->mutations, MAT::Mutations_Collection::KEEP_SELF);
        node = node->parent;
    }
    for (auto &m : temp) {
        if (m.get_ref_one_hot() != m.get_mut_one_hot()) {
            mutations.push_back(m);
        }
    }
}
#endif

bool get_new_mut_binary(MAT::Mutation &base, nuc_one_hot left_branch,
                        nuc_one_hot right_branch);
//Because moving nodes around breaks dfs order, and doing a dfs after each move is too expensive, so 
// use a slower method of tracing back from src to root, and find the earliest intersection with path
// from dst to root
static bool find_path_no_dfs(std::vector<MAT::Node *> &dst_to_root_path,
                      std::vector<MAT::Node *> &src_to_root_path,
                      MAT::Node *src_ancestor, MAT::Node *dst_ancestor) {
    std::unordered_map<size_t, int> dst_to_root_idx_map;
    int idx=0;
    while (dst_ancestor) {
        if (dst_ancestor==src_ancestor) {
            return false;
        }
        dst_to_root_path.push_back(dst_ancestor);
        dst_to_root_idx_map.emplace((size_t)dst_ancestor,idx);
        idx++;
        dst_ancestor=dst_ancestor->parent;
    }
    while (src_ancestor) {
        auto iter=dst_to_root_idx_map.find((size_t)src_ancestor);
        if (iter!=dst_to_root_idx_map.end()) {
            dst_to_root_path.erase(dst_to_root_path.begin()+iter->second,dst_to_root_path.end());
            return true;
        }
        src_to_root_path.push_back(src_ancestor);
        src_ancestor=src_ancestor->parent;
    }
    assert(false);
    return false;
}

bool merge_mutation_single_child(MAT::Node *node,
                            const MAT::Mutations_Collection &merge_with) {
    bool have_inconsistent_mut = false;
    auto iter = merge_with.begin();
    auto end = merge_with.end();
    MAT::Mutations_Collection mutations;
    for (const auto &mut : node->mutations) {
        while (iter != end && iter->get_position() < mut.get_position()) {
            //mutation in parent of node, need to be merged to mutation of node to remove the parent of node
            if (iter->is_valid()) {
                mutations.push_back(*iter);
                //set left and right children both to the parent state as they inherits
                mutations.back().set_children(0, iter->get_mut_one_hot(),
                                              iter->get_mut_one_hot(),
                                              iter->get_mut_one_hot());
            }
            iter++;
        }
        mutations.push_back(mut);
        auto &new_mut = mutations.back();
        if (iter != end && iter->get_position() == mut.get_position()) {
            //match, adjust parent state
            new_mut.set_par_one_hot(iter->get_par_one_hot());
            have_inconsistent_mut |= (mut.is_valid() ^ new_mut.is_valid());
            iter++;
        }
        //remove unnecessary mutation
        if (new_mut.get_all_major_allele() == new_mut.get_par_one_hot() &&
            (!new_mut.get_boundary1_one_hot())) {
            mutations.mutations.pop_back();
        }
    }
    while (iter != end) {
        if (iter->is_valid()) {
            mutations.push_back(*iter);
            mutations.back().set_children(0, iter->get_mut_one_hot(),
                                          iter->get_mut_one_hot(),
                                          iter->get_mut_one_hot());
        }
        iter++;
    }
    node->mutations.swap(mutations);
    return have_inconsistent_mut;
}
/**
 * @brief Remove internal nodes that lost all its children or having only one child after removing src
 * @param[in] node parent of src node just removed
 * @param[out] deleted record nodes to free
 * @param[out] nodes_to_clean add sibling of src , if any, 
    a forward pass on them is needed, as their parent state have changed
 * @param[out] tree to remove removed node from all_nodes hashmap
 * @return the last surviving parent of src, need to do backward pass on it as its children set have changed
 */
static MAT::Node *
clean_up_after_remove(MAT::Node *node, std::unordered_set<size_t> &deleted,
                      std::vector<MAT::Node *> &nodes_to_clean,MAT::Tree& tree) {
    MAT::Node *parent_node = node->parent;
    //root node, cannot remove
    if (!parent_node) {
        return node;
    }
    //internal node lost all children 
    if (node->children.empty()) {
        auto &parent_children = parent_node->children;
        auto iter =
            std::find(parent_children.begin(), parent_children.end(), node);
        parent_children.erase(iter);
        deleted.insert((size_t)node);
        tree.all_nodes.erase(node->identifier);
        //delete node;
        //this node is removed, its parental node may have similar situation
        return clean_up_after_remove(parent_node, deleted, nodes_to_clean,tree);
    } else if (node->children.size() == 1) {
        //node with one child left
        auto &parent_children = parent_node->children;
        //detach from parent
        auto iter =
            std::find(parent_children.begin(), parent_children.end(), node);
        MAT::Node *child = node->children[0];
        //child may have its parental state changed
        if (merge_mutation_single_child(child, node->mutations)) {
            nodes_to_clean.push_back(child);
        }
        //fix boundary one alleles if child have only one child (parent of condensed node)
        if (child->children.size() <= 1) {
            auto &child_mut = child->mutations;
            for (auto &mut : child_mut) {
                mut.set_boundary_one_hot(0xf & (~mut.get_all_major_allele()));
            }
            child_mut.mutations.erase(
                std::remove_if(child_mut.begin(), child_mut.end(),
                               [](const MAT::Mutation &mut) {
                                   return mut.get_all_major_allele() ==
                                          mut.get_par_one_hot();
                               }),
                child_mut.end());
        }
        //reattach , this replace node with child in parent_children vector
        *iter = child;
        child->parent = parent_node;
        deleted.insert((size_t)node);
        /*if (node==sibling_node) {
            moved_mut.merge(node->mutations,
                        MAT::Mutations_Collection::KEEP_SELF);
        }*/
        tree.all_nodes.erase(node->identifier);
        //delete node;
        return parent_node;
    }
    return node;
}
MAT::Node *replace_with_internal_node(MAT::Node *to_replace,
                                             MAT::Tree &tree) {
    MAT::Node *new_node = new MAT::Node();
    new_node->identifier = std::to_string(++tree.curr_internal_node);
    new_node->parent = to_replace->parent;
    auto &to_replace_parent_children = to_replace->parent->children;

    auto iter = std::find(to_replace_parent_children.begin(),
                          to_replace_parent_children.end(), to_replace);
    tree.all_nodes.emplace(new_node->identifier, new_node);
    *iter = new_node;

    new_node->children.push_back(to_replace);
    to_replace->parent = new_node;
    return new_node;
}
//given par_nuc and all_major_alleles in mut, output one of the major allele, preferably the parental allele
static nuc_one_hot set_state(MAT::Mutation &mut) {
    nuc_one_hot state = mut.get_par_one_hot() & mut.get_all_major_allele();
    if (state) {
        mut.set_mut_one_hot(state);
        return state;
    }
    state = mut.get_all_major_allele().choose_first();
    mut.set_mut_one_hot(state);
    return state;
}

/**
 * @brief add mutation originating from either the dst node, or moved src node, not both to the shared node
    formed by splitting the branch
 * @param[out] shared_node_mutations_out shared node mutation vector to append 
 * @param[in] sibling_mut the non-parent allele from either sibling (dst) or new node(src)
 * @param is_sibling whether it is the new node(src) (on right) or the sibling node (on left)
 * @return whether the allele on sibling node cannot be shared with the new node
 */
static bool
new_internal_single_node(MAT::Mutations_Collection &shared_node_mutations_out,
                         const MAT::Mutation &sibling_mut, bool is_sibling) {
    MAT::Mutation temp = sibling_mut;
    nuc_one_hot par_allele = sibling_mut.get_par_one_hot();
    nuc_one_hot sibling_nuc = sibling_mut.get_all_major_allele();
    //set anxxilary state for temp
    get_new_mut_binary(temp, is_sibling ? sibling_nuc : par_allele,
                       is_sibling ? par_allele : sibling_nuc);
    assert(temp.get_par_one_hot() & temp.get_all_major_allele());
    //as one of the child at shared node is in parental state at this loci, so it is always possible to follow parental state
    temp.set_mut_one_hot(temp.get_par_one_hot());
    //append if necessary
    if (temp.get_par_one_hot() != temp.get_all_major_allele() ||
        temp.get_boundary1_one_hot()) {
        shared_node_mutations_out.push_back(temp);
    }
    return sibling_mut.is_valid();
}
//Filter away unnecessary mutation before adding
static bool add_mut(const MAT::Mutation &to_add, nuc_one_hot par_nuc,
                    MAT::Mutations_Collection &out) {
    if (to_add.get_all_major_allele() != par_nuc ||
        to_add.get_boundary1_one_hot()) {
        out.push_back(to_add);
        out.back().set_par_one_hot(par_nuc);
        return out.back().is_valid() &&(par_nuc&out.back().get_all_major_allele());
    }
    return false;
}

// sibling node left, new node right
char merge_new_node_mutations(
    const MAT::Mutations_Collection &new_node_mutations,
    const MAT::Mutations_Collection &sibling_node_mutations,
    MAT::Mutations_Collection &shared_node_mutations_out,
    MAT::Mutations_Collection &sibling_node_mutations_out,
    MAT::Mutations_Collection &new_node_mutations_out) {
    char flags = 0;
    auto sibling_iter = sibling_node_mutations.begin();
    auto sibling_end = sibling_node_mutations.end();
    for (const auto &mut : new_node_mutations) {
        while (sibling_iter != sibling_end &&
               sibling_iter->get_position() < mut.get_position()) {
            //unique to sibling node
            flags |= (new_internal_single_node(shared_node_mutations_out,
                                               *sibling_iter, true)
                      << SIBLING_UNIQUE_SHAMT);
            sibling_node_mutations_out.push_back(*sibling_iter);
            sibling_iter++;
        }
        if (sibling_iter != sibling_end &&
            sibling_iter->get_position() == mut.get_position()) {
            //sibling node and new node have mutation at the same loci
            shared_node_mutations_out.push_back(mut);
            auto &shared_node_output_mut = shared_node_mutations_out.back();
            nuc_one_hot sibling_mut = sibling_iter->get_all_major_allele();
            //set anxuliary states of mut
            bool shared =
                get_new_mut_binary(shared_node_output_mut, sibling_mut,
                                   mut.get_all_major_allele());
            flags |= (shared << HAVE_SHARED_SHAMT);
            flags |= ((!shared) << SIBLING_UNIQUE_SHAMT);
            //set par_nuc and mut_nuc
            nuc_one_hot shared_state = set_state(shared_node_output_mut);
            shared_node_output_mut.set_par_mut(sibling_iter->get_par_one_hot(),
                                               shared_state);
            //remove if this shared_node_mut is unnecessary
            if (shared_node_output_mut.get_par_one_hot() ==
                    shared_node_output_mut.get_all_major_allele() &&
                (!shared_node_output_mut.get_boundary1_one_hot())) {
                shared_node_mutations_out.mutations.pop_back();
            }
            //update parent state of both sibling node and new node
            flags |= (add_mut(*sibling_iter, shared_state,
                              sibling_node_mutations_out)
                      << SIBLING_INCONSISTENT_SHAMT);
            flags |= (add_mut(mut, shared_state, new_node_mutations_out)
                      << NEW_NODE_INCONSISTENT_SHAMT);
            //keep track of shared muts which make the mutation on the edge from shared node to sibling node unnecessary
            //They need to merge back to sibling node if splitting the edge between dst and its parent is unnecessary
            sibling_iter++;
        } else {
            //unique to new node, similar treatment as mutations unique to sibling node
            new_internal_single_node(shared_node_mutations_out, mut, false);
            if (mut.get_par_one_hot() != mut.get_all_major_allele() ||
                mut.get_boundary1_one_hot()) {
                new_node_mutations_out.push_back(mut);
            }
        }
    }
    while (sibling_iter != sibling_end) {
        //mutation unique to sibling node
        flags |= (new_internal_single_node(shared_node_mutations_out,
                                           *sibling_iter, true)
                  << SIBLING_UNIQUE_SHAMT);
        sibling_node_mutations_out.push_back(*sibling_iter);
        sibling_iter++;
    }
    // assert(have_shared);

    return flags;
}
//Set mutation vector of a node, with post processing boundary1_allele if that node only have one child
void update_src_mutation(MAT::Node *src,
                         MAT::Mutations_Collection &this_unique) {
    if (src->children.size() <= 1) {
        for (auto &mut : this_unique) {
            mut.set_boundary_one_hot(0xf & (~mut.get_all_major_allele()));
        }
        this_unique.remove_boundary_only();
    }
    // assert(!common.empty());
    src->mutations.swap(this_unique);
}

MAT::Node *add_as_sibling(MAT::Node *&src, MAT::Node *&dst,
                          MAT::Mutations_Collection &other_unique,
                          MAT::Mutations_Collection &common, MAT::Tree &tree,
                          char flag, std::vector<MAT::Node *> &nodes_to_clean) {
    //split branch
    MAT::Node *new_node = replace_with_internal_node(dst, tree);
    //set mutations
    new_node->mutations.swap(common);
    if (dst->is_leaf()) {
        other_unique.remove_boundary_only();
    }
    dst->mutations.swap(other_unique);
    new_node->children.push_back(src);
    src->parent = new_node;

    if (flag & (1 << SIBLING_INCONSISTENT_SHAMT)) {
        nodes_to_clean.push_back(dst);
    }
    return new_node;
}

static MAT::Node *place_node_LCA(MAT::Node *&src, MAT::Node *parent,
                                 MAT::Tree &tree,
                                 MAT::Mutations_Collection &mutations,
                                 MAT::Node *sibling,
                                 std::vector<MAT::Node *> &nodes_to_clean) {
    MAT::Mutations_Collection this_unique;
    MAT::Mutations_Collection other_unique;
    MAT::Mutations_Collection common;
    auto flags = merge_new_node_mutations(mutations, sibling->mutations, common,
                                          other_unique, this_unique);
        //otherwise, they have shared mutation, split the branch
        update_src_mutation(src, this_unique);
        MAT::Node *new_node = add_as_sibling(src, sibling, other_unique, common,
                                             tree, flags, nodes_to_clean);
        return new_node->parent;
}

static MAT::Node *place_node(MAT::Node *&src, MAT::Node *dst, MAT::Tree &tree,
                             MAT::Mutations_Collection &mutations,
                             std::vector<MAT::Node *> &nodes_to_clean) {
    MAT::Mutations_Collection this_unique;
    MAT::Mutations_Collection other_unique;
    MAT::Mutations_Collection common;
    auto flags = merge_new_node_mutations(mutations, dst->mutations, common,
                                          other_unique, this_unique);
    update_src_mutation(src, this_unique);
        //split branch
        MAT::Node *new_node = add_as_sibling(src, dst, other_unique, common,
                                             tree, flags, nodes_to_clean);
        return new_node->parent;
}

void move_node(MAT::Node *src, MAT::Node *dst,
               std::vector<MAT::Node *> &altered_node, MAT::Tree &tree,
               std::unordered_set<size_t> &deleted,
               std::vector<MAT::Node *> &nodes_to_clean
#ifdef CHECK_PRIMARY_MOVE
               //,
               //Original_State_t original_state
#endif
) {
    MAT::Mutations_Collection mutations;
    std::vector<MAT::Node *> dst_to_root_path;
    std::vector<MAT::Node *> src_to_root_path;
    //find the path from src to dst again, and make sure dst didn't got moved under src by some previous moves
    if(!find_path_no_dfs(dst_to_root_path, src_to_root_path, src, dst)){return;}
    //perserve ambiguous/boundary mutations of src
    for (const auto &mut : src->mutations) {
        if (mut.get_par_one_hot() != mut.get_all_major_allele() ||
            mut.get_boundary1_one_hot()) {
            mutations.push_back(mut);
        }
    }
    //merge mutation up to src_branch_node
    for (size_t idx = 1; idx < src_to_root_path.size(); idx++) {
        mutations.merge(src_to_root_path[idx]->mutations,
                        MAT::Mutations_Collection::KEEP_SELF);
    }

#ifdef CHECK_PRIMARY_MOVE
    MAT::Mutations_Collection before_move;
    get_mutation_set_from_root(src, before_move);
#endif

    if (!dst_to_root_path.empty()) {
        //not moving to LCA,merge mutation down to parent of dst
        for (int idx = dst_to_root_path.size() - 1; idx > 0; idx--) {
            mutations.merge(dst_to_root_path[idx]->mutations,
                            MAT::Mutations_Collection::INVERT_MERGE);
        }
    }
    if (src->is_leaf()) {
        mutations.remove_boundary_only();
    }
    MAT::Node *src_parent = src->parent;
    auto &src_parent_children = src_parent->children;
    //remove src
    auto iter =
        std::find(src_parent_children.begin(), src_parent_children.end(), src);
    src_parent_children.erase(iter);

    nodes_to_clean.push_back(src);
    src->changed=true;
    nodes_to_clean.push_back(dst);
    dst->changed=true;

    MAT::Node *dst_altered = dst;
    //actually placing the node
    if (dst_to_root_path.empty()) {
        dst_altered = place_node_LCA(src, dst, tree, mutations,
                                     src_to_root_path.back(), nodes_to_clean);
        src_to_root_path.back()->changed=true;
    } else {
        dst_altered = place_node(src, dst, tree, mutations, nodes_to_clean);
    }
    //push nodes with altered children for backward pass
    altered_node.push_back(
        clean_up_after_remove(src_parent, deleted, nodes_to_clean,tree));
    altered_node.back()->changed=true;
    altered_node.push_back(dst_altered);
    altered_node.back()->changed=true;
#ifdef CHECK_PRIMARY_MOVE
    MAT::Mutations_Collection after_move;
    get_mutation_set_from_root(src, after_move);
    auto size = before_move.size();
    for (size_t idx = 0; idx < size; idx++) {
        assert(before_move[idx].get_position() ==
               after_move[idx].get_position());
        assert(before_move[idx].get_mut_one_hot() ==
               after_move[idx].get_mut_one_hot());
    }
    //check_samples(tree.root, original_state, &tree);
#endif
}