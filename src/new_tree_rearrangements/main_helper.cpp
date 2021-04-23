#include "tree_rearrangement_internal.hpp"
#include "check_samples.hpp"
#include "Fitch_Sankoff.hpp"
#include "src/new_tree_rearrangements/Twice_Bloom_Filter.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
namespace MAT=Mutation_Annotated_Tree;
void add_root(MAT::Tree *tree) {
    MAT::Node *old_root = tree->root;
    MAT::Node *new_root = new MAT::Node();
    new_root->tree=tree;
    new_root->identifier = std::to_string(++tree->curr_internal_node);
    std::string &node_name = new_root->identifier;
    new_root->children.push_back(old_root);
    tree->all_nodes.emplace(node_name, new_root);
    old_root->parent=new_root;
    tree->root=new_root;
}
void fix_condensed_nodes(MAT::Tree *tree) {
    std::vector<MAT::Node *> nodes_to_fix;
    for (auto iter : tree->all_nodes) {
        if (tree->condensed_nodes.count(iter.first) &&
            (!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for (auto node : nodes_to_fix) {
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier,
                          std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier, node);
    }
}


static void find_nodes_with_recurrent_mutations(
    const std::vector<MAT::Node *> &all_nodes,
    std::vector<MAT::Node *> &output) {
    Twice_Bloom_Filter filter;
    for (MAT::Node *n : all_nodes) {
        for (const MAT::Mutation &m : n->mutations) {
            filter.insert(m.get_position());
        }
    }
    for (MAT::Node *n : all_nodes) {
        for (const MAT::Mutation &m : n->mutations) {
            if (filter.query(m.get_position())) {
                output.push_back(n);
                break;
            }
        }
    }
}
 void
find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                   tbb::concurrent_vector<MAT::Node*> &output) {
    std::vector<MAT::Node *> nodes_with_recurrent_mutations;
    find_nodes_with_recurrent_mutations(bfs_ordered_nodes,
                                        nodes_with_recurrent_mutations);

    std::unordered_set<size_t> pushed_nodes;
    for (MAT::Node *this_node : nodes_with_recurrent_mutations) {
        while (this_node->parent) {
            auto result = pushed_nodes.insert(this_node->bfs_index);
            if (result.second) {
                output.push_back(this_node);
                this_node = this_node->parent;
            } else {
                break;
            }
        }
    }
}
