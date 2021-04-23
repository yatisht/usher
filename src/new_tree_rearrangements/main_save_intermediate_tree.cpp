#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
namespace MAT=Mutation_Annotated_Tree;
static MAT::Node* replicate_nodes(MAT::Node* node,MAT::Tree* tree){
    MAT::Node* new_node=new MAT::Node(node->identifier,node->branch_length);
    new_node->tree=tree;
    tree->all_nodes.emplace(node->identifier,new_node);
    for(MAT::Node* child:node->children){
        MAT::Node* new_child=replicate_nodes(child, tree);
        new_child->parent=new_node;
        new_node->children.push_back(new_child);
    }
    for(const MAT::Mutation& mut:node->mutations){
        if (mut.is_valid()) {
            new_node->mutations.push_back(mut);
        }
    }
    return new_node;
}
void save_intermediate_tree(MAT::Tree& tree_to_save,const std::string& path){
    MAT::Tree new_tree(tree_to_save);
    new_tree.all_nodes.clear();
    new_tree.root=replicate_nodes(new_tree.root, &new_tree);
    fix_condensed_nodes(&new_tree);
    Mutation_Annotated_Tree::save_mutation_annotated_tree(new_tree, path);
}