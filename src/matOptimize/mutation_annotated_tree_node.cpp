#include "mutation_annotated_tree.hpp"
#include <queue>
using namespace Mutation_Annotated_Tree;
bool check_grand_parent(const Mutation_Annotated_Tree::Node* node,const Mutation_Annotated_Tree::Node* grand_parent) {
    const Mutation_Annotated_Tree::Node* cur=node;
    while (cur) {
        if(cur==grand_parent) return true;
        cur=cur->parent;
    }
    return false;
}

/* === Node === */
bool Mutation_Annotated_Tree::Node::is_leaf () const {
    return (children.size() == 0);
}

bool Mutation_Annotated_Tree::Node::is_root() {
    return (parent == NULL);
}

Mutation_Annotated_Tree::Node::Node() {
    identifier = "";
    parent = NULL;
    branch_length = -1.0;
    clade_annotations.clear();
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, float len) {
    identifier = id;
    parent = NULL;
    branch_length = len;
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, Node* p, float len) {
    identifier = id;
    parent = p;
    branch_length = len;
    mutations.clear();
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Node::add_child(Node *new_child,Mutation_Annotated_Tree::Tree* tree) {
    Mutation_Annotated_Tree::Node* ret=nullptr;
    if (is_leaf()) {
        std::string old_name(identifier);
        tree->rename_node(identifier, std::to_string(++tree->curr_internal_node));
        tree->create_node(old_name,this);
        ret= this;
    }
    new_child->parent = this;
    children.push_back(new_child);
    return ret;
}
Mutation_Annotated_Tree::Node::Node(const Node &other, Node *parent, Tree *tree,bool copy_mutations)
    :  branch_length(other.branch_length),
       identifier(other.identifier), parent(parent) {
    children.reserve(other.children.size());
    if (copy_mutations) {
        mutations=other.mutations;
    }
    for (auto c : other.children) {
        children.push_back(new Node(*c, this, tree,copy_mutations));
    }
    tree->all_nodes.emplace(other.identifier,this);
}
Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, float branch_len, size_t num_annotations) {
    all_nodes.clear();
    Node* n = new Node(identifier, branch_len);
    for (size_t k=0; k < num_annotations; k++) {
        n->clade_annotations.emplace_back("");
    }
    root = n;
    all_nodes[identifier] = root;
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, Node* par, float branch_len) {
    if (all_nodes.find(identifier) != all_nodes.end()) {
        fprintf(stderr, "Error: %s already in the tree!\n", identifier.c_str());
        exit(1);
    }
    Node* n = new Node(identifier, par, branch_len);
    size_t num_annotations = get_num_annotations();
    for (size_t k=0; k < num_annotations; k++) {
        n->clade_annotations.emplace_back("");
    }
    all_nodes[identifier] = n;
    par->children.push_back(n);
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, std::string const& parent_id, float branch_len) {
    Node* par = all_nodes[parent_id];
    return create_node(identifier, par, branch_len);
}

Node* Mutation_Annotated_Tree::Tree::get_node_c_str (char* identifier) const {
    return get_node(std::string(identifier));
}

void Mutation_Annotated_Tree::Tree::rename_node(std::string old_nid, std::string new_nid) {
    auto n = get_node(old_nid);
    if (n != NULL) {
        n->identifier = new_nid;
        all_nodes.erase(old_nid);
        all_nodes[new_nid] = n;
    } else {
        fprintf(stderr, "ERROR: %s not found in the Tree!\n", old_nid.c_str());
        exit(1);
    }
}
nuc_one_hot get_parent_state(Node* ancestor,int position) {
    auto iter = ancestor->mutations.find(position);
    if (iter == ancestor->mutations.end()) {
        ancestor = ancestor->parent;
    } else {
        return iter->get_par_one_hot();
    }
    while (ancestor) {
        auto iter = ancestor->mutations.find(position);
        if (iter == ancestor->mutations.end()) {
            ancestor = ancestor->parent;
        } else {
            return iter->get_mut_one_hot();
        }
    }
    return Mutation::refs[position];
}