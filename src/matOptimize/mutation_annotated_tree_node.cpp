#include "mutation_annotated_tree.hpp"
#include <cstddef>
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


Mutation_Annotated_Tree::Node::Node (size_t id) {
    node_id = id;
    have_masked=false;
    parent = NULL;
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node(const Node &other, Node *parent, Tree *tree,bool copy_mutations)
    :  branch_length(other.branch_length),
       node_id(other.node_id), parent(parent) {
    children.reserve(other.children.size());
    if (copy_mutations) {
        mutations=other.mutations;
    }
    have_masked=other.have_masked;
    for (auto c : other.children) {
        children.push_back(new Node(*c, this, tree,copy_mutations));
    }
    tree->register_node_serial(this);
}
Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier) {
    Node* n = create_node();
    node_names.emplace(n->node_id,identifier);
    node_name_to_idx_map.emplace(identifier,n->node_id);
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node () {
    auto new_node_id=node_idx++;
    Node* n = new Node(new_node_id);
    size_t num_annotations = get_num_annotations();
    n->clade_annotations.resize(num_annotations,"");
    register_node_serial(n);
    return n;
}

Node* Mutation_Annotated_Tree::Tree::get_node_c_str (const char* identifier) const {
    return get_node(std::string(identifier));
}

void Mutation_Annotated_Tree::Tree::rename_node(size_t old_nid, std::string new_nid) {
    if (new_nid=="") {
        auto iter=node_names.find(old_nid);
        if (iter!=node_names.end()) {
            node_name_to_idx_map.erase(iter->second);
            node_names.erase(iter);            
        }
        return;
    }
    auto ins_result=node_names.emplace(old_nid,new_nid);
    if (ins_result.second==false) {
        auto ori_name=ins_result.first->second;
        ins_result.first->second=new_nid;
        node_name_to_idx_map.erase(ori_name);
        node_name_to_idx_map.emplace(new_nid,old_nid);
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
void Mutation_Annotated_Tree::Node::populate_ignored_range() {
    ignore.clear();
    for (const auto& mut : mutations) {
        if (mut.get_all_major_allele()==0xf) {
            if (ignore.empty()||ignore.back().second!=(mut.get_position()-1)) {
                ignore.emplace_back(mut.get_position(),mut.get_position());
            } else {
                ignore.back().second=mut.get_position();
            }
        }
    }
    if (!ignore.empty()) {
        ignore.emplace_back(INT_MAX,INT_MAX);
    }
}
size_t Mutation_Annotated_Tree::Node::get_num_leaves() const{
    if (children.empty()) {
        return 1;
    }
    size_t leaf_count=0;
    for (auto child : children) {
        leaf_count+=child->get_num_leaves();
    }
    return leaf_count;
}