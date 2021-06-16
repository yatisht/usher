#include "mutation_annotated_tree.hpp"
#include <queue>
using namespace Mutation_Annotated_Tree;
bool check_grand_parent(const Mutation_Annotated_Tree::Node* node,const Mutation_Annotated_Tree::Node* grand_parent){
    const Mutation_Annotated_Tree::Node* cur=node;
    while (cur) {
        if(cur==grand_parent) return true;
        cur=cur->parent;
    }
    return false;
}
static void remove_child_helper(Node* child_to_remove,std::vector<Node*>& removed_nodes, bool delete_child){
        auto parent=child_to_remove->parent;
        auto iter=std::find(parent->children.begin(),parent->children.end(),child_to_remove);
        assert(iter!=parent->children.end());
        parent->children.erase(iter);
        if (parent->children.size()==0) {
            remove_child_helper(parent,removed_nodes,true);
        }
        if (delete_child) {
            removed_nodes.push_back(child_to_remove);
            delete child_to_remove;
        }
}

void Mutation_Annotated_Tree::remove_child(Mutation_Annotated_Tree::Node* child_to_remove,std::vector<Node*>& removed_nodes){
        remove_child_helper(child_to_remove, removed_nodes, false);
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
      identifier(other.identifier), parent(parent){
    children.reserve(other.children.size());
    if (copy_mutations) {
        mutations=other.mutations;
    }
    for (auto c : other.children) {
        children.push_back(new Node(*c, this, tree,copy_mutations));
    }
    tree->all_nodes.emplace(other.identifier,this);
}

// Get the last common ancestor of two node identifiers. Return NULL if does not
// exist
Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::LCA (const Mutation_Annotated_Tree::Tree& tree, const std::string& nid1, const std::string& nid2) {
    TIMEIT();
    
    if ((tree.get_node(nid1) == NULL) || (tree.get_node(nid2) == NULL)) {
        return NULL;
    }

    auto n2_ancestors = tree.rsearch(nid2, true);

    for (auto anc1: tree.rsearch(nid1, true)) {
        for (auto anc2: n2_ancestors) {
            if (anc1 == anc2) {
                return anc1;
            }
        }
    }

    return NULL;
}
void Mutation_Annotated_Tree::Tree::remove_node_helper (std::string nid, bool move_level) { 
    auto it = all_nodes.find(nid);
    if (it == all_nodes.end()) {
        fprintf(stderr, "ERROR: Tried to remove node identifier %s but it was not found!\n", nid.c_str());
        exit(1);
    }
    Node* source = it->second;
    Node* curr_parent = source->parent;
    
    if (curr_parent != NULL) {
        // Remove source from curr_parent
        auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), source);
        assert (iter != curr_parent->children.end());
        curr_parent->children.erase(iter);

        // Remove parent if it no longer has any children
        if (curr_parent->children.size() == 0) {
            if (curr_parent == root) {
                fprintf(stderr, "ERROR: Tree empty!\n");
                exit(1);
            }
            remove_node_helper (curr_parent->identifier, move_level);
        }
        // Move the remaining child one level up if it is the only child of its parent 
        else if (move_level && (curr_parent->children.size() == 1)) {
            auto child = curr_parent->children[0];
            if (curr_parent->parent != NULL) {
                for (size_t k=0; k < curr_parent->clade_annotations.size(); k++) {
                    if (child->clade_annotations[k] == "") {
                        child->clade_annotations[k] = curr_parent->clade_annotations[k];
                    }
                }
                child->parent = curr_parent->parent;
                child->branch_length += curr_parent->branch_length;

                child->mutations.merge(curr_parent->mutations, 1);

                curr_parent->parent->children.push_back(child);
                
                iter = std::find(curr_parent->parent->children.begin(), curr_parent->parent->children.end(), curr_parent);
                assert(iter != curr_parent->parent->children.end());
                curr_parent->parent->children.erase(iter);
                
                // Update levels of source descendants
                std::queue<Node*> remaining_nodes;
                remaining_nodes.push(child);
                while (remaining_nodes.size() > 0) {
                    Node* curr_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    for (auto c: curr_node->children) {
                        remaining_nodes.push(c);
                    }
                }
            }

            auto par_it = all_nodes.find(curr_parent->identifier);
            assert (par_it != all_nodes.end());
            all_nodes.erase(par_it);
            delete curr_parent;
        }
    }

    //Remove source and descendants from all_nodes
    std::queue<Node*> desc;
    desc.push(source);
    while (desc.size() > 0) {
        Node* curr_node = desc.front();
        desc.pop();
        for (auto c: curr_node->children) {
            desc.push(c);
        }
        it = all_nodes.find(curr_node->identifier);
        all_nodes.erase(it);
        delete curr_node;
    }
}

void Mutation_Annotated_Tree::Tree::remove_node (std::string nid, bool move_level) { 
    TIMEIT();
    remove_node_helper (nid, move_level);
}

void Mutation_Annotated_Tree::Tree::move_node (std::string source_id, std::string dest_id) {
    Node* source = all_nodes[source_id];
    Node* destination = all_nodes[dest_id];
    Node* curr_parent = source->parent;

    source->parent = destination;
    source->branch_length = -1.0; // Invalidate source branch length

    destination->children.push_back(source);

    // Remove source from curr_parent
    auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), source);
    curr_parent->children.erase(iter);
    if (curr_parent->children.size() == 0) {
        remove_node(curr_parent->identifier, true);
    }
    
    // Update levels of source descendants
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(source);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
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

bool Mutation_Annotated_Tree::Tree::is_ancestor (std::string anc_id, std::string nid) const {
    Node* node = get_node(nid);
    while (node->parent != NULL) {
        node = node->parent;
        if (node->identifier == anc_id) {
            return true;
        }
    }
    return false; 
}
Node* Mutation_Annotated_Tree::Tree::get_node_c_str (char* identifier) const{
                return get_node(std::string(identifier));
            }

void Mutation_Annotated_Tree::Tree::rename_node(std::string old_nid, std::string new_nid) {
    auto n = get_node(old_nid);
    if (n != NULL) {
        n->identifier = new_nid;
        all_nodes.erase(old_nid);
        all_nodes[new_nid] = n;
    }
    else {
        fprintf(stderr, "ERROR: %s not found in the Tree!\n", old_nid.c_str());
        exit(1);
    }
}
nuc_one_hot get_parent_state(Node* ancestor,int position){
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