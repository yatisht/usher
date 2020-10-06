#include "tree.hpp"
#include <iomanip>

bool Node::is_leaf () {
    return (children.size() == 0);
}

bool Node::is_root() {
    return (parent == NULL);
}

Node::Node() {
    level = 0;
    identifier = "";
    parent = NULL;
    branch_length = -1.0;
}

Node::Node (std::string id, float len) {
    identifier = id;
    parent = NULL;
    level = 1;
    branch_length = len;
}

Node::Node (std::string id, Node* p, float len) {
    identifier = id;
    parent = p;
    level = p->level + 1;
    branch_length = len;
}

size_t Tree::get_max_level () {
    return max_level;
}
        
void Tree::rename_node(std::string old_nid, std::string new_nid) {
    Node* n = get_node(old_nid);
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

std::vector<Node*> Tree::get_leaves() {
    std::vector<Node*> leaves;
    for (auto x: all_nodes) {
        auto node = x.second;
        if (node->is_leaf()) {
            leaves.push_back(node);
        }
    }
    return leaves;
}

size_t Tree::get_num_leaves(Node* node) {
    if (node->is_leaf()) {
        return 1;
    }
    size_t num_leaves = 0;
    for (auto c: node->children) {
        num_leaves += get_num_leaves(c);
    }
    return num_leaves;
}

size_t Tree::get_num_leaves() {
    return get_num_leaves(root);
}

std::vector<Node*> Tree::get_leaves(std::string nid) {
    std::vector<Node*> leaves;
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves;
}

void Tree::create_node (std::string identifier, float branch_len) {
    all_nodes.clear();
    max_level = 1;
    Node* n = new Node(identifier, branch_len);
    root = n;
    all_nodes[identifier] = root;
}

void Tree::create_node (std::string identifier, std::string parent_id, float branch_len) {
    Node* par = all_nodes[parent_id];
    Node* n = new Node(identifier, par, branch_len);
    if (all_nodes.find(identifier) != all_nodes.end()) {
        fprintf(stderr, "Error: %s already in the tree!\n", identifier.c_str());
        exit(1);
    }
    all_nodes[identifier] = n;
    par->children.push_back(n);
    if (n->level > max_level) {
        max_level = n->level;
    }
}

Node* Tree::get_node (std::string nid) {
    if (all_nodes.find(nid) != all_nodes.end()) {
        return all_nodes[nid];
    }
    return NULL;

}

bool Tree::is_ancestor (std::string anc_id, std::string nid) {
    Node* node = all_nodes[nid];
    while (node->parent != NULL) {
        node = node->parent;
        if (node->identifier == anc_id) {
            return true;
        }
    }
    return false; 
}

std::vector<Node*> Tree::rsearch (std::string nid) {
    Node* node = all_nodes[nid];
    std::vector<Node*> ancestors;
    while (node->parent != NULL) {
        ancestors.push_back(node->parent);
        node = node->parent;
    }
    return ancestors;
}

void Tree::remove_node_helper (std::string nid, bool move_level) { 
    auto it = all_nodes.find(nid);
    Node* source = it->second;
    Node* curr_parent = source->parent;
    
    if (curr_parent != NULL) {
        // Remove source from curr_parent
        auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), source);
        assert (iter != curr_parent->children.end());
        curr_parent->children.erase(iter);

        // Remove parent if it no longer has any children
        if (curr_parent->children.size() == 0) {
            remove_node_helper (curr_parent->identifier, move_level);
        }
        // Move the remaining child one level up if it is the only child of its parent 
        else if (move_level && (curr_parent->children.size() == 1)) {
            auto child = curr_parent->children[0];
            if (curr_parent->parent != NULL) {
                child->parent = curr_parent->parent;
                child->level = curr_parent->parent->level + 1;
                child->branch_length += curr_parent->branch_length;

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
                    curr_node->level = curr_node->parent->level + 1;
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

void Tree::remove_node (std::string nid, bool move_level) { 
    remove_node_helper (nid, move_level);

    // Update max level
    size_t new_max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > new_max_level) {
            new_max_level = x.second->level;
        }
    }
    max_level = new_max_level;
}

void Tree::move_node (std::string source_id, std::string dest_id) {
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
        curr_node->level = curr_node->parent->level + 1;
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    
    // Update max level
    size_t new_max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > new_max_level) {
            new_max_level = x.second->level;
        }
    }
    max_level = new_max_level;
}

std::vector<Node*> Tree::breadth_first_expansion() {
    return breadth_first_expansion(root->identifier);
}

std::vector<Node*> Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        traversal.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }

    return traversal;
}

void Tree::depth_first_expansion_helper(Node* node, std::vector<Node*>& vec) {
    vec.push_back(node);
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec);
    }
}

std::vector<Node*> Tree::depth_first_expansion() {
    std::vector<Node*> traversal;
    if (root != NULL) {
        depth_first_expansion_helper(root, traversal);
    }
    return traversal;
}

std::vector<Node*> Tree::depth_first_expansion(Node* node) {
    std::vector<Node*> traversal;
    depth_first_expansion_helper(node, traversal);
    return traversal;
}

std::string TreeLib::get_newick_string (Tree& T, Node* node, bool print_internal, bool print_branch_len) {
    std::string newick_string = "";

    std::vector<Node*> traversal = T.depth_first_expansion(node);
    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->branch_length;
        if (curr_level < level) {
            if (!prev_open) {
                newick_string += ",";
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                newick_string += "(";
                prev_open = true;
            }
            if (n->is_leaf()) {
                newick_string += n->identifier;
                if ((print_branch_len) && (branch_length >= 0)) {
                    std::stringstream ss;
                    ss << branch_length;
                    newick_string += ":" + ss.str();
                }
                prev_open = false;
            }
            else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        }
        else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                newick_string += ")";
                if (print_internal){
                    newick_string += node_stack.top();
                }
                if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
                    std::stringstream ss;
                    ss << branch_length_stack.top();
                    newick_string += ":" + ss.str();
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->is_leaf()) {
                newick_string += "," + n->identifier;
                if ((print_branch_len) && (branch_length >= 0)) {
                    std::stringstream ss;
                    ss << branch_length;
                    newick_string += ":" + ss.str(); 
                }
            }
            else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        }
        else {
            prev_open = false;
            if (n->is_leaf()) {
                newick_string += "," + n->identifier;
                if ((print_branch_len) && (branch_length >= 0)) {
                    std::stringstream ss;
                    ss << branch_length;
                    newick_string += ":" + ss.str();
                }
            }
            else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        }
        curr_level = level;
    }
    size_t remaining = node_stack.size();
    for (size_t i = 0; i < remaining; i++) {
        newick_string += ")";
        if (print_internal) {
            newick_string += node_stack.top();
        }
        if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
            std::stringstream ss;
            ss << branch_length_stack.top();
            newick_string += ":" + ss.str(); 
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    newick_string += ";";
    return newick_string;
}

std::string TreeLib::get_newick_string (Tree& T, bool print_internal, bool print_branch_len) {
    return get_newick_string(T, T.root, print_internal, print_branch_len);
}

void TreeLib::split (std::string s, char delim, std::vector<std::string>& words) {
    std::string curr = "";
    for (auto c: s) {
        if (c == delim) {
            words.push_back(std::move(curr));
            curr = "";
        }
        else {
            curr += c;
        }
    }
    if (curr != "") {
        words.push_back(std::move(curr));
    }
}

void TreeLib::split (std::string s, std::vector<std::string>& words) {
    std::string curr = "";
    std::vector<std::string> ret;
    
    // Used to split string around spaces.
    std::istringstream ss(s);

    std::string word;
    // Traverse through all words
    while (ss >> word) {
        words.push_back(std::move(word));
    };
}

Tree TreeLib::create_tree_from_newick_string (std::string newick_string) {
    Tree T;

    std::vector<std::string> leaves;
    std::vector<size_t> num_open;
    std::vector<size_t> num_close;
    std::stack<float> branch_len;

    std::vector<std::string> s1;
    split(newick_string, ',', s1);

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        bool stop = false;
        bool branch_start = false;
        std::string leaf = "";
        std::string branch = "";
        for (auto c: s) {
            if (c == ':') {
                stop = true;
                branch = "";
                branch_start = true;
            }
            else if (c == '(') {
                no++;
            }
            else if (c == ')') {
                stop = true;
                nc++;
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                branch_len.push(len);
                branch_start = false;
            }
            else if (!stop) {
                leaf += c;
                branch_start = false;
            }
            else if (branch_start) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        leaves.push_back(std::move(leaf));
        num_open.push_back(no);
        num_close.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        branch_len.push(len);
    }

    if (num_open.size() != num_close.size()) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    T.curr_internal_node = 0;
    std::stack<std::string> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = std::to_string(++T.curr_internal_node);
            if (parent_stack.size() == 0) {
                T.create_node(nid, branch_len.top());
                branch_len.pop();
            }
            else {
                T.create_node(nid, parent_stack.top(), branch_len.top());
                branch_len.pop();
            }
            parent_stack.push(nid);
        }
        T.create_node(leaf, parent_stack.top(), branch_len.top());
        branch_len.pop();
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
        }
    }

    return T;
}

Tree TreeLib::create_tree_from_newick (std::string filename) {
    std::ifstream infile(filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the tree file: %s!\n", filename.c_str());
        exit(1);
    }
    std::string newick_string;
    std::getline(infile, newick_string);

    return create_tree_from_newick_string(newick_string);
}
