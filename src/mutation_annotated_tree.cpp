#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <string>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8
using Mutation_Annotated_Tree::Node;

int8_t Mutation_Annotated_Tree::get_nuc_id(char nuc) {
    int8_t ret = 0b1111;
    switch(nuc) {
        case 'a':
        case 'A': ret = 0b1;
                  break;
        case 'c':
        case 'C': ret = 0b10;
                  break;
        case 'g':
        case 'G': ret = 0b100; 
                  break;
        case 't':
        case 'T': ret = 0b1000; 
                  break;
        case 'R': ret = 0b101;
                  break;
        case 'Y': ret = 0b1010;
                  break;
        case 'S': ret = 0b110;
                  break;
        case 'W': ret = 0b1001;
                  break;
        case 'K': ret = 0b1100;
                  break;
        case 'M': ret = 0b11;
                  break;
        case 'B': ret = 0b1110;
                  break;
        case 'D': ret = 0b1101;
                  break;
        case 'H': ret = 0b1011;
                  break;
        case 'V': ret = 0b111;
        case 'n':
        case 'N': 
        default: ret = 0b1111;
                 break;
    }
    return ret;
}

// Sets bits at positions specified by nuc_vec to 1 in int8
int8_t Mutation_Annotated_Tree::get_nuc_id (std::vector<int8_t> nuc_vec) {
    int8_t ret = 0;
    int8_t one = 1;
    for (auto nuc: nuc_vec) {
        assert((nuc >= 0) && (nuc <=3));
        ret += (one << nuc);
    }
    return ret;
}

// Convert nuc_id back to IUPAC base 
char Mutation_Annotated_Tree::get_nuc (int8_t nuc_id) {
    char ret = 'N';
    //assert ((nuc_id >= 1) && (nuc_id <= 15));
    switch(nuc_id) {
        case 1: ret = 'A';
                break;
        case 2: ret = 'C';
                break;
        case 3: ret = 'M';
                break;
        case 4: ret = 'G';
                break;
        case 5: ret = 'R';
                break;
        case 6: ret = 'S';
                break;
        case 7: ret = 'V';
                break;
        case 8: ret = 'T';
                break;
        case 9: ret = 'W';
                break;
        case 10: ret = 'Y';
                break;
        case 11: ret = 'H';
                break;
        case 12: ret = 'K';
                break;
        case 13: ret = 'D';
                break;
        case 14: ret = 'B';
                break;
        default: ret = 'N';
                 break;
    }
    return ret;
}

// A:0, C:1, G:2, T:3 
int8_t Mutation_Annotated_Tree::get_nt (int8_t nuc_id) {
    int8_t ret = 0;
    switch(nuc_id) {
        case 1: ret = 0;
                break;
        case 2: ret = 1;
                break;
        case 4: ret = 2;
                break;
        case 8: ret = 3;
                break;
        default: ret = -1;
                 break;
    }
    return ret;
}

std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec (char c) {
    switch (c) {
        case 'a':
        case 'A': return std::vector<int8_t>{0};
        case 'c':
        case 'C': return std::vector<int8_t>{1};
        case 'g':
        case 'G': return std::vector<int8_t>{2};
        case 't':
        case 'T': return std::vector<int8_t>{3};
        case 'R': return std::vector<int8_t>{0,2};
        case 'Y': return std::vector<int8_t>{1,3};
        case 'S': return std::vector<int8_t>{1,2};
        case 'W': return std::vector<int8_t>{0,3};
        case 'K': return std::vector<int8_t>{2,3};
        case 'M': return std::vector<int8_t>{0,1};
        case 'B': return std::vector<int8_t>{1,2,3};
        case 'D': return std::vector<int8_t>{0,2,3};
        case 'H': return std::vector<int8_t>{0,1,3};
        case 'V': return std::vector<int8_t>{0,1,2};
        case 'n':
        case 'N': return std::vector<int8_t>{0,1,2,3};
        default: return std::vector<int8_t>{0,1,2,3};
    }
}

std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec_from_id (int8_t nuc_id) {
    return get_nuc_vec(get_nuc(nuc_id));
}

// Get newick stringstream for the input subtree rooted at some node (node) in 
// the input tree T. Boolean arguments decide whether
// internal node ids and branch lengths are printed. If last boolean argument is
// set, branch lengths from input tree are retained, otherwise, branch length
// for a branch is equal to the number of mutations annotated on that branch 
void Mutation_Annotated_Tree::write_newick_string (std::stringstream& ss, const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    TIMEIT();

    std::vector<Node*> traversal = T.depth_first_expansion(node);
    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->branch_length;
        if (!retain_original_branch_len) {
            branch_length = static_cast<float>(n->mutations.size());
        }
        if (curr_level < level) {
            if (!prev_open) {
                ss << ',';
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                ss << '(';
                prev_open = true;
            }
            if (n->is_leaf()) {
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                }
                else {
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
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
                ss << ')';
                if (print_internal){
                    ss << node_stack.top();
                }
                if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
                    ss << ':';
                    ss << branch_length_stack.top();
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->is_leaf()) {
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    ss << ',';
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                }
                else {
                    ss << ',';
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
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
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    ss << ',';
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                }
                else {
                    ss << ',';
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
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
        ss << ')';
        if (print_internal) {
            ss << node_stack.top();
        }
        if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
            ss << ':';
            ss << branch_length_stack.top();
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    ss << ';';
}

std::string Mutation_Annotated_Tree::get_newick_string (const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    std::stringstream newick_ss;
    write_newick_string(newick_ss, T, node, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
    return newick_ss.str();
}

std::string Mutation_Annotated_Tree::get_newick_string (const Tree& T, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    return get_newick_string(T, T.root, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
}

// Split string into words for a specific delimiter delim
void Mutation_Annotated_Tree::string_split (std::string const& s, char delim, std::vector<std::string>& words) {
    TIMEIT();
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if ((end_pos == start_pos) || end_pos >= s.length()) {
            break;
        }
        words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
    
}

// Split string into words (delimited by space, tabs etc.)
void Mutation_Annotated_Tree::string_split (std::string s, std::vector<std::string>& words) {
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

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::create_tree_from_newick_string (std::string newick_string) {
    TIMEIT();
    Tree T;

    std::vector<std::string> leaves;
    std::vector<size_t> num_open;
    std::vector<size_t> num_close;
    std::stack<float> branch_len;

    std::vector<std::string> s1;
    string_split(newick_string, ',', s1);

    num_open.reserve(s1.size());
    num_close.reserve(s1.size());
    
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
    std::stack<Node*> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = std::to_string(++T.curr_internal_node);
            Node* new_node = NULL;
            if (parent_stack.size() == 0) {
                new_node = T.create_node(nid, branch_len.top());
                branch_len.pop();
            }
            else {
                new_node = T.create_node(nid, parent_stack.top(), branch_len.top());
                branch_len.pop();
            }
            parent_stack.push(new_node);
        }
        T.create_node(leaf, parent_stack.top(), branch_len.top());
        branch_len.pop();
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
        }
    }

    if (T.root == NULL) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    return T;
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::create_tree_from_newick (std::string filename) {
    std::ifstream infile(filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the tree file: %s!\n", filename.c_str());
        exit(1);
    }
    std::string newick_string;
    std::getline(infile, newick_string);

    return create_tree_from_newick_string(newick_string);
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::load_mutation_annotated_tree (std::string filename) {
    TIMEIT();
    Tree tree;

    Parsimony::data data;

    std::ifstream inpfile(filename, std::ios::in | std::ios::binary);
    if (!inpfile) {
        fprintf(stderr, "ERROR: Could not load the mutation-annotated tree object from file: %s!\n", filename.c_str());
        exit(1);
    }
    data.ParseFromIstream(&inpfile);
    inpfile.close();
    //check if the pb has a metadata field
    bool hasmeta = (data.metadata_size()>0);
    if (!hasmeta) {
        fprintf(stderr, "WARNING: This pb does not include any metadata. Filling in default values\n");
    }
    tree = create_tree_from_newick_string(data.newick());
    auto dfs = tree.depth_first_expansion();
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, dfs.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t idx = r.begin(); idx < r.end(); idx++) {
               auto node = dfs[idx];
               auto mutation_list = data.node_mutations(idx);
               if (hasmeta) {
                   node->clade = data.metadata(idx).clade(); 
               } else {
                   node->clade = ""; 
               }
               for (int k = 0; k < mutation_list.mutation_size(); k++) {
                  auto mut = mutation_list.mutation(k);
                  Mutation m;
                  m.chrom = mut.chromosome();
                  m.position = mut.position();
                  if (!m.is_masked()) {
                     m.ref_nuc = (1 << mut.ref_nuc());
                     m.par_nuc = (1 << mut.par_nuc());
                     m.is_missing = false;
                     std::vector<int8_t> nuc_vec(mut.mut_nuc_size());
                     for (int n = 0; n < mut.mut_nuc_size(); n++) {
                        nuc_vec[n] = mut.mut_nuc(n);
                     }
                     m.mut_nuc = get_nuc_id(nuc_vec);
                  }
                  else {
                      // Mutation masked
                      m.ref_nuc = 0;
                      m.par_nuc = 0;
                      m.mut_nuc = 0;
                  }
                  node->add_mutation(m);
               }
            }
        }, ap);

    size_t num_condensed_nodes = static_cast<size_t>(data.condensed_nodes_size());
    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_condensed_nodes),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t idx = r.begin(); idx < r.end(); idx++) {
               auto cn = data.condensed_nodes(idx);
               tree.condensed_nodes.emplace(std::pair<std::string, std::vector<std::string>>(cn.node_name(), std::vector<std::string>(cn.condensed_leaves_size())));
               for (int k = 0; k < cn.condensed_leaves_size(); k++) {
                  tree.condensed_nodes[cn.node_name()][k] = cn.condensed_leaves(k);
                  tree.condensed_leaves.emplace(cn.condensed_leaves(k));
               }
            }
    }, ap);

    return tree;
}

void Mutation_Annotated_Tree::save_mutation_annotated_tree (Mutation_Annotated_Tree::Tree tree, std::string filename) {
    TIMEIT();
    Parsimony::data data;
    data.set_newick(get_newick_string(tree, false, true));

    auto dfs = tree.depth_first_expansion();

    for (size_t idx = 0; idx < dfs.size(); idx++) {
        auto meta = data.add_metadata();
        meta->set_clade(dfs[idx]->clade);
        auto mutation_list = data.add_node_mutations();
        for (auto m: dfs[idx]->mutations) {
            auto mut = mutation_list->add_mutation();
            mut->set_chromosome(m.chrom);
            mut->set_position(m.position);
            
            if (m.is_masked()) {
                mut->set_ref_nuc(-1);
                mut->set_par_nuc(-1);
            }
            else {
                int8_t j = get_nt(m.ref_nuc);
                assert (j >= 0);
                mut->set_ref_nuc(j);

                j = get_nt(m.par_nuc);
                assert(j >= 0);
                mut->set_par_nuc(j);

                mut->clear_mut_nuc();
                for (auto nuc: get_nuc_vec_from_id(m.mut_nuc)) {
                    mut->add_mut_nuc(nuc);
                }
            }
        }
    }

    // Add condensed nodes
    for (auto cn: tree.condensed_nodes) {
        auto cn_ptr = data.add_condensed_nodes();
        cn_ptr->set_node_name(cn.first);
        for (auto lid: cn.second) {
            cn_ptr->add_condensed_leaves(lid);
        }
    }

    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    data.SerializeToOstream(&outfile);
    outfile.close();
}

/* === Node === */
bool Mutation_Annotated_Tree::Node::is_leaf () const {
    return (children.size() == 0);
}

bool Mutation_Annotated_Tree::Node::is_root() {
    return (parent == NULL);
}

Mutation_Annotated_Tree::Node::Node() {
    level = 0;
    identifier = "";
    clade = "";
    parent = NULL;
    branch_length = -1.0;
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, float len) {
    identifier = id;
    parent = NULL;
    level = 1;
    branch_length = len;
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, Node* p, float len) {
    identifier = id;
    parent = p;
    level = p->level + 1;
    branch_length = len;
    mutations.clear();
}
#ifdef matToVCF

// Assumes mutations are added in chronological order. If a new mutation occurs
// at the same position, it should either be updated to the new allele or
// removed entirely (in case of reversal mutation)
void Mutation_Annotated_Tree::Node::add_mutation (Mutation mut) {
    auto iter = std::lower_bound(mutations.begin(), mutations.end(), mut);
    // check if mutation at the same position has occured before
    if ((iter != mutations.end()) && (iter->position == mut.position)) {
        // update to new allele
        if (iter->par_nuc != mut.mut_nuc) {
            iter->mut_nuc = mut.mut_nuc;
        }
        //reversal mutation
        else {
            std::vector<Mutation> tmp;
            for (auto m: mutations) {
                if (m.position != iter->position) {
                    tmp.emplace_back(m.copy());
                }
            }
            mutations.clear();
            for (auto m: tmp) {
                mutations.emplace_back(m.copy());
            }
        }
    }
    // new mutation
    else {
        mutations.insert(iter, mut);
    }
}

void Mutation_Annotated_Tree::Node::clear_mutations() {
    mutations.clear();
}
#endif
/* === Tree === */
size_t Mutation_Annotated_Tree::Tree::get_max_level () {
    size_t max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > max_level) {
            max_level = x.second->level;
        }
    }
    return max_level;
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

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::get_leaves(std::string nid) {
    std::vector<Node*> leaves;
    if (nid == "") {
        if (root == NULL) {
            return leaves;
        }
        nid = root->identifier;
    }
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

std::vector<std::string> Mutation_Annotated_Tree::Tree::get_leaves_ids(std::string nid) {
    std::vector<std::string> leaves_ids;
    if (nid == "") {
        if (root == NULL) {
            return leaves_ids;
        }
        nid = root->identifier;
    }
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves_ids.push_back(curr_node->identifier);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves_ids;
}

size_t Mutation_Annotated_Tree::Tree::get_num_leaves(Node* node) {
    if (node == NULL) {
        node = root;
    }

    if (node->is_leaf()) {
        return 1;
    }
    size_t num_leaves = 0;
    for (auto c: node->children) {
        num_leaves += get_num_leaves(c);
    }
    return num_leaves;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, float branch_len) {
    all_nodes.clear();
    Node* n = new Node(identifier, branch_len,this);
    root = n;
    all_nodes[identifier] = root;
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, Node* par, float branch_len) {
    Node* n = new Node(identifier, par, branch_len,this);
    if (all_nodes.find(identifier) != all_nodes.end()) {
        fprintf(stderr, "Error: %s already in the tree!\n", identifier.c_str());
        exit(1);
    }
    all_nodes[identifier] = n;
    par->children.push_back(n);
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, std::string const& parent_id, float branch_len) {
    Node* par = all_nodes[parent_id];
    return create_node(identifier, par, branch_len);
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::get_node (std::string nid) const {
    auto iter=all_nodes.find(nid);
    if (iter != all_nodes.end()) {
        if(iter->second->identifier!=nid){
            fprintf(stderr,"looking for %s but get %s \n",nid.c_str(),iter->second->identifier.c_str());
        }
        return iter->second;
    }
    return NULL;

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

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::rsearch (const std::string& nid, bool include_self) const {
    std::vector<Node*> ancestors;
    Node* node = get_node(nid);
    if (node==NULL) {
        return ancestors;
    }    
    if (include_self) {
        ancestors.push_back(node);
    }
    while (node->parent != NULL) {
        ancestors.push_back(node->parent);
        node = node->parent;
    }
    return ancestors;
}

void Mutation_Annotated_Tree::Tree::remove_node_helper (std::string nid, bool move_level) { 
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
                if (child->clade == "") {
                    child->clade = curr_parent->clade;
                }
                child->parent = curr_parent->parent;
                child->level = curr_parent->parent->level + 1;
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
static void reassign_level_helper(Mutation_Annotated_Tree::Node* root){
    for(auto c:root->children){
        c->level=root->level+1;
        reassign_level_helper(c);
    }
}
void Mutation_Annotated_Tree::Tree::reassign_level(){
    root->level=1;
    reassign_level_helper(root);
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
        curr_node->level = curr_node->parent->level + 1;
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;
    
    if (nid == "") {
        if (root == NULL) {
            return traversal;
        }
        nid = root->identifier;
    }

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

static void depth_first_expansion_helper(Mutation_Annotated_Tree::Node* node, std::vector<Mutation_Annotated_Tree::Node*>& vec, size_t& index) {
    #ifdef DETAIL_DEBUG_NO_LOOP
    assert(std::find(vec.begin(),vec.end(),node)==vec.end());
    #endif
    vec.push_back(node);
    assert(vec.size()-1==index);
    node->index=index;
    index++;
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec,index);
    }
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::depth_first_expansion(Mutation_Annotated_Tree::Node* node) const {
    TIMEIT();
    std::vector<Node*> traversal;
    if (node == NULL) {
        node = root;
    }
    size_t index=0;
    depth_first_expansion_helper(node, traversal,index);
    return traversal;
}

size_t Mutation_Annotated_Tree::Tree::get_parsimony_score() {
    size_t score = 0;
    auto dfs = depth_first_expansion();
    for (auto n: dfs) {
        score += n->mutations.size();
    }
    return score;
}

void Mutation_Annotated_Tree::Tree::condense_leaves(std::vector<std::string> missing_samples) {
    if (condensed_nodes.size() > 0) {
        fprintf(stderr, "WARNING: tree contains condensed nodes. It may be condensed already!\n");
    }

    auto tree_leaves = get_leaves_ids();
    for (auto l1_id: tree_leaves) {
        std::vector<Node*> polytomy_nodes;

        auto l1 = get_node(l1_id);
        if (l1 == NULL) {
            continue;
        }
        if (std::find(missing_samples.begin(), missing_samples.end(), l1->identifier) != missing_samples.end()) {
            continue;
        }
        if (l1->mutations.size() > 0) {
            continue;
        }

        for (auto l2: l1->parent->children) {
                if (std::find(missing_samples.begin(), missing_samples.end(), l2->identifier) != missing_samples.end()) {
                    continue;
                }
            if (l2->is_leaf() && (get_node(l2->identifier) != NULL) && (l2->mutations.size() == 0)) {
                polytomy_nodes.push_back(l2);
            }
        }

        if (polytomy_nodes.size() > 1) {
            std::string new_node_name = "node_" + std::to_string(1+condensed_nodes.size()) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
            
            auto curr_node = get_node(l1->identifier);
            auto new_node = create_node(new_node_name, curr_node->parent, l1->branch_length);

            new_node->clear_mutations();
            
            condensed_nodes[new_node_name] = std::vector<std::string>(polytomy_nodes.size());

            for (size_t it = 0; it < polytomy_nodes.size(); it++) {
                condensed_nodes[new_node_name][it] = polytomy_nodes[it]->identifier;
                remove_node(polytomy_nodes[it]->identifier, false);
            }
        }
    }
}

void Mutation_Annotated_Tree::Tree::uncondense_leaves() {
    for (size_t it = 0; it < condensed_nodes.size(); it++) {
        auto cn = condensed_nodes.begin();
        std::advance(cn, it);

        auto n = get_node(cn->first);
        auto par = (n->parent != NULL) ? n->parent : n;

        size_t num_samples = cn->second.size();

        if (num_samples > 0) {
            rename_node(n->identifier, cn->second[0]);
        }

        for (size_t s = 1; s < num_samples; s++) {
            create_node(cn->second[s], par, n->branch_length);
        }
    }
    condensed_nodes.clear();
    condensed_leaves.clear();
}
// Merge nodes that have no mutations comparing to parent into parent node 
void Mutation_Annotated_Tree::Tree::collapse_tree() {
    auto bfs = breadth_first_expansion();

    for (size_t idx = 1; idx < bfs.size(); idx++) {
        auto node = bfs[idx];
        auto mutations = node->mutations;        
        if (mutations.size() == 0) {
            auto parent = node->parent;
            auto children = node->children;
            for (auto child: children) {
                move_node(child->identifier, parent->identifier);
            }
        }
        //If internal node has one child, the child can be moved up one level
        else if (node->children.size() == 1) {
            auto child = node->children.front();
            auto parent = node->parent;
            for (auto m: mutations) {
                child->add_mutation(m);
            }
            move_node(child->identifier, parent->identifier);
        }
    }
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_tree_copy(const Mutation_Annotated_Tree::Tree& tree, const std::string& identifier) {
    TIMEIT();
    auto root = tree.root;
    if (identifier != "") {
        root = tree.get_node(identifier);
    }
    
    Tree copy = create_tree_from_newick_string (get_newick_string(tree, root, true, true));

    std::vector<Node*> dfs1;
    std::vector<Node*> dfs2;

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, 2),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
              if (k==0) {
                dfs1 = tree.depth_first_expansion(root);
              }
              else {
                dfs2 = copy.depth_first_expansion();
              }
            }
            }, ap);


    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs1.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
              auto n1 = dfs1[k];
              auto n2 = dfs2[k];
              n2->clade = n1->clade;
              for (auto m: n1->mutations) {
                Mutation m2 = m.copy();
                n2->add_mutation(m2);
                }
              }
            }, ap);

    size_t num_condensed_nodes = static_cast<size_t>(tree.condensed_nodes.size());
    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_condensed_nodes),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t idx = r.begin(); idx < r.end(); idx++) {
               auto cn = tree.condensed_nodes.begin(); 
               std::advance(cn, idx);
               copy.condensed_nodes.insert(std::pair<std::string, std::vector<std::string>>(cn->first, std::vector<std::string>(cn->second.size())));
               for (size_t k = 0; k < cn->second.size(); k++) {
                  copy.condensed_nodes[cn->first][k] = cn->second[k];
                  copy.condensed_leaves.insert(cn->second[k]);
               }
            }
    }, ap);

    return copy;
}



void Mutation_Annotated_Tree::exchange(Node *branch1, Node *branch2){
    //Make sure they are not root
    assert(!branch1->is_root());
    assert(!branch2->is_root());
    Node* const branch1_old_parent=branch1->parent;
    Node* const branch2_old_parent=branch2->parent;

    //locate branch 1 among the children of branch 1, and replace it with branch2
    auto iter=std::find(branch1_old_parent->children.begin(),branch1_old_parent->children.end(),branch1);
    assert(iter!=branch1_old_parent->children.end());//ehh, why not here...
    *iter=branch2;
    //change its parent
    branch2->parent=branch1_old_parent;

    //the same for branch2
    iter=std::find(branch2_old_parent->children.begin(),branch2_old_parent->children.end(),branch2);
    assert(iter!=branch2_old_parent->children.end());
    *iter=branch1;
    branch1->parent=branch2_old_parent;
}
// Get the last common ancestor of two node identifiers. Return NULL if does not
// exist
Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::LCA (const Mutation_Annotated_Tree::Tree& tree, const std::string& nid1, const std::string& nid2) {
    TIMEIT();
    Node* ret = NULL;

    if ((tree.get_node(nid1) == NULL) || (tree.get_node(nid2) == NULL)) {
        return ret;
    }

    for (auto anc: tree.rsearch(nid1)) {
        ret = anc;
        if (tree.is_ancestor(anc->identifier, nid2)) {
            return ret;
        }
    }

    return ret;
}

// Extract the subtree consisting of the specified set of samples. This routine
// maintains the internal node names of the input tree. Mutations are copied
// from the tree such that the path of mutations from root to the sample is
// same as the original tree.
Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_subtree (const Mutation_Annotated_Tree::Tree& tree, const std::vector<std::string>& samples) {
    TIMEIT();
    Tree subtree;

    // Set of leaf and internal nodes corresponding to the subtree
    tbb::concurrent_unordered_set<Node*> subtree_nodes;
    // Maintain a set of all ancestors of a sample for each sample
    std::vector<tbb::concurrent_unordered_set<Node*>> all_ancestors(samples.size());

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k){
               subtree_nodes.insert(tree.get_node(samples[k]));
               for (auto anc: tree.rsearch(samples[k], true)) {
                   all_ancestors[k].insert(anc);
               }
            }
    }, ap);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
            [&](tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
               for (size_t j=i+1; j<samples.size(); ++j){
                   for (auto anc: tree.rsearch(samples[i], true)) {
                      if (all_ancestors[j].find(anc) != all_ancestors[j].end()) {
                         subtree_nodes.insert(anc);
                         break;
                      }
                   }
               }
            }
    }, ap);
    
    auto dfs = tree.depth_first_expansion();
    std::stack<Node*> last_subtree_node;
    for (auto n: dfs) {
        // If the node is in subtree_nodes, it should be added to the subtree
        if (subtree_nodes.find(n) != subtree_nodes.end()) {
            Node* subtree_parent = NULL;
            if (last_subtree_node.size() > 0) {
                while (!tree.is_ancestor(last_subtree_node.top()->identifier, n->identifier)) {
                    last_subtree_node.pop();
                }
                subtree_parent = last_subtree_node.top();
            }
            // Add as root of the subtree
            if (subtree_parent == NULL) {
                Node* new_node = subtree.create_node(n->identifier);
                
                std::vector<Node*> root_to_node = tree.rsearch(n->identifier); 
                std::reverse(root_to_node.begin(), root_to_node.end());
                root_to_node.emplace_back(n);

                for (auto curr: root_to_node) {
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            // Add to the parent identified
            else {
                Node* new_node = subtree.create_node(n->identifier, subtree_parent->identifier);

                auto par_to_node = tree.rsearch(n->identifier, true);
                std::reverse(par_to_node.begin(), par_to_node.end());

                for (auto curr: par_to_node) {
                    if (curr->clade != "") {
                        new_node->clade = curr->clade;
                    }
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            last_subtree_node.push(n);
        }
    }

    subtree.curr_internal_node = tree.curr_internal_node;

    return subtree;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Node::add_child(Node *new_child) {
    if (is_leaf()) {
        std::string old_name(identifier);
        tree->rename_node(identifier, std::to_string(++tree->curr_internal_node));
        Mutation_Annotated_Tree::Node* sample_node=tree->create_node(old_name,this);
        children.push_back(new_child);
        return sample_node;
    }
    new_child->parent = this;
    children.push_back(new_child);
    return nullptr;
}
Mutation_Annotated_Tree::Node::Node(const Node &other, Node *parent, Tree *tree)
    : level(other.level), branch_length(other.branch_length),
      identifier(other.identifier), parent(parent),
      mutations(other.mutations), tree(tree) {
    children.reserve(other.children.size());
    for (auto c : other.children) {
        children.push_back(new Node(*c, this, tree));
    }
    tree->all_nodes.emplace(other.identifier,this);
}

static void write_newick_with_mutations_helper(FILE *f,Mutation_Annotated_Tree::Node* node){
    if(!node->children.empty()){
        fputc('(',f);
        bool first=true;
        for(auto child:node->children){
            if (first) {
                first=false;
            }else{
                fputc(',',f);
            }
            write_newick_with_mutations_helper(f, child);
        }
        fputc(')',f);
    }
    fputs(node->identifier.c_str(),f);
    for(const Mutation_Annotated_Tree::Mutation& m:node->mutations){
        fprintf(f, "_%c%d%c",Mutation_Annotated_Tree::get_nuc(m.par_nuc),m.position,Mutation_Annotated_Tree::get_nuc(m.mut_nuc));
    }
    fprintf(f,":%zu",node->mutations.size());
}
void Mutation_Annotated_Tree::Tree::write_newick_with_mutations(FILE *f){
    write_newick_with_mutations_helper(f, root);
    fclose(f);
}
#ifdef MEMDEBUG
void Node::delete_this(){
    for(Node* n:children){
        n->delete_this();
    }
    delete this;
}
void Mutation_Annotated_Tree::Tree::delete_nodes(){
    root->delete_this();
}
#endif