#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>

class Node {
    public:
        size_t level;
        float branch_length;
        std::string identifier;
        Node* parent;
        std::vector<Node*> children;

        bool is_leaf();
        bool is_root();
        Node();
        Node(std::string id, float l);
        Node(std::string id, Node* p, float l);
};

class Tree {
    private:
        void remove_node_helper (std::string nid, bool move_level);
        void depth_first_expansion_helper(Node* node, std::vector<Node*>& vec);
        std::unordered_map <std::string, Node*> all_nodes;
    public:
        Tree() {
            max_level = 0;
            root = NULL;
            all_nodes.clear();
        }
        
        Tree (Node* n);

        size_t max_level;
        Node* root;

        size_t curr_internal_node;
        size_t get_max_level ();
        void rename_node(std::string old_nid, std::string new_nid);
        std::vector<Node*> get_leaves();
        std::vector<Node*> get_leaves(std::string nid);
        size_t get_num_leaves(Node* node);
        size_t get_num_leaves();
        void create_node (std::string identifier, float branch_length = -1.0);
        void create_node (std::string identifier, std::string parent_id, float branch_length = -1.0);
        Node* get_node (std::string identifier);
        bool is_ancestor (std::string anc_id, std::string nid);
        std::vector<Node*> rsearch (std::string nid);
        void remove_node (std::string nid, bool move_level);
        void move_node (std::string source, std::string destination);
        std::vector<Node*> breadth_first_expansion();
        std::vector<Node*> breadth_first_expansion(std::string nid);
        std::vector<Node*> depth_first_expansion();
        std::vector<Node*> depth_first_expansion(Node* node);
};

namespace TreeLib {
std::string get_newick_string(Tree& T, bool b1, bool b2);
std::string get_newick_string(Tree& T, Node* node, bool b1, bool b2);
Tree create_tree_from_newick (std::string filename);
Tree create_tree_from_newick_string (std::string newick_string);
void split(std::string s, char delim, std::vector<std::string>& words);
void split(std::string s, std::vector<std::string>& words);
}

