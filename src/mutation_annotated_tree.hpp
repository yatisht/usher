#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>
#include "parsimony.pb.h"

namespace Mutation_Annotated_Tree {
    int8_t get_nuc_id (char nuc);
    int8_t get_nuc_id (std::vector<int8_t> nuc_vec);
    char get_nuc (int8_t nuc_id);
    std::vector<int8_t> get_nuc_vec (char nuc);
    std::vector<int8_t> get_nuc_vec_from_id (int8_t nuc_id);

    // WARNING: chrom is currently ignored!
    struct Mutation {
        std::string chrom;
        int position;
        int8_t ref_nuc;
        int8_t par_nuc;
        int8_t mut_nuc;
        bool is_missing;
        inline bool operator< (const Mutation& m) const {
            return ((*this).position < m.position);
        }
    };

    class Node {
        public:
            size_t level;
            float branch_length;
            std::string identifier;
            Node* parent;
            std::vector<Node*> children;
            std::vector<Mutation> mutations;

            bool is_leaf();
            bool is_root();

            Node();
            Node(std::string id, float l);
            Node(std::string id, Node* p, float l);
            
            void add_mutation(Mutation mut);
            void clear_mutations();
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
            std::unordered_map<std::string, std::vector<std::string>> condensed_nodes;
            std::unordered_set<std::string> condensed_leaves;

            size_t curr_internal_node;
            size_t get_max_level ();
            void rename_node(std::string old_nid, std::string new_nid);
            std::vector<Node*> get_leaves(std::string nid="");
            size_t get_num_leaves(Node* node=NULL);
            void create_node (std::string identifier, float branch_length = -1.0);
            void create_node (std::string identifier, std::string parent_id, float branch_length = -1.0);
            Node* get_node (std::string identifier);
            bool is_ancestor (std::string anc_id, std::string nid);
            std::vector<Node*> rsearch (std::string nid);
            void remove_node (std::string nid, bool move_level);
            void move_node (std::string source, std::string destination);
            std::vector<Node*> breadth_first_expansion(std::string nid="");
            std::vector<Node*> depth_first_expansion(Node* node=NULL);

            size_t get_parsimony_score();
            void condense_leaves(std::vector<std::string> = std::vector<std::string>());
            void uncondense_leaves();
            void collapse_tree();
    };
    
    std::string get_newick_string(Tree& T, bool b1, bool b2);
    std::string get_newick_string(Tree& T, Node* node, bool b1, bool b2);
    Tree create_tree_from_newick (std::string filename);
    Tree create_tree_from_newick_string (std::string newick_string);
    void string_split(std::string s, char delim, std::vector<std::string>& words);
    void string_split(std::string s, std::vector<std::string>& words);

    Tree load_mutation_annotated_tree (std::string filename);
    void save_mutation_annotated_tree (Tree tree, std::string filename);

    Tree get_tree_copy(Tree tree);
}

