#pragma once
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <cassert>
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/task_group.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include "parsimony.pb.h"
#include "Instrumentor.h"

// Forward declaration of structs from usher_graph
struct Missing_Sample;

#if SAVE_PROFILE == 1
#  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__);
#else
#  define TIMEIT()
#endif

namespace Mutation_Annotated_Tree {
int8_t get_nuc_id (char nuc);
int8_t get_nuc_id (std::vector<int8_t> nuc_vec);
char get_nuc (int8_t nuc_id);
int8_t get_nt (int8_t nuc_id);
std::vector<int8_t> get_nuc_vec (char nuc);
std::vector<int8_t> get_nuc_vec_from_id (int8_t nuc_id);

// WARNING: chrom is currently ignored!
// position < 0 implies masked mutations i.e. mutations that exist but
// details are unknown
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
    inline Mutation copy() const {
        Mutation m;
        m.chrom = chrom;
        m.position = position;
        m.ref_nuc = ref_nuc;
        m.par_nuc = par_nuc;
        m.mut_nuc = mut_nuc;
        m.is_missing = is_missing;
        return m;
    }
    Mutation () {
        chrom = "";
        is_missing = false;
    }
    inline bool is_masked() const {
        return (position < 0);
    }
    inline std::string get_string() const {
        if (is_masked()) {
            return "MASKED";
        } else {
            return get_nuc(par_nuc) + std::to_string(position) + get_nuc(mut_nuc);
        }
    }
};

class Node {
  public:
    size_t level;
    float branch_length;
    std::string identifier;
    std::vector<std::string> clade_annotations;
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
    void clear_annotations();
};

class Tree {
  private:
    void remove_node_helper (std::string nid, bool move_level);
    void depth_first_expansion_helper(Node* node, std::vector<Node*>& vec) const;
    std::unordered_map <std::string, Node*> all_nodes;
  public:
    Tree() {
        root = NULL;
        curr_internal_node = 0;
        all_nodes.clear();
    }

    std::string new_internal_node_id() {
        return "node_" + std::to_string(++curr_internal_node);

    }

    Node* root;
    tbb::concurrent_unordered_map<std::string, std::vector<std::string>> condensed_nodes;
    tbb::concurrent_unordered_set<std::string> condensed_leaves;

    size_t curr_internal_node;
    size_t get_max_level () const;
    size_t get_num_annotations () const;
    void rename_node(std::string old_nid, std::string new_nid);
    std::vector<Node*> get_leaves(std::string nid="");
    std::vector<std::string> get_leaves_ids(std::string nid="");
    size_t get_num_leaves(Node* node=NULL);
    Node* create_node (std::string const& identifier, float branch_length = -1.0, size_t num_annotations=0);
    Node* create_node (std::string const& identifier, Node* par, float branch_length = -1.0);
    Node* create_node (std::string const& identifier, std::string const& parent_id, float branch_length = -1.0);
    Node* get_node (std::string identifier) const;
    bool is_ancestor (std::string anc_id, std::string nid) const;
    std::vector<Node*> rsearch (const std::string& nid, bool include_self = false) const;
    std::string get_clade_assignment (const Node* n, int clade_id, bool include_self = true) const;
    void remove_node (std::string nid, bool move_level);
    void remove_single_child_nodes();
    void move_node (std::string source, std::string destination, bool move_level=true);
    std::vector<Node*> breadth_first_expansion(std::string nid="");
    std::vector<Node*> depth_first_expansion(Node* node=NULL) const;

    size_t get_parsimony_score();

    void condense_leaves(std::vector<std::string> = std::vector<std::string>());
    void uncondense_leaves();
    void collapse_tree();
    void rotate_for_display(bool reverse = false);
    void rotate_for_consistency();
};

std::string get_newick_string(const Tree& T, bool b1, bool b2, bool b3=false, bool b4=false);
std::string get_newick_string(const Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
void write_newick_string (std::stringstream& ss, const Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
Tree create_tree_from_newick (std::string filename);
Tree create_tree_from_newick_string (std::string newick_string);
void string_split(std::string const& s, char delim, std::vector<std::string>& words);
void string_split(std::string s, std::vector<std::string>& words);
Mutation* mutation_from_string(const std::string& mut_string);

Tree load_mutation_annotated_tree (std::string filename);
void save_mutation_annotated_tree (Tree tree, std::string filename);

Tree get_tree_copy(const Tree& tree, const std::string& identifier="");

Node* LCA (const Tree& tree, const std::string& node_id1, const std::string& node_id2);
Tree get_subtree (const Tree& tree, const std::vector<std::string>& samples, bool keep_clade_annotations=false);
void get_random_single_subtree (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string outdir, size_t subtree_size, size_t tree_idx = 0, bool use_tree_idx = false, bool retain_original_branch_len = false);
void get_random_sample_subtrees (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string outdir, size_t subtree_size, size_t tree_idx = 0, bool use_tree_idx = false, bool retain_original_branch_len = false);
void get_sample_mutation_paths (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string mutation_paths_filename);
void clear_tree(Tree& tree);

void read_vcf (Mutation_Annotated_Tree::Tree* T, std::string &vcf_filename, std::vector<Missing_Sample>& missing_samples, bool create_new_mat);
}

