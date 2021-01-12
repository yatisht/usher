#ifndef MUTATION_ANNOTATED_TREE
#define MUTATION_ANNOTATED_TREE
#include <cstddef>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
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
#include <tbb/tbb.h>
#include "parsimony.pb.h"
namespace Mutation_Annotated_Tree {
    int8_t get_nuc_id (char nuc);
    int8_t get_nuc_id (std::vector<int8_t> nuc_vec);
    char get_nuc (int8_t nuc_id);
    int8_t get_nt (int8_t nuc_id);
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
    };

#ifndef matToVCF
    class Mutations_Collection{
        std::vector<Mutation> mutations;
        public:
        size_t size() const {return mutations.size();}
        std::vector<Mutation>::const_iterator begin() const {
            return mutations.begin();
        }
        std::vector<Mutation>::const_iterator end() const{
            return mutations.end();
        }
        void clear(){
            mutations.clear();
        }
        std::vector<Mutation>::const_iterator find_next(const int pos) const{
            #ifndef NDEBUG
            int lastPos=0;
            #endif
            for(auto iter=mutations.begin();iter<mutations.end();iter++){
                //check sorting in debug mode
                #ifndef NDEBUG
                assert(lastPos<iter->position);
                lastPos=iter->position;
                #endif
                if (iter->position>=pos) {
                    return iter;
                }
            }
            return mutations.end();
        }
        std::vector<Mutation>::const_iterator find(const Mutation& mut) const{
            auto iter=find_next(mut.position);
            if(iter!=mutations.end()&&iter->position>mut.position){
                return mutations.end();
            }
            return iter;
        }

        /**
         * @brief Merge other mutations into this mutation set
         * 
         * @param other sorted vector of mutations to merge in
         * @param keep_self if two mutations are at the same position 0: keep other, 1: keep self, -1: impossible, throw an error.
         */
        void merge(const Mutations_Collection other, char keep_self){
//used for checking whether the two vectors are sorted while merging
#ifndef NDEBUG
#define mutation_vector_check_order(newly_inserted) \
assert(last_pos_inserted<(newly_inserted));\
last_pos_inserted=(newly_inserted);
            int last_pos_inserted=-1; 
#else
#define mutation_vector_check_order(newly_inserted)
#endif
            std::vector<Mutation> new_set;
            new_set.reserve(other.mutations.size()+mutations.size());
            auto other_iter=other.mutations.begin();
            for(auto this_mutation:mutations){
                while (other_iter->position<this_mutation.position) {
                    mutation_vector_check_order(other_iter->position);
                    new_set.push_back(*other_iter);
                    other_iter++;
                }
                if (other_iter==other.mutations.end()||
                    this_mutation.position<other_iter->position) {
                    mutation_vector_check_order(this_mutation.position);
                    new_set.push_back(this_mutation);
                }else{
                    mutation_vector_check_order(this_mutation.position);

                    assert(this_mutation.position==other_iter->position);
                    assert(keep_self!=-1);
                    if (keep_self) {
                        new_set.push_back(this_mutation);
                    }else {
                        new_set.push_back(*other_iter);
                    }
                    other_iter++;
                }
            }
            while (other_iter<other.mutations.end()) {
                mutation_vector_check_order(other_iter->position);
                new_set.push_back(*other_iter);
            }
            mutations.swap(new_set);
        }

        void insert(const Mutation& mut){
            auto iter=find_next(mut.position);
            //Not supposed to insert 2 mutations at the same position
            assert(iter==mutations.end()||iter->position>mut.position);
            mutations.insert(iter,mut);
        }
        void finalize(){
            mutations.shrink_to_fit();
        }
    };
#else
    typedef std::vector<Mutation> Mutations_Collection;
#endif
    class Node {
        public:
            size_t level;
            float branch_length;
            std::string identifier;
            Node* parent;
            std::vector<Node*> children;
            Mutations_Collection mutations;
            bool is_new; //whether this node is a new sample

            bool is_leaf();
            bool is_root();

            Node();
            Node(std::string id, float l);
            Node(std::string id, Node* p, float l);
#ifndef matToVCF   
            void add_mutation(Mutation& mut){
                mutations.insert(mut);
            }
            void clear_mutations(){
                mutations.clear();
            }
#endif
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
            Node* create_node (std::string identifier, std::string parent_id, float branch_length = -1.0);
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
    
    std::string get_newick_string(Tree& T, bool b1, bool b2, bool b3=false);
    std::string get_newick_string(Tree& T, Node* node, bool b1, bool b2, bool b3=false);
    Tree create_tree_from_newick (std::string filename);
    Tree create_tree_from_newick_string (std::string newick_string);
    void string_split(std::string s, char delim, std::vector<std::string>& words);
    void string_split(std::string s, std::vector<std::string>& words);

    Tree load_mutation_annotated_tree (std::string filename);
    void save_mutation_annotated_tree (Tree tree, std::string filename);

    Tree get_tree_copy(Tree tree);
    // Exchange 2 branches of the same tree that are not root (not checked to be the same tree)
    void exchange(Node*  branch1, Node* branch2);
}

#endif