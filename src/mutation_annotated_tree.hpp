#ifndef MUTATION_ANNOTATED_TREE
#define MUTATION_ANNOTATED_TREE
#include <cstddef>
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
#include <tbb/tbb.h>
#include "parsimony.pb.h"
#include "Instrumentor.h"

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
        mutable int8_t mut_nuc;
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

        bool operator==(const Mutation& other) const{
            if(other.mut_nuc!=mut_nuc) return false;
            if(other.par_nuc!=par_nuc) return false;
            if(other.is_missing!=is_missing) return false;
            if(other.position!=position) return false;
            assert(other.chrom==chrom);
            assert(other.ref_nuc==ref_nuc);
            return true;
        }
        inline bool is_masked() const {
            return (position < 0);
        }
        inline std::string get_string() const {
            if (is_masked()) {
                return "MASKED";
            }
            else {
                return get_nuc(par_nuc) + std::to_string(position) + get_nuc(mut_nuc);
            }
        }
    };

#ifndef matToVCF
    class Mutations_Collection{
        std::vector<Mutation> mutations;
        public:
        static const char NO_DUPLICATE=-1;
        static const char KEEP_SELF=1;
        static const char KEEP_OTHER=0;
        static const char MERGE=2;
        static const char INVERT_MERGE=3;
        typedef std::vector<Mutation>::iterator iterator;
        size_t size() const {return mutations.size();}
        void swap(Mutations_Collection& in){
            mutations.swap(in.mutations);
        }
        Mutation& operator[] (size_t idx){
            return mutations[idx];
        }
        void operator=(const Mutations_Collection & other){
            mutations=other.mutations;
        }
        iterator begin()  {
            return mutations.begin();
        }
        iterator end(){
            return mutations.end();
        }
        void clear(){
            mutations.clear();
        }
        void reserve(size_t n){
            mutations.reserve(n);
        }
        void push_back(Mutation& m){
            assert(m.position>mutations.back().position);
            mutations.push_back(m);
        }
        iterator find_next( int pos) {
            #ifndef NDEBUG
            int lastPos=0;
            #endif
            auto iter = mutations.begin();
            for (; iter < mutations.end(); iter++) {
                //check sorting in debug mode
                #ifndef NDEBUG
                assert(lastPos<iter->position);
                lastPos=iter->position;
                #endif
                if (iter->position>=pos) {
                    break;
                }
            }
            assert(iter == mutations.begin()||(iter-1)->position<pos);
            assert(iter==mutations.end()||(iter)->position>=pos);
            return iter;
        }

        iterator find(int position) {
            auto iter=find_next(position);
            if(iter!=mutations.end()&&iter->position>position){
                return mutations.end();
            }
            return iter;
        }
        
        Mutations_Collection reverse() const{
            Mutations_Collection result;
            result.mutations.reserve(mutations.size());
            for(Mutation m:mutations){
                auto t=m.mut_nuc;
                m.mut_nuc=m.par_nuc;
                m.par_nuc=t;
                result.mutations.push_back(m);
            }
            return result;
        }
        iterator find(const Mutation& mut) {
            return find(mut.position);
        }
        bool remove(int pos){
            auto iter=find(pos);
            if(iter==mutations.end()){
                return false;
            }
            mutations.erase(iter);
            return true;
        }
        /**
         * @brief Merge other mutations into this mutation set
         * 
         * @param other sorted vector of mutations to merge in
         * @param keep_self if two mutations are at the same position 0: keep other, 1: keep self, -1: impossible, throw an error. 2: Merge as successor, 3. Inverted merge
         */
        void merge_out(const Mutations_Collection& other,Mutations_Collection& out, char keep_self) const{
//used for checking whether the two vectors are sorted while merging
#ifndef NDEBUG
#define mutation_vector_check_order(newly_inserted) \
assert(last_pos_inserted<(newly_inserted));\
last_pos_inserted=(newly_inserted);
            int last_pos_inserted=-1; 
#else
#define mutation_vector_check_order(newly_inserted)
#endif
            out.mutations.reserve(other.mutations.size()+mutations.size());
            auto other_iter=other.mutations.begin();
            for(auto this_mutation:mutations){
                while (other_iter != other.mutations.end()&&other_iter->position<this_mutation.position) {
                    mutation_vector_check_order(other_iter->position);
                    out.mutations.push_back(*other_iter);
                    if(keep_self==INVERT_MERGE){
                        auto temp=out.mutations.back().mut_nuc;
                        out.mutations.back().mut_nuc=out.mutations.back().par_nuc;
                        out.mutations.back().par_nuc=temp;
                    }
                    other_iter++;
                }
                if (other_iter==other.mutations.end()||
                    this_mutation.position<other_iter->position) {
                    mutation_vector_check_order(this_mutation.position);
                    out.mutations.push_back(this_mutation);
                }else{
                    mutation_vector_check_order(this_mutation.position);

                    assert(this_mutation.position==other_iter->position);
                    switch(keep_self){
                        case NO_DUPLICATE: assert(false);
                        case KEEP_OTHER: out.mutations.push_back(*other_iter); break;
                        case KEEP_SELF:out.mutations.push_back(this_mutation); break;
                        case MERGE:assert(other_iter->par_nuc==this_mutation.mut_nuc);
                            if(other_iter->mut_nuc!=this_mutation.ref_nuc){
                                out.mutations.push_back(*other_iter);
                                out.mutations.back().par_nuc=this_mutation.par_nuc;
                            }
                            break;
                        default:
                            assert(keep_self==INVERT_MERGE);
                            assert(other_iter->mut_nuc==this_mutation.mut_nuc);
                            if(other_iter->par_nuc!=this_mutation.ref_nuc){
                                out.mutations.push_back(this_mutation);
                                out.mutations.back().mut_nuc=other_iter->par_nuc;
                            }
                            break;
                    }
                    other_iter++;
                }
            }
            while (other_iter<other.mutations.end()) {
                mutation_vector_check_order(other_iter->position);
                out.mutations.push_back(*other_iter);
                other_iter++;
            }
        }

        void merge(const Mutations_Collection& other, char keep_self){
            Mutations_Collection new_set;
            merge_out(other, new_set, keep_self);
            mutations.swap(new_set.mutations);
        }

        void set_difference(const Mutations_Collection &other,
                            Mutations_Collection &this_unique,
                            Mutations_Collection &other_unique,
                            Mutations_Collection &common) {
            this_unique.mutations.reserve(mutations.size());
            other_unique.mutations.reserve(other.mutations.size());
            common.mutations.reserve(std::min(mutations.size(), other.mutations.size()));
            auto other_iter = other.mutations.begin();
#ifndef NDEBUG
            int last_pos_inserted = -1;
#endif
            // merge sort again
            for (auto this_mutation : mutations) {
                while (other_iter != other.mutations.end()&&other_iter->position < this_mutation.position) {
                    mutation_vector_check_order(other_iter->position);
                    other_unique.mutations.push_back(*other_iter);
                    other_iter++;
                }
                if (other_iter == other.mutations.end() ||
                    this_mutation.position < other_iter->position) {
                    mutation_vector_check_order(this_mutation.position);
                    this_unique.mutations.push_back(this_mutation);
                } else {
                    mutation_vector_check_order(this_mutation.position);
                    assert(this_mutation.position == other_iter->position);

                    if (other_iter->mut_nuc == this_mutation.mut_nuc) {
                        assert(other_iter->par_nuc == this_mutation.par_nuc);
                        common.mutations.push_back(this_mutation);
                    } else {
                        this_unique.mutations.push_back(this_mutation);
                        other_unique.mutations.push_back(*other_iter);
                    }
                    other_iter++;
                }
            }
            while (other_iter < other.mutations.end()) {
                mutation_vector_check_order(other_iter->position);
                other_unique.mutations.push_back(*other_iter);
                other_iter++;
            }
            assert(this_unique.size() + other_unique.size() +
                       2 * common.size() ==
                   mutations.size() + other.mutations.size());
        }
        void batch_find(Mutations_Collection &target) {
            // yet another merge sort
#ifndef NDEBUG
            int last_pos_inserted = -1;
#endif
            auto target_iter = target.mutations.begin();
            for (auto this_mutation : mutations) {
                if (target_iter == target.mutations.end()) {
                    return;
                }
                while (target_iter->position < this_mutation.position) {
                    mutation_vector_check_order(target_iter->position);
                    target_iter++;
                    if (target_iter == target.mutations.end()) {
                        return;
                    }
                }
                if (this_mutation.position < target_iter->position) {
                    mutation_vector_check_order(this_mutation.position);
                } else {
                    mutation_vector_check_order(this_mutation.position);
                    assert(this_mutation.position == target_iter->position);
                    *target_iter = this_mutation;
                    target_iter++;
                }
            }
            if (target_iter != target.mutations.end()) {
                mutation_vector_check_order(target_iter->position);
            }
        }
        bool insert(const Mutation &mut, char keep_self = -1) {
            assert(mut.par_nuc!=mut.mut_nuc);
            auto iter=find_next(mut.position);
            if(iter!=mutations.end()&&iter->position==mut.position){
                assert(keep_self!=-1);
                if(!keep_self){
                    if(*iter==mut){
                        return false;
                    }
                    *iter=mut;
                }
                return false;
            }
            mutations.insert(iter,mut);
            return true;
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
            size_t index; //index in dfs pre-order
            bool is_leaf();
            bool is_root();

            Node();
            Node(std::string id, float l);
            Node(std::string id, Node* p, float l);
#ifndef matToVCF   
            void add_mutation(Mutation& mut){
                mutations.insert(mut);
            }
            // HACK: try to identify whether a node is a sample from name
            bool not_sample() const{
                for (auto c : this->identifier) {
                    if (!std::isdigit(c))
                        return false;
                }
                return true;
            }
            void clear_mutations(){
                mutations.clear();
            }
#else
            void add_mutation(Mutation mut);
            void clear_mutations();
#endif
    };

    class Tree {
        private:
            void remove_node_helper (std::string nid, bool move_level);
            std::unordered_map <std::string, Node*> all_nodes;
        public:
            Tree() {
                root = NULL;
                all_nodes.clear();
            }

            Tree (Node* n);
            std::vector<Node*> new_nodes;
            size_t max_level;

            Node* root;
            tbb::concurrent_unordered_map<std::string, std::vector<std::string>> condensed_nodes;
            tbb::concurrent_unordered_set<std::string> condensed_leaves;

            size_t curr_internal_node;
            size_t get_max_level ();
            void reassign_level();
            void rename_node(std::string old_nid, std::string new_nid);
            std::vector<Node*> get_leaves(std::string nid="");
            std::vector<std::string> get_leaves_ids(std::string nid="");
            size_t get_num_leaves(Node* node=NULL);
            Node* create_node (std::string const& identifier, float branch_length = -1.0);
            Node* create_node (std::string const& identifier, Node* par, float branch_length = -1.0);
            Node* create_node (std::string const& identifier, std::string const& parent_id, float branch_length = -1.0);
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
    
    std::string get_newick_string(Tree& T, bool b1, bool b2, bool b3=false, bool b4=false);
    std::string get_newick_string(Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
    void write_newick_string (std::stringstream& ss, Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
    Tree create_tree_from_newick (std::string filename);
    Tree create_tree_from_newick_string (std::string newick_string);
    void string_split(std::string const& s, char delim, std::vector<std::string>& words);
    void string_split(std::string s, std::vector<std::string>& words);

    Tree load_mutation_annotated_tree (std::string filename);
    void save_mutation_annotated_tree (Tree tree, std::string filename);
    
    Tree get_tree_copy(Tree tree, std::string identifier="");

    // Exchange 2 branches of the same tree that are not root (not checked to be the same tree)
    void exchange(Node*  branch1, Node* branch2);
}
static bool check_grand_parent(const Mutation_Annotated_Tree::Node* node,const Mutation_Annotated_Tree::Node* grand_parent){
    const Mutation_Annotated_Tree::Node* cur=node;
    while (cur) {
        if(cur==grand_parent) return true;
        cur=cur->parent;
    }
    return false;
}




static char one_hot_to_two_bit(char arg) {return 31-__builtin_clz((unsigned int)arg);}
static char two_bit_to_one_hot(char arg) {return 1<<(arg);}
#endif