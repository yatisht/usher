#ifndef MUTATION_ANNOTATED_TREE
#define MUTATION_ANNOTATED_TREE
#include <cstddef>
#include <mutex>
#include <sys/types.h>
#include <unordered_map>
#include <string>
#include <utility>
#include <vector>
#include <cassert>
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_unordered_map.h"

#if SAVE_PROFILE == 1
#  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__);
#else
#  define TIMEIT()
#endif

static uint8_t one_hot_to_two_bit(uint8_t arg) {return 31-__builtin_clz((unsigned int)arg);}
static uint8_t two_bit_to_one_hot(uint8_t arg) {return 1<<(arg);}

class nuc_one_hot;
class nuc_2bit{
    uint8_t nuc;
    public:
    nuc_2bit(uint8_t nuc):nuc(nuc){}
    operator uint8_t() const{
        return nuc;
    }
    operator nuc_one_hot() const;
};
class nuc_one_hot{
    uint8_t nuc;
    public:
    nuc_one_hot():nuc(0xff){}

    nuc_one_hot(uint8_t nuc,bool skip_check=false):nuc(nuc){
        assert(!(nuc&0xf0)||skip_check);
    }
    bool is_invalid() const{
        return nuc==0xff;
    }
    operator uint8_t() const{
        assert(!(nuc&0xf0));
        return nuc;
    }

    uint8_t get_nuc_no_check()const{
        return nuc;
    }
    nuc_2bit to_2bit() const{
        assert(!(nuc&0xf0));
        assert(nuc);
        assert(__builtin_popcount(nuc)==1);
        return one_hot_to_two_bit(nuc);
    }

    bool is_ambiguous()const{
        assert(nuc);
        assert(!(nuc&0xf0));
        return __builtin_popcount(nuc)!=1;
    }

    nuc_one_hot choose_first()const{
        assert(!(nuc&0xf0));
        uint8_t ret=1<<__builtin_ctz(nuc);
        assert(ret&nuc);
        return ret;
    }
};

inline nuc_2bit::operator nuc_one_hot() const{
    return two_bit_to_one_hot(nuc);
}

namespace Mutation_Annotated_Tree {
    class Tree;
    class Node;
    int8_t get_nuc_id (char nuc);
    int8_t get_nuc_id (std::vector<int8_t> nuc_vec);
    char get_nuc (int8_t nuc_id);
    int8_t get_nt (int8_t nuc_id);
    std::vector<int8_t> get_nuc_vec (char nuc);
    std::vector<int8_t> get_nuc_vec_from_id (int8_t nuc_id);

    // WARNING: chrom is currently ignored!
    // position < 0 implies masked mutations i.e. mutations that exist but
    // details are unknown
    class Mutation;
    class Mutation {
        #ifdef LOAD
        public:
        #endif
        int position;
        uint8_t chrom_idx;
        uint8_t par_mut_nuc;
        uint8_t boundary1_tie;
        uint8_t boundary2_flag;
        static tbb::concurrent_unordered_map<std::string, uint8_t> chromosome_map;
        static std::mutex ref_lock;
        public:
        static std::vector<std::string> chromosomes;
        static std::vector<nuc_one_hot> refs;
        void set_boundary_one_hot(nuc_one_hot boundary1){
            boundary1_tie=(boundary1_tie&0xf)|(boundary1<<4);
        }
        
        Mutation(const std::string& chromosome,int position,nuc_one_hot mut,nuc_one_hot par,nuc_one_hot tie,nuc_one_hot boundary1, nuc_one_hot boundary2,nuc_one_hot ref=0);

        Mutation(int pos):position(pos),chrom_idx(0),par_mut_nuc(0),boundary1_tie(0),boundary2_flag(0){}
        bool same_chrom(const Mutation& other) const{
            return chrom_idx==other.chrom_idx;
        }

        const std::string& get_chromosome() const{
            return chromosomes[chrom_idx];
        }

        const uint8_t get_chromIdx()const{
            return chrom_idx;
        }

        bool get_boundary_by_2()const{
            return boundary2_flag>>4;
        }

        inline bool operator< (const Mutation& m) const {
            return ((*this).position < m.position);
        }

        inline bool operator<= (const Mutation& m) const {
            return ((*this).position <= m.position);
        }

        inline Mutation copy() const {
            Mutation m(*this);
            return m;
        }
        
        Mutation():position(-1),chrom_idx(0),par_mut_nuc(0),boundary1_tie(0),boundary2_flag(0){
        }

/*
        void set_minor_allele(nuc_one_hot minor_allele){
            ref_minor&=0xf0;
            ref_minor|=minor_allele;
        }
        uint8_t get_flag() const{
            return flag2;
        }

        void set_flag(uint8_t flag){
            this->flag2=flag;
        }
        */
        nuc_one_hot get_ref_one_hot() const{
            return refs[position];
        }

        nuc_one_hot get_mut_one_hot() const{
            return par_mut_nuc&0xf;
        }

        void set_mut_one_hot(nuc_one_hot new_mut) {
            par_mut_nuc&=0xf0;
            par_mut_nuc|=new_mut;
        }

        nuc_one_hot get_par_one_hot() const{
            return par_mut_nuc>>4;
        }

        nuc_one_hot get_tie_one_hot() const{
            return boundary1_tie&0xf;
        }

        nuc_one_hot get_boundary1_one_hot() const{
            return boundary1_tie>>4;
        }

        nuc_one_hot get_boundary2_one_hot() const{
            return boundary2_flag>>4;
        }

        void set_auxillary(nuc_one_hot tie,nuc_one_hot boundary1,nuc_one_hot boundary2){
            boundary1_tie=tie|(boundary1<<4);
            boundary2_flag=(boundary2_flag&0xf)|(boundary2<<4);
        }

        void set_par_mut(nuc_one_hot new_par,nuc_one_hot new_mut){
            par_mut_nuc=(new_par<<4)|new_mut;
        }
        void set_par_one_hot(nuc_one_hot new_par){
            par_mut_nuc&=0xf;
            par_mut_nuc|=(new_par<<4);
        }


        bool is_valid() const{
            return (par_mut_nuc^(par_mut_nuc<<4))&0xf0;
        }


        int get_position() const {
            return position;
        }

        bool operator==(const Mutation& other) const{
            if(other.par_mut_nuc!=par_mut_nuc) return false;
            if(other.position!=position) return false;
            if (other.boundary1_tie!=boundary1_tie) {
                return false;
            }
            if (other.boundary2_flag!=boundary2_flag) {
                return false;
            }
            assert(other.chrom_idx==chrom_idx);
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
                return get_nuc(get_par_one_hot()) + std::to_string(position) + get_nuc(get_mut_one_hot());
            }
        }
    };

    class Mutations_Collection{
        public:
        std::vector<Mutation> mutations;
        std::vector<char> remove_or_replace;
        std::vector<Mutation> new_inserts;

        typedef std::mutex mutex_type;
        //mutex_type mutex;
        static const char NO_DUPLICATE=-1;
        static const char KEEP_SELF=1;
        static const char KEEP_OTHER=0;
        static const char MERGE=2;
        static const char INVERT_MERGE=3;
        typedef std::vector<Mutation>::iterator iterator;
        typedef std::vector<Mutation>::const_iterator const_iterator;
        size_t size() const {return mutations.size();}
        Mutations_Collection(){}
        void swap(Mutations_Collection& in){
            mutations.swap(in.mutations);
        }
        Mutation& operator[] (size_t idx){
            return mutations[idx];
        }
        const Mutation& operator[] (size_t idx) const{
            return mutations[idx];
        }
        void refill(std::vector<Mutation>& in){
            mutations=std::move(in);
            std::sort(mutations.begin(),mutations.end());
        }
        void operator=(const Mutations_Collection & other){
            mutations=other.mutations;
        }
        Mutation& back(){
            return mutations.back();
        }
        iterator begin()  {
            return mutations.begin();
        }
        iterator end(){
            return mutations.end();
        }
        const_iterator begin() const {
            return mutations.begin();
        }
        const_iterator end() const{
            return mutations.end();
        }
        void clear(){
            mutations.clear();
        }
        bool empty() const{
            return mutations.empty();
        }
        void reserve(size_t n){
            mutations.reserve(n);
        }
        void push_back(const Mutation& m){
            assert(mutations.empty()||m.get_position()>mutations.back().get_position());
            mutations.push_back(m);
        }
        iterator find_next(int pos);

        iterator find(int position) {
            auto iter=find_next(position);
            if(iter!=mutations.end()&&iter->get_position()>position){
                return mutations.end();
            }
            return iter;
        }
        Mutations_Collection reverse() const{
            Mutations_Collection result;
            result.mutations.reserve(mutations.size());
            for(Mutation m:mutations){
                auto t=m.get_mut_one_hot();
                m.set_mut_one_hot(m.get_par_one_hot());
                m.set_par_one_hot(t);
                result.mutations.push_back(m);
            }
            return result;
        }
        iterator find(const Mutation& mut) {
            return find(mut.get_position());
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
        void merge_out(const Mutations_Collection& other,Mutations_Collection& out, char keep_self) const;
        void merge(const Mutations_Collection& other, char keep_self){
            Mutations_Collection new_set;
            merge_out(other, new_set, keep_self);
            mutations.swap(new_set.mutations);
        }
        void remove_invalid();
        void set_difference(const Mutations_Collection &other,
                            Mutations_Collection &this_unique,
                            Mutations_Collection &other_unique,
                            Mutations_Collection &common) const;
        void batch_find(Mutations_Collection &target);
        bool insert(const Mutation &mut, char keep_self = -1) {
            //assert(mut.get_par_one_hot()!=mut.get_mut_one_hot());
            auto iter=find_next(mut.get_position());
            if(iter!=mutations.end()&&iter->get_position()==mut.get_position()){
                assert(keep_self!=NO_DUPLICATE);
                if(keep_self==KEEP_OTHER){
                    *iter=mut;
                }
                assert(keep_self==KEEP_SELF||keep_self==KEEP_OTHER);
                return false;
            }
            mutations.insert(iter,mut);
            return true;
        }
        void finalize();
    };

    class Node {
        public:
            float branch_length;
            std::string identifier;
            std::vector<std::string> clade_annotations;
            Node* parent;
            std::vector<Node*> children;
            Mutations_Collection mutations;
            //Mutations_Collection boundary_mutations;
            size_t dfs_index; //index in dfs pre-order
            size_t bfs_index; //index in bfs
            bool is_leaf() const;
            bool is_root();
            Tree* tree;
            Node();
            Node(std::string id, float l);
            Node(std::string id, Node* p, float l);
            Node(Tree* tree):Node("",-1,tree){
            }
            
            Node(const std::string& id, float l,Tree* tree):Node(id,nullptr,l,tree){}
            
            Node(const std::string& id, Node* p,float l,Tree* tree):branch_length(l),identifier(id),parent(p),tree(tree){}
            Node(const Node& other, Node* parent,Tree* tree,bool copy_mutation=true);
            bool add_mutation(Mutation& mut){
                return mutations.insert(mut,Mutations_Collection::KEEP_OTHER);
            }
            void finalize(){
                mutations.finalize();
            }
            void clear_mutations(){
                mutations.clear();
            }
            void clear_annotations() {
                clade_annotations.clear();
            }
            template<typename iter_t>
            void refill(iter_t begin,iter_t end,size_t size=0,bool retain_invalid=true){
                std::vector<Mutation> mutations;
                mutations.reserve(size);
                for(;begin<end;begin++){
                    if (retain_invalid||begin->is_valid()) {
                        mutations.push_back(*begin);
                    }
                }
                mutations.shrink_to_fit();
                this->mutations.refill(mutations);
            }
            Node* add_child(Node* new_child);
            #ifdef MEMDEBUG
            void delete_this();
            #endif
    };

    class Tree {
        private:
            void remove_node_helper (std::string nid, bool move_level);
        public:
            typedef  tbb::concurrent_unordered_map<std::string, std::vector<std::string>> condensed_node_t;
            std::unordered_map <std::string, Node*> all_nodes;
            Tree() {
                root = NULL;
                all_nodes.clear();
            }

            Tree (Node* n);

            std::vector<Node*> new_nodes;
            //size_t max_level;

            Node* root;
            condensed_node_t condensed_nodes;
            tbb::concurrent_unordered_set<std::string> condensed_leaves;

            size_t curr_internal_node;

            void rename_node(std::string old_nid, std::string new_nid);
            std::vector<Node*> get_leaves(std::string nid="");
            std::vector<std::string> get_leaves_ids(std::string nid="");
            size_t get_num_leaves(Node* node=NULL);
            Node* create_node (std::string const& identifier, float branch_length = -1.0, size_t num_annotations=0); 
            Node* create_node (std::string const& identifier, Node* par, float branch_length = -1.0);
            Node* create_node (std::string const& identifier, std::string const& parent_id, float branch_length = -1.0);
            Node* get_node (const std::string& identifier) const {
                auto iter=all_nodes.find(identifier);
                if (iter != all_nodes.end()) {
                    return iter->second;
                }
                return NULL;
            }
            Node* get_node_c_str (char* identifier) const;
            void remove_node (std::string nid, bool move_level);
            void move_node (std::string source, std::string destination);
            std::vector<Node*> breadth_first_expansion(std::string nid="");
            std::vector<Node*> depth_first_expansion(Node* node=NULL) const;

            size_t get_parsimony_score();

            size_t get_num_annotations() const {
                size_t ret = 0;
                if (root != NULL) {
                    ret = root->clade_annotations.size();
                }
                return ret;
            }
            void save_detailed_mutations(const std::string& path)const;
            void load_detatiled_mutations(const std::string& path);
            bool is_ancestor (std::string anc_id, std::string nid) const;
            std::vector<Mutation_Annotated_Tree::Node*> rsearch (const std::string& nid, bool include_self) const ;
            void load_from_newick(const std::string& newick_string,bool use_internal_node_label=false);
            //void write_newick_with_mutations(FILE* f);
            void condense_leaves(std::vector<std::string> = std::vector<std::string>());
            void uncondense_leaves();
            void collapse_tree();
            
            void finalize();
            friend class Node;
            #ifdef MEMDEBUG
            void delete_nodes();
            #endif
    };
    
    std::string get_newick_string(const Tree& T, bool b1, bool b2, bool b3=false, bool b4=false);
    std::string get_newick_string(const Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
    void write_newick_string (std::stringstream& ss, const Tree& T, Node* node, bool b1, bool b2, bool b3=false, bool b4=false);
    Tree create_tree_from_newick (std::string filename);
    Tree create_tree_from_newick_string (std::string newick_string);
    void string_split(std::string const& s, char delim, std::vector<std::string>& words);
    void string_split(std::string s, std::vector<std::string>& words);

    Tree load_mutation_annotated_tree (std::string filename);
    void save_mutation_annotated_tree (const Tree& tree, std::string filename);
    
    Tree get_tree_copy(const Tree& tree, const std::string& identifier="");
    
    Node* LCA (const Tree& tree, const std::string& node_id1, const std::string& node_id2);
    Tree get_subtree (const Tree& tree, const std::vector<std::string>& samples);
    void remove_child(Node* child_to_remove,std::vector<Node*>& removed_nodes);
}
bool check_grand_parent(const Mutation_Annotated_Tree::Node* node,const Mutation_Annotated_Tree::Node* grand_parent);

template<>
struct std::hash<Mutation_Annotated_Tree::Mutation> {
    size_t operator()(const Mutation_Annotated_Tree::Mutation &in) const {
        return in.get_position();
    }
};
nuc_one_hot get_parent_state(Mutation_Annotated_Tree::Node* ancestor,int position);
#endif