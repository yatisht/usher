#ifndef MUTATION_ANNOTATED_TREE
#define MUTATION_ANNOTATED_TREE
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
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
static uint8_t one_hot_to_two_bit(uint8_t arg) {
    return 31-__builtin_clz((unsigned int)arg);
}
static uint8_t two_bit_to_one_hot(uint8_t arg) {
    return 1<<(arg);
}
class nuc_one_hot;
class nuc_2bit {
    uint8_t nuc;
  public:
    nuc_2bit(uint8_t nuc):nuc(nuc) {}
    operator uint8_t() const {
        return nuc;
    }
    operator nuc_one_hot() const;
};
//wrapper of one hot encoded state
class nuc_one_hot {
    uint8_t nuc;
  public:
    nuc_one_hot() {
#ifndef NDEBUG
        nuc=0xff;
#endif
    }

    nuc_one_hot(uint8_t nuc,bool skip_check=false):nuc(nuc) {
        //assert(!(nuc&0xf0)||skip_check);
    }
    bool is_invalid() const {
        return nuc==0xff;
    }
    operator uint8_t() const {
        //assert(!(nuc&0xf0));
        return nuc;
    }

    uint8_t get_nuc_no_check()const {
        return nuc;
    }
    nuc_2bit to_2bit() const {
        assert(!(nuc&0xf0));
        assert(nuc);
        assert(__builtin_popcount(nuc)==1);
        return one_hot_to_two_bit(nuc);
    }

    bool is_ambiguous()const {
        //assert(nuc);
        //assert(!(nuc&0xf0));
        return __builtin_popcount(nuc)!=1;
    }

    nuc_one_hot choose_first()const {
        //assert(!(nuc&0xf0));
        uint8_t ret=1<<__builtin_ctz(nuc);
        //assert(ret&nuc);
        return ret;
    }
};

inline nuc_2bit::operator nuc_one_hot() const {
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
    uint8_t boundary1_all_major_allele; //boundary 1 alleles are alleles with allele count one less than major allele count
    uint8_t decrement_increment_effect;//Decrement which allele will increase parsimony score, then increment which allele may decrease parsimony score
    //uint8_t child_muts;//left child state then right child state for binary nodes
    static tbb::concurrent_unordered_map<std::string, uint8_t> chromosome_map;
    static std::mutex ref_lock;//reference nuc are stored in a separate vector, need to be locked when adding new mutations
  public:
    static std::vector<std::string> chromosomes;// chromosome index to name map
    static std::vector<nuc_one_hot> refs;
    void set_boundary_one_hot(nuc_one_hot boundary1) {
        boundary1_all_major_allele=(boundary1_all_major_allele&0xf)|(boundary1<<4);
    }
    nuc_one_hot get_sensitive_decrement()const {
        return decrement_increment_effect>>4;
    }
    nuc_one_hot get_sensitive_increment()const {
        return decrement_increment_effect&0xf;
    }
    void set_sensitive_change(nuc_one_hot decrement,nuc_one_hot increment) {
        decrement_increment_effect=(decrement<<4)|increment;
    }
    Mutation(const std::string& chromosome,int position,nuc_one_hot mut,nuc_one_hot par,nuc_one_hot tie,nuc_one_hot ref=0);

    Mutation(int pos):position(pos),chrom_idx(0),par_mut_nuc(0),boundary1_all_major_allele(0) {}

    Mutation(uint8_t chrom_idx,int pos,uint8_t par_nuc,uint8_t mut_nuc):position(pos),chrom_idx(chrom_idx),par_mut_nuc((par_nuc<<4)|mut_nuc),boundary1_all_major_allele(mut_nuc) {}
    struct ignored {};
    Mutation(uint8_t chrom_idx,int pos,uint8_t par_nuc,uint8_t mut_nuc,ignored):position(pos),chrom_idx(chrom_idx),par_mut_nuc((par_nuc<<4)|mut_nuc),boundary1_all_major_allele(0xf) {}

    bool same_chrom(const Mutation& other) const {
        return chrom_idx==other.chrom_idx;
    }

    const std::string& get_chromosome() const {
        return chromosomes[chrom_idx];
    }

    uint8_t get_chromIdx()const {
        return chrom_idx;
    }

    //order defined only by position, not chromidx yet
    inline bool operator< (const Mutation& m) const {
        return ((*this).position < m.position);
    }

    inline bool operator<= (const Mutation& m) const {
        return ((*this).position <= m.position);
    }
    /*
            inline Mutation copy() const {
                Mutation m(*this);
                return m;
            }
     */
    Mutation():position(-1),chrom_idx(0),par_mut_nuc(0),boundary1_all_major_allele(0) {
    }

    nuc_one_hot get_ref_one_hot() const {
        return refs[position];
    }

    nuc_one_hot get_mut_one_hot() const {
        return par_mut_nuc&0xf;
    }

    void set_mut_one_hot(nuc_one_hot new_mut) {
        par_mut_nuc&=0xf0;
        par_mut_nuc|=new_mut;
    }

    nuc_one_hot get_par_one_hot() const {
        return par_mut_nuc>>4;
    }

    nuc_one_hot get_all_major_allele() const {
        return boundary1_all_major_allele&0xf;
    }

    nuc_one_hot get_boundary1_one_hot() const {
        return boundary1_all_major_allele>>4;
    }

    void set_auxillary(nuc_one_hot all_major_allele,nuc_one_hot boundary1) {
        //assert(all_major_allele&get_mut_one_hot());
        boundary1_all_major_allele=all_major_allele|(boundary1<<4);
    }

    void set_par_mut(nuc_one_hot new_par,nuc_one_hot new_mut) {
        par_mut_nuc=(new_par<<4)|new_mut;
    }
    void set_par_one_hot(nuc_one_hot new_par) {
        par_mut_nuc&=0xf;
        par_mut_nuc|=(new_par<<4);
    }
    //valid when it par_nuc!=mut_nuc
    bool is_valid() const {
        return (par_mut_nuc^(par_mut_nuc<<4))&0xf0;
    }

    int get_position() const {
        return position;
    }
    //position and all states equal
    bool operator==(const Mutation& other) const {
        if(other.par_mut_nuc!=par_mut_nuc) return false;
        if(other.position!=position) return false;
        if (other.boundary1_all_major_allele!=boundary1_all_major_allele) {
            return false;
        }
        //assert(other.chrom_idx==chrom_idx);
        return true;
    }
    inline bool is_masked() const {
        return (position < 0);
    }
    inline std::string get_string() const {
        if (is_masked()) {
            return "MASKED";
        } else {
            return get_nuc(get_par_one_hot()) + std::to_string(position) + get_nuc(get_mut_one_hot());
        }
    }
};
//Wraper of a vector of mutation, sort of like boost::flatmap
class Mutations_Collection {
  public:
    std::vector<Mutation> mutations;
    //Falgs for merging two mutation vector
    static const char NO_DUPLICATE=-1;
    static const char KEEP_SELF=1;
    static const char KEEP_OTHER=0;
    static const char MERGE=2;
    static const char INVERT_MERGE=3;
    typedef std::vector<Mutation>::iterator iterator;
    typedef std::vector<Mutation>::const_iterator const_iterator;
    typedef Mutation value_type;
    size_t size() const {
        return mutations.size();
    }
    Mutations_Collection() {}
    void swap(Mutations_Collection& in) {
        mutations.swap(in.mutations);
    }
    Mutation& operator[] (size_t idx) {
        return mutations[idx];
    }
    const Mutation& operator[] (size_t idx) const {
        return mutations[idx];
    }
    bool no_valid_mutation()const;
    //For finishing up fitch sankoff, fill itself with a sorted input
    void refill(std::vector<Mutation>& in) {
        mutations=std::move(in);
        std::sort(mutations.begin(),mutations.end());
    }
    void operator=(const Mutations_Collection & other) {
        mutations=other.mutations;
    }
    Mutation& back() {
        return mutations.back();
    }
    const Mutation* data() const {
        return mutations.data();
    }
    iterator begin()  {
        return mutations.begin();
    }
    iterator end() {
        return mutations.end();
    }
    const_iterator begin() const {
        return mutations.begin();
    }
    const_iterator end() const {
        return mutations.end();
    }
    void clear() {
        mutations.clear();
    }
    bool empty() const {
        return mutations.empty();
    }
    void reserve(size_t n) {
        mutations.reserve(n);
    }
    void push_back(const Mutation& m) {
        assert(mutations.empty()||m.get_position()>mutations.back().get_position());
        mutations.push_back(m);
    }
    //Find next mutation with position greater or equal to pos
    iterator find_next(int pos);

    iterator find(int position) {
        auto iter=find_next(position);
        if(iter!=mutations.end()&&iter->get_position()>position) {
            return mutations.end();
        }
        return iter;
    }
    //invert par_nuc and mut_nuc for all mutations inside the vector
    Mutations_Collection reverse() const {
        Mutations_Collection result;
        result.mutations.reserve(mutations.size());
        for(Mutation m:mutations) {
            auto t=m.get_mut_one_hot();
            m.set_mut_one_hot(m.get_par_one_hot());
            m.set_par_one_hot(t);
            result.mutations.push_back(m);
        }
        return result;
    }
    int count_valid_mutations()const;
    iterator find(const Mutation& mut) {
        return find(mut.get_position());
    }
    bool remove(int pos) {
        auto iter=find(pos);
        if(iter==mutations.end()) {
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
    void merge(const Mutations_Collection& other, char keep_self) {
        Mutations_Collection new_set;
        merge_out(other, new_set, keep_self);
        mutations.swap(new_set.mutations);
    }
    //For outputing usher compatible pb, remove all mutations with par_nuc==mut_nuc
    void remove_invalid();
    //for leaf nodes and nodes with one children, their boundary1 allele are inferred
    void remove_boundary_only() {
        mutations.erase(std::remove_if(mutations.begin(), mutations.end(), [](const Mutation& mut) {
            return mut.get_par_one_hot()==mut.get_all_major_allele();
        }),mutations.end());
    }
    void set_difference(const Mutations_Collection &other,
                        Mutations_Collection &this_unique,
                        Mutations_Collection &other_unique,
                        Mutations_Collection &common) const;
    bool insert(const Mutation &mut, char keep_self = -1) {
        //assert(mut.get_par_one_hot()!=mut.get_mut_one_hot());
        auto iter=find_next(mut.get_position());
        if(iter!=mutations.end()&&iter->get_position()==mut.get_position()) {
            assert(keep_self!=NO_DUPLICATE);
            if(keep_self==KEEP_OTHER) {
                *iter=mut;
            }
            assert(keep_self==KEEP_SELF||keep_self==KEEP_OTHER);
            return false;
        }
        mutations.insert(iter,mut);
        return true;
    }
};
typedef std::vector<std::pair<int,int>> ignored_t;
class Node {
#ifdef LOAD
  public:
#endif
    std::atomic_uint8_t changed;
    static const uint8_t SELF_CHANGED_MASK=1;
    static const uint8_t ANCESTOR_CHANGED_MASK=2;
    static const uint8_t DESCENDENT_CHANGED_MASK=4;
    static const uint8_t SELF_MOVED_MASK=8;
  public:
    float branch_length;
    std::string identifier;
    std::vector<std::string> clade_annotations;
    Node* parent;
    std::vector<Node*> children;
    Mutations_Collection mutations;
    ignored_t ignore;
    //Mutations_Collection boundary_mutations;
    size_t dfs_index; //index in dfs pre-order
    size_t dfs_end_index; //index in dfs pre-order
    size_t bfs_index; //index in bfs
    size_t level;
    //size_t last_searched_arcs;
    bool have_masked;
    bool is_leaf() const;
    bool is_root();
    Node();
    Node(std::string id, float l);
    Node(std::string id, Node* p, float l);

    Node(const Node& other, Node* parent,Tree* tree,bool copy_mutation=true);
    bool add_mutation(Mutation& mut) {
        return mutations.insert(mut,Mutations_Collection::KEEP_OTHER);
    }
    void set_self_changed() {
        changed|=SELF_CHANGED_MASK;
    }
    void set_ancestor_changed() {
        changed|=ANCESTOR_CHANGED_MASK;
    }
    void set_descendent_changed() {
        changed|=DESCENDENT_CHANGED_MASK;
    }
    void set_self_moved() {
        changed|=SELF_MOVED_MASK;
    }
    void clear_changed() {
        changed=0;
    }
    bool get_self_changed() const {
        return (changed&SELF_CHANGED_MASK);
    }
    bool get_ancestor_changed() const {
        return (changed&ANCESTOR_CHANGED_MASK);
    }
    bool get_descendent_changed()const {
        return (changed&DESCENDENT_CHANGED_MASK);
    }
    bool get_self_moved()const {
        return changed&SELF_MOVED_MASK;
    }
    bool have_change_in_neighbor()const {
        return changed;
    }
    void clear_mutations() {
        mutations.clear();
    }
    void clear_annotations() {
        clade_annotations.clear();
    }
    bool no_valid_mutation()const {
        return mutations.no_valid_mutation();
    }
    void populate_ignored_range();
    //refill mutation with the option of filtering invalid mutations
    template<typename iter_t>
    void refill(iter_t begin,iter_t end,size_t size=0,bool retain_invalid=true) {
        std::vector<Mutation> mutations;
        mutations.reserve(size);
        for(; begin<end; begin++) {
            assert(begin->get_position()>0);
            if (retain_invalid||begin->is_valid()) {
                mutations.push_back(*begin);
            }
        }
        this->mutations.refill(mutations);
    }
    Node* add_child(Node* new_child,Mutation_Annotated_Tree::Tree* tree);
    void delete_this();
};

class Tree {
  public:
    typedef  tbb::concurrent_unordered_map<std::string, std::vector<std::string>> condensed_node_t;
    std::unordered_map <std::string, Node*> all_nodes;
    Tree() {
        root = NULL;
        all_nodes.clear();
    }

    Tree (Node* n);

    std::vector<Node*> new_nodes;
    size_t max_level;

    Node* root;
    condensed_node_t condensed_nodes;

    size_t curr_internal_node;

    void rename_node(std::string old_nid, std::string new_nid);
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
    void load_from_newick(const std::string& newick_string,bool use_internal_node_label=false);
    void condense_leaves(std::vector<std::string> = std::vector<std::string>());
    void uncondense_leaves();
    std::vector<Node*> get_leaves() const;
    void populate_ignored_range();
    size_t get_max_level();
    friend class Node;
    void MPI_send_tree() const;
    void MPI_receive_tree();
    void delete_nodes();
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
}
bool check_grand_parent(const Mutation_Annotated_Tree::Node* node,const Mutation_Annotated_Tree::Node* grand_parent);

nuc_one_hot get_parent_state(Mutation_Annotated_Tree::Node* ancestor,int position);
#endif
