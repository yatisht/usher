#ifndef PROFITABLE_MOVES_ENUMERATOR
#define PROFITABLE_MOVES_ENUMERATOR
#include "../tree_rearrangement_internal.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <cstdint>
#include <vector>
#include <array>
#include "../stack_allocator.hpp"
typedef std::vector<std::array<std::vector<unsigned int>, 4>> mutated_node_dfs_idx_t;
extern mutated_node_dfs_idx_t mutated_node_dfs_idx;
namespace MAT = Mutation_Annotated_Tree;
//Class for recording change in major allele set
class Mutation_Count_Change {
    int position;
    uint8_t chromIdx;
    nuc_one_hot decremented_allele;
    nuc_one_hot incremented_allele;
    nuc_one_hot par_state;

  public:
    //bool from_invalid_mutation;
    static const char VALID_MASK=1;
    static const char END_MASK=2;
    Mutation_Count_Change() {
        position=INT_MAX;
        decremented_allele=0;
        incremented_allele=0;
        //from_invalid_mutation=false;
    }
    //copy over position, original state, par state from a mutation
    Mutation_Count_Change(const Mutation_Count_Change& child_mut_count,nuc_one_hot new_major_allele) {
        position=child_mut_count.position;
        chromIdx=child_mut_count.chromIdx;
        auto ori_state=child_mut_count.par_state;
        par_state=child_mut_count.par_state;
        incremented_allele=new_major_allele&(~ori_state);
        decremented_allele=ori_state&(~new_major_allele);
        //from_invalid_mutation=false;
    }
    Mutation_Count_Change(const MAT::Mutation &pos,nuc_one_hot decremented,nuc_one_hot incremented) {
        position = pos.get_position();
        chromIdx = pos.get_chromIdx();
        par_state=pos.get_par_one_hot();
        set_change(decremented, incremented);
        //from_invalid_mutation=false;
    }
    template<typename T>
    Mutation_Count_Change(const T &pos,nuc_one_hot decremented,nuc_one_hot incremented) {
        position = pos.get_position();
        chromIdx = pos.get_chromIdx();
        par_state=pos.get_par_state();
        set_change(decremented, incremented);
        //from_invalid_mutation=false;
    }
    Mutation_Count_Change(const Mutation_Count_Change &pos,nuc_one_hot decremented,nuc_one_hot incremented):Mutation_Count_Change(pos) {
        set_change(decremented, incremented);
        //from_invalid_mutation=false;
    }
    int get_position() const {
        return position;
    }
    nuc_one_hot get_decremented() const { //assert(decremented_allele!=0xff);
        return decremented_allele;
    }
    nuc_one_hot get_incremented() const { //assert(incremented_allele!=0xff);
        return incremented_allele;
    }
    void set_change(nuc_one_hot decremented, nuc_one_hot incremented) {
        decremented_allele = decremented;
        incremented_allele = incremented;
    }
    nuc_one_hot get_par_state() const {
        return par_state;
    }
    uint8_t get_chromIdx()const {
        return chromIdx;
    }
    //parsimony score change (number of children able to follow major allele) if the parent node is not sensitive at this loci
    int get_default_change_internal()const {
        if(par_state&incremented_allele) {
            //able to follow major allele
            //assert(was_valid);
            return -1;
        }
        if(par_state&decremented_allele) {
            //no longer able to follow major allele
            //assert(!was_valid);
            return 1;
        }
        return 0;
    }
    //Also account for change in number of children, if is terminal
    int get_default_change_terminal()const {
        if(incremented_allele&&(!(par_state&incremented_allele))) {
            //add new children, and not of parent allele state
            return 1;
        }
        if((!(decremented_allele&par_state))&&decremented_allele) {
            //removed children with valid mutation
            return -1;
        }
        return 0;
    }
    bool operator<(const Mutation_Count_Change &rhs) const {
        return position < rhs.position;
    }
    void set_par_nuc(nuc_one_hot par_nuc) {
        par_state=par_nuc;
    }
};

static bool operator<(const MAT::Mutation &lhs, const Mutation_Count_Change &rhs) {
    return lhs.get_position() < rhs.get_position();
}
typedef std::vector<Mutation_Count_Change,stack_allocator<Mutation_Count_Change>> Mutation_Count_Change_Collection_FIFO;

typedef std::vector<Mutation_Count_Change> Mutation_Count_Change_Collection;
typedef Mutation_Count_Change_Collection::const_iterator Mut_Change_Iter;
//convience class for tracking iterator and its end position
template <typename value_type> class range {
    const value_type* curr;
    const value_type* const end;

  public:
    template<typename T>
    range(const T &container) : curr(container.data()), end(container.data()+container.size()) {}
    operator bool() const {
        return curr != end;
    }
    const value_type* operator->() const {
        return curr;
    }
    const value_type &operator*() const {
        return *curr;
    }
    void operator++() {
        curr++;
    }
    void operator++(int) {
        curr++;
    }
    size_t size()const {
        return end-curr;
    }
};
struct Reachable {
    bool reachable_change;
    bool always_search;
};
void find_moves_bounded(MAT::Node* src,output_t& out,int search_radius,bool do_drift,Reachable
#ifdef CHECK_BOUND
                        ,counters& count
#endif
                       );
struct node_info {
    size_t dfs_idx;
    size_t level;
    bool operator<(const node_info& other)const {
        return dfs_idx<other.dfs_idx;
    }
    uint8_t base;
};
struct range_tree_node {
    uint32_t dfs_start_idx;
    uint32_t dfs_end_idx;
    std::array<LEVEL_T, 4> min_level;
    IDX_TREE_IDX_T children_start_idx;
    IDX_TREE_IDX_T parent_idx;
    IDX_TREE_IDX_T level;
    //uint8_t child_size;
    //uint8_t dist_from_end;
};
struct range_tree {
    std::vector<uint32_t> end_idxes;
    std::vector<range_tree_node> nodes;
    IDX_TREE_IDX_T find_idx(const MAT::Node* node,IDX_TREE_IDX_T& last_probe_idx,IDX_TREE_IDX_T probe_start_idx=0) const ;
    IDX_TREE_IDX_T find_idx(const MAT::Node* node) const {
        IDX_TREE_IDX_T ignored;
        return find_idx(node,ignored,0);
    }
};
bool output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,int radius_left);
extern std::vector<range_tree> addable_idxes;
void check_parsimony_score_change_above_LCA(MAT::Node *start_node, int &parsimony_score_change,
        Mutation_Count_Change_Collection &parent_added,
        Mutation_Count_Change_Collection &parent_of_parent_added);
bool dst_branch(const MAT::Node *LCA,
                const range<Mutation_Count_Change> &mutations,
                int &parsimony_score_change, MAT::Node *this_node,
                Mutation_Count_Change_Collection &parent_added,int src_side_max_improvement
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,std::vector<Mutation_Count_Change_Collection> &debug_from_dst
#endif
               );
int check_move_profitable_LCA(
    MAT::Node *src, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change,
    const MAT::Node* last_src_branch_node_below_LCA,
    output_t &output,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> &debug_above_LCA,
    const MAT::Tree* tree
#endif
);
int check_move_profitable_dst_not_LCA(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const range<Mutation_Count_Change>  &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change, output_t &output,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> debug_from_src,
    const MAT::Tree* tree
#endif
);
#endif
