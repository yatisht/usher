#include <boost/program_options.hpp>
#include <src/mutation_annotated_tree.hpp>
#include <signal.h>
#include <iostream>

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
struct Mapper_Info {
    const MAT::Mutation* begin;
    const MAT::Mutation* end;
    unsigned int sibling_start_idx;
    unsigned short level;
    bool is_leaf;
};
namespace MAT = Mutation_Annotated_Tree;
struct Ripples_Mapper_Mut {
    static const unsigned short NULL_MUT_IDX = -1;
    int position;
    unsigned short mut_idx;
    char curr_mut;
    char dest_mut;
    inline bool valid() const {
        return curr_mut != dest_mut;
    }
    Ripples_Mapper_Mut()
        : position(INT_MAX), mut_idx(NULL_MUT_IDX), curr_mut(0), dest_mut(0) {}
    Ripples_Mapper_Mut(const MAT::Mutation &mut, size_t idx)
        : position(mut.position), mut_idx(idx), curr_mut(mut.ref_nuc),
          dest_mut(mut.mut_nuc) {
    }
    Ripples_Mapper_Mut(const MAT::Mutation &mut)
        : position(mut.position), mut_idx(NULL_MUT_IDX), curr_mut(mut.mut_nuc),
          dest_mut(mut.ref_nuc) {
        //assert(mut.par_nuc == mut.ref_nuc);
    }
    Ripples_Mapper_Mut(const Ripples_Mapper_Mut &mut, char new_mut)
        : position(mut.position), mut_idx(mut.mut_idx), curr_mut(new_mut),
          dest_mut(mut.dest_mut) {}
};
class Mut_Count_t {
    unsigned short internal;

  public:
    Mut_Count_t() {}
    void set(unsigned short count_before_exclusive, bool is_self) {
        internal = (is_self << 15) | count_before_exclusive;
    }
    unsigned short count_before_exclusive() const {
        return internal & 0x7fff;
    }
    unsigned short is_self_counted() const {
        return (internal >> 15);
    }
    unsigned short count_before_inclusive() const {
        return count_before_exclusive() + is_self_counted();
    }
};
typedef std::vector<Mut_Count_t> Mut_Count_Out_t;
typedef std::vector<std::vector<Ripples_Mapper_Mut>> Mut_Out_t;

struct Ripples_Mapper_Output_Interface {
    Mut_Count_Out_t mut_count_out;
    //Mut_Out_t mut_out;
    std::vector<char> is_sibling;
};

struct Pruned_Sample {
    MAT::Node* sample_name;
    std::vector<MAT::Mutation> sample_mutations;
    std::unordered_set<uint32_t> positions;

    // Assumes mutations are added in reverse chrono order
    void add_mutation(MAT::Mutation mut);
    Pruned_Sample() {}

    Pruned_Sample(MAT::Node* name);
};

struct Recomb_Node {
    const MAT::Node* node;
    int node_parsimony;
    int parsimony;
    bool is_sibling;
    Recomb_Node() {
        node_parsimony = -1;
        parsimony = -1;
        is_sibling = false;
    }
    Recomb_Node(const MAT::Node* node, int np, int p, char s)
        : node(node), node_parsimony(np), parsimony(p), is_sibling(s) {}
    inline bool operator<(const Recomb_Node &n) const {
        return (
                   ((*this).parsimony < n.parsimony) ||
                   ((this->node->identifier < n.node->identifier) && ((*this).parsimony == n.parsimony)));
    }
};

struct Recomb_Interval {
    Recomb_Node d; // donor
    Recomb_Node a; // acceptor
    int start_range_low;
    int start_range_high;
    int end_range_low;
    int end_range_high;
    Recomb_Interval(Recomb_Node donor, Recomb_Node acceptor, int sl, int sh,
                    int el, int eh)
        : d(donor), a(acceptor), start_range_low(sl), start_range_high(sh),
          end_range_low(el), end_range_high(eh) {}
    bool
    operator<(const Recomb_Interval &other) const { // compare second interval
        return end_range_low < other.end_range_low;
    }
};

struct Comp_First_Interval {
    inline bool operator()(const Recomb_Interval &one,
                           const Recomb_Interval &other) {
        return one.start_range_low < other.start_range_low;
    }
};
std::vector<Recomb_Interval>
combine_intervals(std::vector<Recomb_Interval> pair_list);
po::variables_map check_options(int argc, char **argv);
void ripples_mapper(const Pruned_Sample &sample,
                    Ripples_Mapper_Output_Interface &out,
                    size_t node_size,
                    const std::vector<int>& idx_map,
                    const std::vector<bool>& do_parallel,
                    const std::vector<Mapper_Info> &traversal_track,
                    const unsigned short tree_height,
                    const MAT::Node *root,
                    const MAT::Node *skip_node) ;
void ripplrs_merger(const Pruned_Sample &pruned_sample,
                    const std::vector<int> & idx_map,
                    const std::vector<MAT::Node *> &nodes_to_search,
                    size_t node_size, int pasimony_threshold,
                    const MAT::Tree &T,
                    tbb::concurrent_vector<Recomb_Interval> &valid_pairs,
                    const Ripples_Mapper_Output_Interface &out_ifc,
                    int nthreads, int branch_len, int min_range,
                    int max_range) ;
size_t check_parallelizable(const MAT::Node *root,
                            std::vector<bool> &do_parallel,
                            size_t parallel_threshold,
                            size_t check_threshold,
                            unsigned short& tree_height,
                            std::vector<Mapper_Info>& traversal_track,unsigned short level);