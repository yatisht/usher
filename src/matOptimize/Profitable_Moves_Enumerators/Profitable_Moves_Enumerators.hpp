#ifndef PROFITABLE_MOVES_ENUMERATOR
#define PROFITABLE_MOVES_ENUMERATOR
#include "../tree_rearrangement_internal.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
#include "../stack_allocator.hpp"

namespace MAT = Mutation_Annotated_Tree;
//Class for recording change in major allele set
class Mutation_Count_Change {
    int position;
    uint8_t chromIdx;
    nuc_one_hot decremented_allele;
    nuc_one_hot incremented_allele;
    nuc_one_hot par_state;

  public:
    static const char VALID_MASK=1;
    static const char END_MASK=2;
    Mutation_Count_Change() {
        decremented_allele=0;
        incremented_allele=0;
    }
    //copy over position, original state, par state from a mutation
    Mutation_Count_Change(const Mutation_Count_Change& child_mut_count,nuc_one_hot new_major_allele) {
        position=child_mut_count.position;
        chromIdx=child_mut_count.chromIdx;
        auto ori_state=child_mut_count.par_state;
        par_state=child_mut_count.par_state;
        incremented_allele=new_major_allele&(~ori_state);
        decremented_allele=ori_state&(~new_major_allele);
    }
    Mutation_Count_Change(const MAT::Mutation &pos,nuc_one_hot decremented,nuc_one_hot incremented) {
        position = pos.get_position();
        chromIdx = pos.get_chromIdx();
        par_state=pos.get_par_one_hot();
        set_change(decremented, incremented);
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
#endif