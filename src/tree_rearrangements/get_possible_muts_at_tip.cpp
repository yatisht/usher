#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <bits/types/FILE.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <sys/mman.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>
#undef NDEBUG
#include <cassert>
//#define NOT_FOUND 3
#define SURE 2
#define POSSIBLE 1
#define IMPOSSIBLE 0
#define INCREMENT 0x1000000

template<typename content_type>
class Temp_File_Map{
    char* filename;
    int fd;
    content_type* ptr;
    size_t mapped_size;
    size_t maxIdx;
    size_t increment;
    void allocate(){
        ftruncate(fd, mapped_size);
        ptr=(content_type*)mmap(0, mapped_size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (ptr==(content_type*)-1) {
            perror("");
            assert(false);
        }
    }
    public:
    Temp_File_Map(size_t size){
        mapped_size=size*sizeof(content_type);
        maxIdx=size;
        increment=size;
        filename=(char*)malloc(15);
        strcpy(filename, "aXXXXXXXXX");
        fd=mkstemp(filename);
        allocate();
    }
    ~Temp_File_Map(){
        munmap(ptr, mapped_size);
        close(fd);
        unlink(filename);
    }
    content_type& operator[](size_t idx){
        if (idx>=maxIdx) {
            resize(maxIdx+increment);
        }
        return ptr[idx];
    }
    void resize(size_t size){
        munmap(ptr, mapped_size);
        allocate();
    }
};

#ifdef TEST
Original_State_t original_state;
class compressed_states {
    boost::dynamic_bitset<> bits;

  public:
    compressed_states(size_t size) : bits(2 * size) {}
    char operator[](size_t idx) {
        return (bits[2 * idx + 1] << 1) | bits[2 * idx];
    }
    void set(size_t idx, char val) {
        assert(val<=3);
        bits[2 * idx] = val & 1;
        bits[2 * idx + 1] = val >> 1;
        assert((*this)[idx]==val);
    }
};
std::unordered_map<int, compressed_states> states;
#endif

static void
fill_state_vector(const Fitch_Sankoff::Scores_Type &score, compressed_states &out,
                  const MAT::Mutation& mutation,
                  const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    char ref_state=one_hot_to_two_bit(mutation.ref_nuc);
    size_t offset = dfs_ordered_nodes.size() - 1;
    out.set(0,
    Fitch_Sankoff::get_child_score_on_par_nuc(ref_state, score.back()).second);
    for (size_t idx = 1; idx < dfs_ordered_nodes.size(); idx++) {
        if (dfs_ordered_nodes[idx]->is_leaf()) {
            auto node_iter=original_state.find(dfs_ordered_nodes[idx]->identifier);
            char nuc=ref_state;
            if (node_iter!=original_state.end()) {
                auto mut_iter=node_iter->second.find(mutation);
                if (mut_iter!=node_iter->second.end()) {
                    nuc=one_hot_to_two_bit(mut_iter->mut_nuc);
                }
            }
            out.set(idx, nuc);
            if(dfs_ordered_nodes[idx]->parent->children.size()==1){
                out.set(dfs_ordered_nodes[idx]->parent->index, nuc);
            }
            
        }else{
        char parent_state = out[dfs_ordered_nodes[idx]->parent->index];
        const auto& this_score=score[offset-idx];
        //char fs_nuc=Fitch_Sankoff::get_child_score_on_par_nuc(out[parent_idx],score[offset-idx]).second;
        char nuc=std::min_element(this_score.begin(),this_score.end())-this_score.begin();
        if (this_score[nuc]==this_score[parent_state]) {
            nuc=parent_state;
        }
        assert(((nuc==parent_state)?this_score[nuc]:(this_score[nuc]+1))==Fitch_Sankoff::get_child_score_on_par_nuc(parent_state,score[offset-idx]).first);
        out.set(idx, nuc);
        }
    }
}
struct Mutation_State {
    const MAT::Mutation *mutation;
    char state;
    Mutation_State(const MAT::Mutation &mutation)
        : mutation(&mutation), state(SURE) {}
    Mutation_State(const Mutation_State &other, char state)
        : mutation(other.mutation), state(state) {}
    int get_pos() const { return mutation->position; }
    char get_mut() const {
        assert(mutation->mut_nuc);
        return mutation->mut_nuc;
    }
    char get_state() const { return state; }
    bool operator<(const Mutation_State &other) const {
        return (get_pos() < other.get_pos()) ||
               ((get_pos() == other.get_pos()) &&
                (get_mut() < other.get_mut()));
    }
};
typedef std::vector<Mutation_State> Tip_Mutation;

static bool add_mutation(Tip_Mutation& add_to,const Mutation_State& to_add,char state){
    if ((!add_to.empty())&&(add_to.back().get_pos()==to_add.get_pos())) {
        if (state==IMPOSSIBLE) {
            assert(add_to.back().get_mut()!=to_add.get_mut());
            return false;
        }else if (add_to.back().get_state()==IMPOSSIBLE) {
            add_to.back()=to_add;
            return true;
        }
    }
    add_to.emplace_back(to_add,state);
    return true;
}

template<typename content_type,typename individual_type>
class Range {
    typename content_type::const_iterator start;
    typename content_type::const_iterator end;
  public:
    Range(const content_type &in) : start(in.begin()), end(in.end()) {}
    Range(typename content_type::const_iterator& start,typename content_type::const_iterator& end) : start(start), end(end) {}
    typename content_type::const_iterator &operator->() { return start; }
    individual_type &operator*() { return *start; }
    const typename content_type::const_iterator &operator->() const { return start; }
    const individual_type &operator*() const { return *start; }
    void operator++(int) { start++; }
    operator bool() { return start != end; }
};
    
struct Merged_Mutation_Comparator {
    bool operator()(const Range<Tip_Mutation,Mutation_State> &lhs, const Range<Tip_Mutation,Mutation_State> rhs) {
        return lhs->get_pos() > rhs->get_pos() ||
               ((lhs->get_pos() == rhs->get_pos()) &&
                (lhs->get_mut() > rhs->get_mut()));
    }
};
struct Nothing{};
class Merged_Mutation_Initer{
    Tip_Mutation* tip_mutations;
    public:
    Merged_Mutation_Initer(Tip_Mutation* start):tip_mutations(start){}
    void operator()(const MAT::Node *this_node,  std::vector<Range<Tip_Mutation,Mutation_State>> & heap,Nothing& ignored)const {
                heap.reserve(this_node->children.size());
        for (MAT::Node *child : this_node->children) {
            if (!tip_mutations[child->index].empty()) {
                heap.emplace_back(tip_mutations[child->index]);
            }
        }
    }
    
};
template<typename content_type,typename individual_type,typename comparator,typename Initer,typename Other>
class Merged_Iterator {
#ifdef TEST
    /*int lastpos;
    char lastmut;
    */
#endif
    typedef Range<content_type,individual_type> Range_Type;
    std::vector<Range_Type> heap;
    Other other;
  public:
    Merged_Iterator(const Initer& initer,
                             const MAT::Node *this_node) {

        initer(this_node,heap,other);
        std::make_heap(heap.begin(), heap.end(), comparator());
/*#ifdef TEST
        if (!heap.empty()) {
            lastmut = heap.front()->get_mut();
            lastpos = heap.front()->get_pos();
        }
#endif*/
    }

    void operator++(int) {
        std::pop_heap(heap.begin(), heap.end(), comparator());
        Range_Type &incremented_iter = heap.back();
        incremented_iter++;
        if (incremented_iter) {
            std::push_heap(heap.begin(), heap.end(), comparator());
        } else {
            heap.pop_back();
        }
/*#ifdef TEST
        if (!heap.empty()) {
            assert(heap.front()->get_pos() > lastpos ||
                   ((heap.front()->get_pos() == lastpos) &&
                    heap.front()->get_mut() >= lastmut));
            lastmut = heap.front()->get_mut();
            lastpos = heap.front()->get_pos();
        }
#endif*/
    }
    const individual_type &operator*() const { return *heap.front(); }

    operator bool() { return !heap.empty(); }
};

typedef Merged_Iterator<Tip_Mutation, Mutation_State, Merged_Mutation_Comparator,Merged_Mutation_Initer,Nothing> Merged_Mutation_Iterator;

static void add_mut(Tip_Mutation &out,
               Mutation_State &last_mutation, size_t &total_count,
               size_t &threshold, bool &is_even) {
    if (total_count > threshold) {
        add_mutation(out, last_mutation, SURE);
    } else if (is_even && total_count == threshold) {
        add_mutation(out, last_mutation, POSSIBLE);
    } else {
        add_mutation(out, last_mutation, IMPOSSIBLE);
    }
}
static Tip_Mutation *compute_possible_mutations_at_tip(
    const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    Tip_Mutation *result = new Tip_Mutation[dfs_ordered_nodes.size()];
    Merged_Mutation_Initer init(result);
    for (size_t node_idx = dfs_ordered_nodes.size() - 1;; node_idx--) {
        MAT::Node *this_node = dfs_ordered_nodes[node_idx];
        if (this_node->is_leaf()) {
            const auto& mutation_set=original_state[this_node->identifier];
            result[node_idx].reserve(mutation_set.size());
            for (const MAT::Mutation &m : mutation_set) {
                result[node_idx].emplace_back(m);
            }
            std::sort(result[node_idx].begin(),result[node_idx].end());
        } else {
            Merged_Mutation_Iterator merged_iter(init, this_node);
            if (!merged_iter) {
                continue;
            }
            auto last_mutation = *merged_iter;
            size_t total_count = 0;
            size_t sure_count = 0;
            size_t threshold = this_node->children.size() / 2;
            bool is_even = !(this_node->children.size() & 1);
            while (merged_iter) {
                auto this_mutation = *merged_iter;
                if (this_mutation.get_pos() != last_mutation.get_pos() ||
                    this_mutation.get_mut() != last_mutation.get_mut()) {
                    add_mut(result[node_idx], last_mutation,
                              total_count, threshold, is_even);
                    sure_count = 0;
                    total_count = 0;
                    last_mutation = this_mutation;
                }

                switch (this_mutation.state) {
                case SURE:
                    sure_count++;
                case POSSIBLE:
                    total_count++;
                }

                merged_iter++;
            }
            add_mut(result[node_idx], last_mutation,total_count, threshold, is_even);
#ifdef TEST
            std::pair<size_t, size_t> range = Fitch_Sankoff::dfs_range(
                dfs_ordered_nodes[node_idx], dfs_ordered_nodes);
            for (auto &m : result[node_idx]) {
                if (m.get_state() == POSSIBLE) {
                    continue;
                }
                Fitch_Sankoff::Scores_Type scores;
                Fitch_Sankoff::sankoff_backward_pass(range, dfs_ordered_nodes,
                                                     scores, original_state,
                                                     *(m.mutation), 1);
                auto &tip_score = scores.back();
                if (m.get_state() == SURE) {
                    char nuc = one_hot_to_two_bit(m.get_mut());
                    for (char nuc_idx = 0; nuc_idx < 4; nuc_idx++) {
                        if (nuc_idx != nuc)
                            assert(tip_score[nuc] < tip_score[nuc_idx]);
                    }
                } else {
                    char nuc = one_hot_to_two_bit(m.mutation->ref_nuc);
                    for (char nuc_idx = 0; nuc_idx < 4; nuc_idx++) {
                        assert(tip_score[nuc] <= tip_score[nuc_idx]);
                    }
                }
            }
#endif
        }
        if (node_idx == 0) {
            break;
        }
    }
    return result;
}

class parent_this_iter {
    Tip_Mutation::const_iterator this_iter;
    Tip_Mutation::const_iterator this_iter_end;
    Tip_Mutation::const_iterator parent_iter;
    Tip_Mutation::const_iterator parent_iter_end;

  public:
    parent_this_iter(const Tip_Mutation &this_states,
                     const Tip_Mutation &parent_states)
        : this_iter(this_states.begin()), this_iter_end(this_states.end()),
          parent_iter(parent_states.begin()),
          parent_iter_end(parent_states.end()) {}
    const Mutation_State &next(char &this_state, char &parent_state) {
#define inc_this \
const Mutation_State &to_return = *this_iter; \
this_state = to_return.state;\
parent_state = IMPOSSIBLE;\
this_iter++;\
return to_return;

#define inc_parent \
const Mutation_State &to_return = *parent_iter;\
this_state = IMPOSSIBLE;\
parent_state = parent_iter->state;\
parent_iter++;\
return to_return;

        if (parent_iter == parent_iter_end) {
            inc_this
        } else if (this_iter == this_iter_end) {
            inc_parent
        } else if (*this_iter < *parent_iter) {
            inc_this
        } else if (this_iter->get_pos() == parent_iter->get_pos() &&
                   this_iter->get_mut() == parent_iter->get_mut()) {
            const Mutation_State &to_return = *this_iter;
            this_state = to_return.state;
            parent_state = parent_iter->state;
            this_iter++;
            parent_iter++;
            return to_return;
        }
        inc_parent
    }
    operator bool() {
        return this_iter < this_iter_end || parent_iter < parent_iter_end;
    }
};
#ifdef TEST
static void
check_parent_state(const Tip_Mutation *start, const MAT::Node *this_node,
                   const std::vector<MAT::Node *> &dfs_ordered_nodes,
                   const Tip_Mutation &state) {
    boost::dynamic_bitset<unsigned long, std::allocator<unsigned long>> visited(
        state.size());
    const MAT::Node *parent_node = this_node;
    while (parent_node) {
        auto parent_node_state_iter = start[parent_node->index].begin();
        auto parent_node_state_iter_end = start[parent_node->index].end();
        for (size_t idx = 0; idx < state.size(); idx++) {
            if (parent_node_state_iter == parent_node_state_iter_end) {
                break;
            }
            auto parent_pos=parent_node_state_iter->get_pos();
            auto this_pos=state[idx].get_pos();
            if (parent_pos == this_pos&&
                parent_node_state_iter->get_mut() == state[idx].get_mut()) {
                if (!visited[idx]) {
                    assert(parent_node_state_iter->get_state() ==
                           state[idx].get_state());
                    visited[idx] = true;
                }
                parent_node_state_iter++;
            } else {
                if(parent_pos==this_pos){
                    //while ((parent_node_state_iter < parent_node_state_iter_end)&&(parent_node_state_iter->get_mut()>state[idx].get_mut()&&parent_node_state_iter->get_pos()==state[idx].get_pos())) {
                    assert(parent_node_state_iter->get_state()==IMPOSSIBLE);
                    parent_node_state_iter++;
                    //}
                }else{
                    assert(parent_node_state_iter->get_pos()>this_pos);
                }
                continue;
            }
        }
        parent_node = parent_node->parent;
    }
    for (size_t idx = 0; idx < state.size(); idx++) {
        if (!visited[idx]&&(state[idx].get_state()!=IMPOSSIBLE)) {
            fprintf(stderr, "Mutation at %d of state %d not found \n",
                    state[idx].get_pos(), state[idx].get_state());
            assert(false);
        }
    }
}

#endif
static void
cleanup_tip_mutations_helper(const std::vector<MAT::Node *> &dfs_ordered_nodes,
                             Tip_Mutation *start, size_t index,
                             const Tip_Mutation &parent_state) {
    parent_this_iter iter(start[index], parent_state);
    Tip_Mutation this_node_state;
    Tip_Mutation new_this_node_change;
    char parent_found;
    char this_found;
    while (iter) {
        const Mutation_State &this_mutation =
            iter.next(this_found, parent_found);

        switch (parent_found) {
        //Both impossible and not found mean this mutation is not at the tip
        //case NOT_FOUND:
        case IMPOSSIBLE:
            if (this_found == SURE) {
                add_mutation(new_this_node_change,this_mutation,SURE);
                add_mutation(this_node_state,this_mutation,SURE);
            } else {
                // when this state is POSSIBLE (same change in parsimony score
                // of child node regardless of the state), better to not have
                // this change to save a parsimony score of 1 at this level
                // this node state is the same as parent (also for the state of
                // POSSIBLE, NOT_FOUND,IMPOSSIBLE), so not push
                add_mutation(this_node_state,this_mutation,IMPOSSIBLE);
            }
            break;
        case POSSIBLE:
        //this gives no additional information on subtree
            add_mutation(this_node_state,this_mutation,this_found);
            add_mutation(new_this_node_change,this_mutation,this_found);
            break;
        case SURE:
            if (this_found == IMPOSSIBLE) {
                add_mutation(new_this_node_change,this_mutation,IMPOSSIBLE);
                add_mutation(this_node_state,this_mutation,IMPOSSIBLE);
            } else {
                add_mutation(this_node_state,this_mutation,SURE);
            }
            break;
        default:
            assert(false);
        }
    }
    start[index].swap(new_this_node_change);
#ifdef TEST
    check_parent_state(start, dfs_ordered_nodes[index], dfs_ordered_nodes,
                       this_node_state);
    for (auto &m : this_node_state) {
        if (m.get_state() == POSSIBLE) {
            continue;
        }
        char ref_nuc = m.mutation->ref_nuc;
        auto iter = states.find(m.get_pos());
        if (iter == states.end()) {
            Fitch_Sankoff::Scores_Type scores;
            Fitch_Sankoff::sankoff_backward_pass(
                std::make_pair(0, dfs_ordered_nodes.size()), dfs_ordered_nodes,
                scores, original_state, *m.mutation, ref_nuc);
            iter = states.emplace(m.get_pos(), dfs_ordered_nodes.size()).first;
            fill_state_vector(scores, iter->second, *m.mutation, dfs_ordered_nodes);
        }
        char actual_nuc = 1 << (iter->second[index]);
        if (m.get_state() == IMPOSSIBLE) {
            assert(actual_nuc == ref_nuc);
        } else {
            assert(m.get_mut() == actual_nuc);
        }
    }
#endif

    for (MAT::Node *child : dfs_ordered_nodes[index]->children) {
        cleanup_tip_mutations_helper(dfs_ordered_nodes, start, child->index,
                                     this_node_state);
    }
}


class Position{
    int content;
    public:
    operator int() const{
        return content&0x7fffffff; 
    }
    explicit operator bool() const{
        return content&0x80000000;
    }
    Position(int pos,bool repeated){
        content=pos|(repeated<<31);
    }
};
class Repeated_Mutations{
    size_t start_offset;
    size_t size;
    public:
    struct const_iterator{
        union{
        size_t offset;
        std::vector<Position>::const_iterator vec_iter;
        };
        bool is_vec;
        void operator++(){
            offset++;
        }
        const Position& operator*() const{
            if(is_vec){
                return *vec_iter;
            }
            return Repeated_Mutations::memory[offset];
        }
        bool operator==(const const_iterator& other){
            return offset==other.offset;
        }
        const_iterator(size_t offset):offset(offset),is_vec(false){}
        const_iterator(const std::vector<Position>::const_iterator& vec_iter):vec_iter(vec_iter),is_vec(true){}
        
    };
    static size_t allocatable_offset;
    static Temp_File_Map<Position>& memory;
    //leaf node
    Repeated_Mutations():start_offset(allocatable_offset),size(0){}
    const_iterator begin()const{
        return const_iterator(start_offset);
    }
    const_iterator end()const{
        return const_iterator(start_offset+size);
    }
    void add_back(int pos,bool repeated){
        assert(start_offset+size==allocatable_offset);
        memory[allocatable_offset]=Position(pos,repeated);
        size++;
        allocatable_offset++;
    }
};

struct Repeated_Mutation_Comparator{
    bool operator()(const Range<Repeated_Mutations,Position>& lhs,const Range<Repeated_Mutations,Position>& rhs){
        return *lhs>*rhs;
    }
};
typedef std::vector<std::vector<Position>> Repeated_Mutations_Temp_Storage;
struct Repeated_Mutation_Initer{
    static std::vector<Repeated_Mutations>& all;
    void operator()(const MAT::Node *this_node,  std::vector<Range<Repeated_Mutations,Position>> & heap,Repeated_Mutations_Temp_Storage& temp){
        temp.emplace_back();
        auto& this_vect=temp.back();
        for(const auto& m:this_node->mutations){
            this_vect.emplace_back(m.position,false);
        }
        Repeated_Mutations::const_iterator start (this_vect.begin());
        Repeated_Mutations::const_iterator end (this_vect.end());
        heap.push_back(Range<Repeated_Mutations,Position>(start,end));
        for(const auto child:this_node->children){
            const auto& temp=all[child->index];
            heap.emplace_back(temp.begin(),temp.end());
        }
    }
};
typedef Merged_Iterator<Repeated_Mutations, int, Repeated_Mutation_Comparator,Repeated_Mutation_Initer,Repeated_Mutations_Temp_Storage> Repeated_Mutation_Iterator;



void get_all_mutation(std::vector<MAT::Node*>& dfs_ordered_nodes,std::vector<int> all_positions){
    Repeated_Mutations::memory=Temp_File_Map<Position>(INCREMENT*sizeof(Position));
    Temp_File_Map<Repeated_Mutations> repeated_mutations(dfs_ordered_nodes.size());
    
}
int main(int argc, char **argv) {
    MAT::Tree tree = MAT::load_mutation_annotated_tree(argv[1]);
    check_samples(tree.root, original_state, &tree);
    std::vector<MAT::Node *> dfs_ordered_nodes = tree.depth_first_expansion();
    Tip_Mutation *result = compute_possible_mutations_at_tip(dfs_ordered_nodes);
    cleanup_tip_mutations_helper(dfs_ordered_nodes, result, 0, Tip_Mutation());
}