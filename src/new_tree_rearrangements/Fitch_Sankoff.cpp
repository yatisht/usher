#include "Fitch_Sankoff.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <smmintrin.h>
#include <string>
#include <array>
#include <unordered_map>
#include <array>
#include <vector>
#include <emmintrin.h>
namespace MAT = Mutation_Annotated_Tree;
//get state of ancestor at position
nuc_one_hot get_this_state(MAT::Node* ancestor,int position){
            auto iter = ancestor->mutations.find(position);
    if (iter == ancestor->mutations.end()) {
                ancestor = ancestor->parent;
            } else {
                return iter->get_mut_one_hot();
            }
    while (ancestor) {
            auto iter = ancestor->mutations.find(position);
            if (iter == ancestor->mutations.end()) {
                ancestor = ancestor->parent;
            } else {
                return iter->get_mut_one_hot();
            }
        }
    return MAT::Mutation::refs[position];
}

static void set_state_2(uint8_t child1,uint8_t child2,uint8_t& boundary1_major_allele){
    uint8_t major_abs=child1&child2&0xf;
    uint8_t nuc_present=(child1|child2)&0xf;
    if (major_abs) {
        uint8_t minor_allele=nuc_present&(~major_abs);
        boundary1_major_allele=(minor_allele<<4)|major_abs;
        //Not seting minor allele if nuc_present is ambiguous, because it cannot be boundary
    }else {
        //This is a tie, Forward pass will break tie with parent, (or arbitrarily if parent agree with neither), taking it as tie case, and boundary case, as alleles not present have a difference of 1 with present allele (can become a 3-allele tie and then follow parent when moving nodes around)
        boundary1_major_allele=((~nuc_present)<<4)|nuc_present;
    }
}
//move the least significant bit of each nibble to the first 4 bits to form one hot state encoding
static uint8_t movemask(int in){
    return (1&in)|(2&(in>>3))|(4&(in>>6))|(8&(in>>9));
}

void set_state_from_cnt(const std::array<int,4>& data, uint8_t& boundary1_major_allele_out){
    //load data
    __m128i ori=_mm_loadu_si128((__m128i*)data.data());
    //shufle it to 1,0,3,2 order
    __m128i suf=_mm_shuffle_epi32(ori, 0x4e);
    //compare, max1 contain max(3,1),max(2,0),max(3,1),max(2,0)
    __m128i max1=_mm_max_epi32(ori, suf);
    //suffle again, get max(2,0),max(3,1),max(2,0),max(3,1)
    __m128i max1_suf=_mm_shuffle_epi32(max1, 0x11);
    //now max_values contain max of 4 on all lanes
    __m128i max_values=_mm_max_epi32(max1, max1_suf);
    __m128i one=_mm_set1_epi32(1);
    //find which allele count is the same as max allele count, which will be set to all 1 on that lane, then extract them by taking most significant bit of each lane
    int is_max_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(ori, max_values));
    //same for boundary one allele, just find which one equals to max count -1
    int is_boundary_1_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_add_epi32(ori,one), max_values));
    uint8_t max_mask=movemask(is_max_raw);
    uint8_t boundary1_mask=movemask(is_boundary_1_raw);
    boundary1_major_allele_out=max_mask|(boundary1_mask<<4);
}

//following are different specialization for get allele count for different number of children
//all of them are shifting a mask for each allele, and then use popcnt to get count
static void get_state_count8(size_t start_idx,size_t size, std::vector<uint8_t>& minor_major_allele,std::array<int,4>& nuc_count){
    uint64_t child_states= *((uint64_t*)(&(minor_major_allele[start_idx])));
    uint64_t mask=0x0101010101010101;
    child_states<<=(8*(8+start_idx-size));
    for (size_t nuc_idx=0; nuc_idx<4; nuc_idx++) {
        nuc_count[nuc_idx]+=__builtin_popcountl(mask&(child_states>>nuc_idx));
    }
}

static void get_state_count16(size_t start_idx,size_t size, std::vector<uint8_t>& minor_major_allele,std::array<int,4>& nuc_count){
    uint64_t child_states= *((uint64_t*)(&(minor_major_allele[start_idx])));
    //as major allele and boundary1 allele are interleaved, possible to do more in a batch by deinterleaving them
    child_states&=0x0f0f0f0f0f0f0f0f;
    uint64_t child_states_after=(0x0f0f0f0f0f0f0f0f&(*((uint64_t*)(&(minor_major_allele[start_idx+8])))));
    child_states_after<<=((8*(16+start_idx-size))+4);
    child_states|=child_states_after;

    uint64_t mask=0x0101010101010101;
    for (size_t nuc_idx=0; nuc_idx<8; nuc_idx++) {
        nuc_count[nuc_idx&3]+=__builtin_popcountl(mask&(child_states>>nuc_idx));
    }
}

static bool get_state_count32(size_t start_idx,size_t size, std::vector<uint8_t>& minor_major_allele,std::array<int,4>& nuc_count){
    bool is_finished=false;
    __m128i child_states_after=_mm_loadu_si128((__m128i*)(&(minor_major_allele[start_idx+16])));
    __m128i child_state_mask=_mm_set1_epi64x(0x0f0f0f0f0f0f0f0fl);
    __m128i mask=_mm_set1_epi64x(0x8080808080808080l);
    __m128i child_states=_mm_loadu_si128((__m128i*)(&(minor_major_allele[start_idx])));
    size_t to_calculate=size-start_idx;
    unsigned long unloadmask1=-1;
    unsigned long unloadmask2=-1;
    if (to_calculate<32) {
        if (to_calculate>24) {
            unloadmask1>>=(8*(32-to_calculate));
        }else {
            #ifdef DETAIL_DEBUG_FITCH_SANKOFF
            assert(to_calculate>16);
            assert(to_calculate<=24);
            #endif
            unloadmask1=0;
            unloadmask2>>=(8*(24-to_calculate));
        }
        is_finished=true;
    }else if (to_calculate==32) {
        is_finished=true;
    }
    __m128i unload_mask=_mm_set_epi64x(unloadmask1, unloadmask2);
    child_states=_mm_and_si128(child_states, child_state_mask);
    child_states_after=_mm_and_si128(child_states_after, child_state_mask);
    child_states_after=_mm_slli_epi64(child_states_after,4 );
    child_states_after=_mm_and_si128(child_states_after, unload_mask);
    child_states=_mm_or_si128(child_states, child_states_after);

    for (int nuc_idx=7; nuc_idx>=0; nuc_idx--) {
        nuc_count[nuc_idx&3]+=__builtin_popcountl(_mm_movemask_epi8(_mm_and_si128(mask, child_states)));
        child_states=_mm_slli_epi32(child_states, 1);
    }

    return is_finished;
}

#ifdef DETAIL_DEBUG_FITCH_SANKOFF
std::array<int,4> count_right(MAT::Node* this_node, std::vector<uint8_t>& minor_major_allele){
    std::array<int,4> to_return{0,0,0,0};
    for(auto child:this_node->children){
        for (char nu_idx=0; nu_idx<4; nu_idx++) {
            if (minor_major_allele[child->bfs_index]&(1<<nu_idx)) {
                to_return[nu_idx]++;
            }
        }
    }
    return to_return;
}
#endif

void FS_backward_pass(const std::vector<backward_pass_range>& child_idx_range, std::vector<uint8_t>& boundary1_major_allele,const std::unordered_map<std::string, nuc_one_hot>& mutated,nuc_one_hot ref_nuc){
    //Using BFS order for memory locality, as children of a node are toghrther in BFS order
    for(long node_idx=child_idx_range.size()-1;node_idx>=0;node_idx--){
        auto child_size=child_idx_range[node_idx].child_size;

        //leaf node
        if (child_size==0) {
            auto iter=mutated.find(*child_idx_range[node_idx].identifier);
            nuc_one_hot allele;
            if (iter==mutated.end()) {
                allele=ref_nuc;
            }else {
                allele=iter->second;
            }
            boundary1_major_allele[node_idx]=allele;
        }else if (child_size==1) {
            auto child_idx=child_idx_range[node_idx].first_child_bfs_idx;
            boundary1_major_allele[node_idx]=boundary1_major_allele[child_idx]&0xf;
        }else if (child_size==2) {
            size_t child_start_idx=child_idx_range[node_idx].first_child_bfs_idx;
            set_state_2(boundary1_major_allele[child_start_idx],boundary1_major_allele[child_start_idx+1],boundary1_major_allele[node_idx]);
        }else {
            size_t child_start_idx=child_idx_range[node_idx].first_child_bfs_idx;
            size_t child_end_idx=child_start_idx+child_size;
           //dispatch on children size
            std::array<int,4> nuc_count{0,0,0,0};
            while (true) {
                auto child_left=child_end_idx-child_start_idx;
                if (child_left<=8) {
                    get_state_count8(child_start_idx,child_end_idx,boundary1_major_allele,nuc_count);
                    break;
                }else if (child_left<=16) {
                    get_state_count16(child_start_idx,child_end_idx,boundary1_major_allele,nuc_count);
                    break;
                }else {
                    if (get_state_count32(child_start_idx,child_end_idx,boundary1_major_allele,nuc_count)) {
                        break;
                    }else {
                        child_start_idx+=32;
                    }
                }
            }
            #ifdef DETAIL_DEBUG_FITCH_SANKOFF
            auto right_count=count_right(this_node, boundary1_major_allele);
            assert(nuc_count==right_count);
            #endif
            set_state_from_cnt(nuc_count, boundary1_major_allele[node_idx]);
        }
        assert(boundary1_major_allele[node_idx]&0xf);
    }
}
//add mutation for non-binary nodes
static nuc_one_hot set_state(const forward_pass_range & this_range,uint8_t boundary1_major_allele,nuc_one_hot par_state,const MAT::Mutation& base,tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>& output,MAT::Tree* try_similar){
    nuc_one_hot this_state;
    bool need_add=false;
    //follow parent if it can
    if(boundary1_major_allele&par_state){
        this_state=par_state;
    }else {
        this_state=boundary1_major_allele&0xf;
        this_state=this_state.choose_first();
#ifdef CHECK_STATE_REASSIGN
        if(try_similar){
            nuc_one_hot ori_state=get_this_state(try_similar->get_node(this_node->identifier), base.get_position());
            if (ori_state&this_state) {
                this_state=ori_state;
            }
        }
#endif
        //assert(this_node->parent->children.size()>1);
        need_add=true;
    }

    nuc_one_hot major_allele(boundary1_major_allele&0xf);

    nuc_one_hot boundary1_allele=boundary1_major_allele>>4;
    assert(major_allele&this_state);
    //add if have allele whose count is exactly one less than major allele count,unless it only have one children to save memory, or major allele count is not the same as parent allele
    if (major_allele.is_ambiguous()||(boundary1_allele)||need_add) {
        MAT::Mutation to_add(base);
        to_add.set_par_mut(par_state, this_state);
        if (this_range.child_size<=1) {
            assert(!boundary1_allele);
            boundary1_allele=(~major_allele)&0xf;
        }
        to_add.set_auxillary(major_allele,boundary1_allele);
        output.push_back(to_add);
    }

    return this_state;
}

static nuc_one_hot set_binary_node_state(uint8_t this_boundary1_major_allele, nuc_one_hot par_state, nuc_one_hot left_child_major_allele,nuc_one_hot right_child_major_allele,const MAT::Mutation& base,tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>& output,MAT::Tree* try_similar){
    nuc_one_hot this_major_allele=this_boundary1_major_allele&0xf;
    if (this_major_allele&par_state) {
        if (this_major_allele!=par_state||(this_boundary1_major_allele>>4)) {
            MAT::Mutation to_add(base);
            to_add.set_par_mut(par_state, par_state);
            //same as the generic node, but also record the state of its two children
            to_add.set_children(this_boundary1_major_allele>>4,this_major_allele, left_child_major_allele, right_child_major_allele);
            output.push_back(to_add);
        }
        return par_state;
    }
    nuc_one_hot this_state=this_major_allele.choose_first();
#ifdef CHECK_STATE_REASSIGN
    if (try_similar) {
    nuc_one_hot ori_state=get_this_state(try_similar->get_node(node->identifier), base.get_position());
    if (ori_state&this_major_allele) {
        this_state=ori_state;
    }
    }
#endif
    MAT::Mutation to_add(base);
    to_add.set_par_mut(par_state, this_state);
    to_add.set_children(this_boundary1_major_allele>>4,this_major_allele, left_child_major_allele, right_child_major_allele);
    output.push_back(to_add);
    return this_state;
}
//just dispatch based on number of children
void forward_per_node(
    const forward_pass_range & this_range,
    const std::vector<uint8_t> &boundary1_major_allele,
    const MAT::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,
    std::vector<nuc_one_hot> &states, size_t node_idx,
    nuc_one_hot parent_state,MAT::Tree* try_similar) {
    if (this_range.child_size == 2) {
        states[node_idx] = set_binary_node_state(
            boundary1_major_allele[node_idx], parent_state,
            boundary1_major_allele[this_range.left_child_idx] & 0xf,
            boundary1_major_allele[this_range.right_child_idx] & 0xf, base,
            output[node_idx],try_similar);
    } else {
        states[node_idx] =
            set_state(this_range, boundary1_major_allele[node_idx], parent_state,
                      base, output[node_idx],try_similar);
    }
}
static void FS_forward_pass(
    const std::vector<forward_pass_range>& forward_pass_idx,
    const std::vector<uint8_t> &boundary1_major_allele,
    const MAT::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,MAT::Tree* try_similar) {
    std::vector<nuc_one_hot> states(forward_pass_idx.size());
    //set state for root node
    forward_per_node(forward_pass_idx[0], boundary1_major_allele, base, output,
                  states, 0,base.get_ref_one_hot(),try_similar);
    //parent is visited before children in bfs order, so iterate in forward order
    for (size_t node_idx=1; node_idx<forward_pass_idx.size(); node_idx++) {
        nuc_one_hot parent_state= states[forward_pass_idx[node_idx].parent_bfs_idx];
        forward_per_node(forward_pass_idx[node_idx], boundary1_major_allele, base, output,
                  states, node_idx, parent_state,try_similar);
    }
}
#if defined CHECK_STATE_REASSIGN || defined DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
static nuc_one_hot get_state_only(uint8_t boundary1_majority_allele,nuc_one_hot parent_allele){
    nuc_one_hot majority_allele=boundary1_majority_allele&0xf;
    nuc_one_hot this_state=majority_allele&parent_allele;
    this_state=this_state?this_state:majority_allele;
    return this_state.choose_first();
}
int FS_forward_assign_states_only(const std::vector<MAT::Node*>& bfs_ordered_nodes,const std::vector<uint8_t>& boundary1_major_allele,const nuc_one_hot parent_state,std::vector<uint8_t>& states_out,std::vector<std::vector<MAT::Node*>>& children_mutation_count){
    states_out[0]=get_state_only(boundary1_major_allele[0], parent_state);
    int mutation_count;
    if ((states_out[0]!=parent_state)) {
        mutation_count=1;
    }else {
        mutation_count=0;
    }
    for (size_t i=1; i<bfs_ordered_nodes.size(); i++) {
        size_t parent_idx=bfs_ordered_nodes[i]->parent->bfs_index;
        nuc_one_hot parent_state=states_out[parent_idx];
        states_out[i]=get_state_only(boundary1_major_allele[i],parent_state);
        assert(__builtin_popcount(states_out[i])==1);
        if (states_out[i]!=parent_state) {
            mutation_count++;
            children_mutation_count[parent_idx].push_back(bfs_ordered_nodes[i]);
        }
    }
    return mutation_count;
}
#endif
void Fitch_Sankoff_Whole_Tree(const std::vector<backward_pass_range>& child_idx_range,const std::vector<forward_pass_range>& parent_idx,const MAT::Mutation & base,const std::unordered_map<std::string, nuc_one_hot>& mutated,std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>& output,MAT::Tree* try_similar){
    std::vector<uint8_t> minor_major_allele(child_idx_range.size()+16);

    FS_backward_pass(child_idx_range,minor_major_allele,mutated,base.get_ref_one_hot());
    FS_forward_pass(parent_idx,minor_major_allele,base,output,try_similar);
}
void Fitch_Sankoff_prep(const std::vector<Mutation_Annotated_Tree::Node*>& bfs_ordered_nodes, std::vector<backward_pass_range>& child_idx_range,std::vector<forward_pass_range>& parent_idx){
    child_idx_range.reserve(bfs_ordered_nodes.size());
    parent_idx.reserve(bfs_ordered_nodes.size());
    for(const auto node:bfs_ordered_nodes){
        child_idx_range.emplace_back();
        child_idx_range.back().child_size=node->children.size();
        if (child_idx_range.back().child_size) {
            child_idx_range.back().first_child_bfs_idx=node->children[0]->bfs_index;
        }else {
            child_idx_range.back().identifier=&node->identifier;
        }
        parent_idx.emplace_back();
        if (node->parent) {
            parent_idx.back().parent_bfs_idx=node->parent->bfs_index;
        }
        parent_idx.back().child_size=node->children.size();
        if(parent_idx.back().child_size==2){
            parent_idx.back().left_child_idx=node->children[0]->bfs_index;
            parent_idx.back().right_child_idx=node->children[1]->bfs_index;
        }
    }
}
