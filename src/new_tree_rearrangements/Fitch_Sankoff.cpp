#include "Fitch_Sankoff.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <smmintrin.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <emmintrin.h>
namespace MAT = Mutation_Annotated_Tree;

static void set_state_2(uint8_t child1,uint8_t child2,uint8_t& boundary1_major_allele,uint8_t& boundary2_allele){
    uint8_t major_abs=child1&child2&0xf;
    uint8_t nuc_present=(child1|child2)&0xf;
    if (major_abs) {
        uint8_t minor_allele=nuc_present&(~major_abs);
        boundary1_major_allele=(minor_allele<<4)|major_abs;
        boundary2_allele=(~nuc_present)&0xf;
        //Not seting minor allele if nuc_present is ambiguous, because it cannot be boundary
    }else {
        //This is a tie, Forward pass will break tie with parent, (or arbitrarily if parent agree with neither), taking it as tie case, and boundary case, as alleles not present have a difference of 1 with present allele (can become a 3-allele tie and then follow parent when moving nodes around)
        boundary1_major_allele=((~nuc_present)<<4)|nuc_present;
        boundary2_allele=0;
    }
}
static uint8_t movemask(int in){
    return (1&in)|(2&(in>>3))|(4&(in>>6))|(8&(in>>9));
}
static void set_state_from_cnt(const std::array<int,4>& data, uint8_t& boundary1_major_allele_out,uint8_t& boundary2_allele){
    __m128i ori=_mm_loadu_si128((__m128i*)data.data());
    __m128i suf=_mm_shuffle_epi32(ori, 0x4e);
    __m128i max1=_mm_max_epi32(ori, suf);
    __m128i max1_suf=_mm_shuffle_epi32(max1, 0x11);
    __m128i max_values=_mm_max_epi32(max1, max1_suf);
    //int32_t max=_mm_extract_epi32(max_values, 0);
    __m128i one=_mm_set1_epi32(1);
    __m128i two=_mm_set1_epi32(2);
    int is_max_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(ori, max_values));
    int is_boundary_1_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_add_epi32(ori,one), max_values));
    int is_boundary_2_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_add_epi32(ori,two), max_values));
    uint8_t max_mask=movemask(is_max_raw);
    uint8_t boundary1_mask=movemask(is_boundary_1_raw);
    boundary2_allele=movemask(is_boundary_2_raw);
    boundary1_major_allele_out=max_mask|(boundary1_mask<<4);
}


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

void FS_backward_pass(const std::vector<MAT::Node*> bfs_ordered_nodes, std::vector<uint8_t>& boundary1_major_allele,std::vector<uint8_t>& boundary2_allele,const std::unordered_map<std::string, nuc_one_hot>& mutated,nuc_one_hot ref_nuc){
    for(long node_idx=bfs_ordered_nodes.size()-1;node_idx>=0;node_idx--){
        auto this_node=bfs_ordered_nodes[node_idx];
        auto child_size=this_node->children.size();

        //leaf node
        if (child_size==0) {
            auto iter=mutated.find(this_node->identifier);
            nuc_one_hot allele;
            if (iter==mutated.end()) {
                allele=ref_nuc;
            }else {
                allele=iter->second;
            }
            boundary1_major_allele[node_idx]=allele;
            boundary2_allele[node_idx]=0;
        }else if (child_size==1) {
            auto child_idx=this_node->children[0]->bfs_index;
            boundary2_allele[node_idx]=0;
            boundary1_major_allele[node_idx]=boundary1_major_allele[child_idx]&0xf;
        }else if (child_size==2) {
            size_t child_start_idx=this_node->children[0]->bfs_index;
            set_state_2(boundary1_major_allele[child_start_idx],boundary1_major_allele[child_start_idx+1],boundary1_major_allele[node_idx],boundary2_allele[node_idx]);
        }else {
            size_t child_start_idx=this_node->children[0]->bfs_index;
            size_t child_end_idx=child_start_idx+this_node->children.size();
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
            set_state_from_cnt(nuc_count, boundary1_major_allele[node_idx], boundary2_allele[node_idx]);
        }
        assert(boundary1_major_allele[node_idx]&0xf);
    }
}

static nuc_one_hot set_state(MAT::Node* this_node,uint8_t boundary1_major_allele,uint8_t boundary2_allele,nuc_one_hot par_state,const MAT::Mutation& base,tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>& output){
    nuc_one_hot this_state;
    bool need_add=false;
    if(boundary1_major_allele&par_state){
        this_state=par_state;
    }else {
        this_state=boundary1_major_allele&0xf;
        this_state=this_state.choose_first();
        assert(this_node->parent->children.size()>1);
        need_add=true;
    }

    nuc_one_hot major_allele(boundary1_major_allele&0xf);
    nuc_one_hot tie_allele=major_allele&(~this_state);
    if (this_node->is_leaf()&&tie_allele) {
        fputc('a', stderr);
    }
    nuc_one_hot boundary1_allele=boundary1_major_allele>>4;
    #ifdef DETAIL_DEBUG_FITCH_SANKOFF
    assert(major_allele&this_state);
    #endif

    if (major_allele.is_ambiguous()||(boundary1_allele|boundary2_allele)||need_add) {
        MAT::Mutation to_add(base);
        to_add.set_par_mut(par_state, this_state);
        if (this_node->children.size()<=1) {
            assert(!boundary1_allele);
            boundary1_allele=(~major_allele)&0xf;
        }
        to_add.set_auxillary(tie_allele,boundary1_allele, boundary2_allele);
        output.push_back(to_add);
    }

    return this_state;
}

static void FS_forward_pass(const std::vector<MAT::Node*> bfs_ordered_nodes,const std::vector<uint8_t>& boundary1_major_allele,const std::vector<uint8_t>& boundary2_allele,const MAT::Mutation & base,std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>& output){
    std::vector<nuc_one_hot> states(bfs_ordered_nodes.size());
    states[0]=set_state(bfs_ordered_nodes[0], boundary1_major_allele[0], boundary2_allele[0], base.get_ref_one_hot(), base,output[0]);
    for (size_t node_idx=1; node_idx<bfs_ordered_nodes.size(); node_idx++) {
        states[node_idx]=set_state(bfs_ordered_nodes[node_idx], boundary1_major_allele[node_idx], boundary2_allele[node_idx], states[bfs_ordered_nodes[node_idx]->parent->bfs_index], base,output[node_idx]);
    }
}

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

void Fitch_Sankoff_Whole_Tree(const std::vector<MAT::Node*> bfs_ordered_nodes,const MAT::Mutation & base,const std::unordered_map<std::string, nuc_one_hot>& mutated,std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>& output){
    std::vector<uint8_t> minor_major_allele(bfs_ordered_nodes.size()+8);
    std::vector<uint8_t> major_minor_count_difference(bfs_ordered_nodes.size()+8);

    FS_backward_pass(bfs_ordered_nodes,minor_major_allele,major_minor_count_difference,mutated,base.get_ref_one_hot());
    FS_forward_pass(bfs_ordered_nodes,minor_major_allele,major_minor_count_difference,base,output);
}

