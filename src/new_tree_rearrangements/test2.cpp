#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <emmintrin.h>
#include <nmmintrin.h>
#include <array>
#include <cassert>
#include <smmintrin.h>

static uint8_t movemask(int in){
    return (1&in)|(2&(in>>3))|(4&(in>>6))|(8&(in>>9));
}
void set_state_from_cnt(const std::array<int,4> in){

    __m128i ori=_mm_loadu_si128((__m128i*)in.data());
    __m128i suf=_mm_shuffle_epi32(ori, 0x4e);
    __m128i max1=_mm_max_epi32(ori, suf);
    __m128i max1_suf=_mm_shuffle_epi32(max1, 0x11);
    __m128i max_values=_mm_max_epi32(max1, max1_suf);
    int32_t max=_mm_extract_epi32(max_values, 0);
    int32_t ref_max=*std::max_element(in.begin(),in.end());
    assert(max==ref_max);
    assert(max==_mm_extract_epi32(max_values, 1));
    assert(max==_mm_extract_epi32(max_values, 2));
    assert(max==_mm_extract_epi32(max_values, 3));
    __m128i one=_mm_set1_epi32(1);
    __m128i two=_mm_set1_epi32(2);
    int is_max_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(ori, max_values));
    int is_boundary_1_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_add_epi32(ori,one), max_values));
    int is_boundary_2_raw=_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_add_epi32(ori,two), max_values));
    uint8_t max_mask=movemask(is_max_raw);
    uint8_t boundary1_mask=movemask(is_boundary_1_raw);
    uint8_t boundary2_mask=movemask(is_boundary_2_raw);
    assert(!(max_mask&boundary1_mask));
    assert(!(max_mask&boundary2_mask));
    assert(!(boundary1_mask&boundary2_mask));
    for (int i=0; i<4; i++) {
        int mask=1<<i;
        if(mask&max_mask){
            assert(in[i]==max);
        }
        else if(mask&boundary1_mask){
            assert(in[i]==max-1);
        }
        else if (mask&boundary2_mask) {
            assert(in[i]==max-2);
        }else{
            assert(in[i]<max-2);
        }
    }
    /*
    __m128i first2=_mm_load_si128((__m128i*)data);
    __m128i second2=_mm_load_si128((__m128i*)(&data[2]));
    __m128i max1mask=_mm_cmpgt_epi64(first2,second2);
    __m128i max1=_mm_or_si128(_mm_and_si128(max1mask, first2), _mm_andnot_si128(max1mask, second2));
    __m128i max1flp=_mm_s
    */
}

void tester(size_t idx,std::array<int,4> in){
    if (idx<4) {
        for (int to_try=0;to_try<4 ; to_try++) {
            in[idx]=to_try;
            tester(idx+1, in);
        }
    }else {
        if (in!=std::array<int,4>({0,0,0,0})){
            set_state_from_cnt(in);
        }
    }
}
int main(int argc,char** argv){
/*    char bytes[]={15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};
    __m128i child_states_after=_mm_loadu_si128((__m128i*)bytes);
    fprintf(stderr, "1:%llx,2:%llx\n",child_states_after[0],child_states_after[1]);
    int to_calculate=23;
    unsigned long unloadmask1=-1;
    unsigned long unloadmask2=-1;
    if (to_calculate<32) {
        if (to_calculate>24) {
            unloadmask1>>=(8*(32-to_calculate));
        }else {
            assert(to_calculate>16);
            assert(to_calculate<=24);
            unloadmask1=0;
            unloadmask2>>=(8*(24-to_calculate));
        }
    }
    __m128i unload_mask=_mm_set_epi64x(unloadmask1, unloadmask2);
    child_states_after=_mm_and_si128(child_states_after, unload_mask);
    fprintf(stderr, "1:%llx,2:%llx\n",child_states_after[0],child_states_after[1]);
*/
    tester(0, std::array<int, 4>({0,0,0,0}));
}