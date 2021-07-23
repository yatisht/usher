#include "transpose_vcf.hpp"
#include <cassert>
#include <string>
#include <vector>
static unsigned int loadVariant(const uint8_t*& in){
    unsigned int out=(*in&0x7f);
    int shamt=7;
    while (*in&0x80) {
        in++;
        out|=((*in&0x7f)<<shamt);
        shamt+=7;
    }
    in++;
    return out;
}
const uint8_t* parse_buffer(const uint8_t* in,std::string& sample_out, std::vector<Pos_Mut>& mutations,std::vector<std::pair<int,int>>& Ns){
    while(*in!=0){
        sample_out.push_back(*in);
        in++;
    }
    in++;
    while(*in){
        int pos1=loadVariant(in);
        assert(*in);
        if (*(in+1)) {
            int pos2=loadVariant(in);
            mutations.push_back(Pos_Mut{pos1,static_cast<uint8_t>((*in)&0xf)});
            mutations.push_back(Pos_Mut{pos2,static_cast<uint8_t>(0xf&((*in)>>4))});
        }else {
            mutations.push_back(Pos_Mut{pos1,(*in)});
        }
        in++;
    }
    in++;
    
    while (*in) {
        int first=loadVariant(in);
        auto after_first=in;
        if (!(*in)) {
            Ns.emplace_back(first,first);
            break;
        }
        int second=loadVariant(in);
        if (first>second) {
            Ns.emplace_back(second,first);
        }else{
            Ns.emplace_back(first,first);
            in=after_first;
        }
    }
    in++;
    return in;
}
