#include <cstdint>
#include <vector>
struct Nuc_Packer{
    std::vector<uint8_t>& content;
    void reserve(size_t size){
        content.reserve(size/2+1);
    }
    void push_back(uint8_t allele){
        if (content.empty()||content.back()&0xf0) {
            content.push_back(allele);
        }else{
            content.back()|=(allele<<4);
        }
    }
    uint8_t operator[](size_t idx){
        uint8_t two=content[idx>>1];
        if (idx&1) {
            return two>>4;
        }else {
            return two&0xf;
        }
    }
};