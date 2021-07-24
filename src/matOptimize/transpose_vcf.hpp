#include <cstdint>
#include <tbb/pipeline.h>
#include <zlib.h>
#include <cassert>
#include <string>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>
//#include <vector>
#define MAX_SIZ 0x30000

struct partitioner{
    const uint8_t*& last_out;
    const uint8_t* const end;
    const uint8_t* operator()(tbb::flow_control& fc)const {
        if (last_out>=end) {
            fc.stop();
        }
        //assert(last_out<end);
        auto out=last_out;
        //fprintf(stderr, "%d\n",*(int*)out);
        unsigned int item_len=*(int*) out;
        last_out+=(item_len+4);
        return out; 
    }
};
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
template<typename output_t>
const uint8_t* parse_buffer(const uint8_t* in,output_t& out_all){
    std::string sample_name;
    while(*in!=0){
        sample_name.push_back(*in);
        in++;
    }
    auto& out=out_all.set_name(std::move(sample_name));
    in++;
    while(*in){
        int pos1=loadVariant(in);
        assert(*in);
        if (*(in+1)) {
            int pos2=loadVariant(in);
            out.add_Not_N(pos1,(*in)&0xf);
            out.add_Not_N(pos2,0xf&((*in)>>4));
        }else {
            out.add_Not_N(pos1,(*in)&0xf);
        }
        in++;
    }
    in++;
    
    while (*in) {
        int first=loadVariant(in);
        auto after_first=in;
        if (!(*in)) {
            out.add_N(first,first);
            break;
        }
        int second=loadVariant(in);
        if (first>second) {
            out.add_N(second,first);
        }else{
            out.add_N(first,first);
            in=after_first;
        }
    }
    in++;
    return in;
}
template<typename output_t>
struct printer{
    output_t& out;
    void operator()(const uint8_t* in) const{
        uint8_t buffer[MAX_SIZ];
        unsigned int item_len=*(int*) in;
        size_t out_len=MAX_SIZ;
        auto uncompress_out=uncompress(buffer, &out_len, in+4, item_len);
        assert(uncompress_out==Z_OK);
        const uint8_t* start=buffer;
        auto end=buffer+out_len;
        while(start!=end){
            start=parse_buffer(start,out);
        }
    }
};
template<typename output_t>
static void load_mutations(char* path,int nthread,output_t& out){
    auto fh=open(path, O_RDONLY);
    struct stat stat_buf;
    fstat(fh, &stat_buf);
    auto size=stat_buf.st_size;
    const uint8_t* file=(uint8_t*) mmap(nullptr, size, PROT_READ, MAP_SHARED, fh, 0);
    const uint8_t* last_out=file;
    auto end=file+size;
    tbb::parallel_pipeline(nthread,tbb::make_filter<void,const uint8_t*>(tbb::filter::serial_in_order,partitioner{last_out,end})&
        tbb::make_filter<const uint8_t*,void>(tbb::filter::parallel,printer<output_t>{out}));

    munmap((void*)file, size);
}