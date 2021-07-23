#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#define LOAD
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "zlib.h"
#include <tbb/task_scheduler_init.h>
#include "tbb/flow_graph.h"
#include <sys/mman.h>
#include <mutex>
#include "transpose_vcf.hpp"
std::mutex print_mutex;
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
#ifdef TEST_DECODER
struct Pos_Mut{
    int position;
    uint8_t mut;
    bool operator<(const Pos_Mut& other)const{
        return position<other.position;
    }
};
#endif
static const uint8_t* parse_buffer(const uint8_t* in,
#ifdef TEST_DECODE
 std::string& sample_out,
 std::vector<Pos_Mut>&
 #else
Original_State_t&
#endif
mutations
){
#ifndef TEST_DECODE
 std::string sample_out;
#endif
    while(*in!=0){
        sample_out.push_back(*in);
        in++;
    }
    in++;
    size_t mut_size=loadVariant(in);
#ifdef TEST_DECODE
    mutations.resize(mut_size);
#else
    auto sample_iter=mutations.find(sample_out);
    bool found=sample_iter!=mutations.end();
    if (found) {
        sample_iter->second.clear();
        sample_iter->second.reserve(mut_size);
    }
    Mutation_Annotated_Tree::Mutation mut1; 
    Mutation_Annotated_Tree::Mutation mut2; 
#endif
    for(size_t idx=0;idx<mut_size;idx+=2){
#ifdef TEST_DECODE
        mutations[idx].position=loadVariant(in);
#else
        mut1.position=loadVariant(in);
#endif
        if (idx+1==mut_size) {
#ifdef TEST_DECODE
            mutations[idx].mut=*in;
#else
            mut1.boundary1_all_major_allele=*in;
            if (found) {
                sample_iter->second.insert(mut1);
            }
#endif
            in++;
        }else{
#ifdef TEST_DECODE
            mutations[idx+1].position=loadVariant(in);
            mutations[idx].mut=*in&0xf;
            mutations[idx+1].mut=*in>>4;
#else
            mut2.position=loadVariant(in);
            mut1.boundary1_all_major_allele=*in&0xf;
            mut2.boundary1_all_major_allele=*in>>4;
            if (found) {
                sample_iter->second.insert(mut1);
                sample_iter->second.insert(mut2);
            }
#endif
            in++;
        }
    }
    return in;
}
struct partitioner{
    const uint8_t*& last_out;
    const uint8_t* const end;
    bool operator()(const uint8_t*& out)const {
        if (last_out>=end) {
            return false;
        }
        //assert(last_out<end);
        out=last_out;
        //fprintf(stderr, "%d\n",*(int*)out);
        unsigned int item_len=*(int*) out;
        last_out+=(item_len+4);
        return true; 
    }
};
#define MAX_SIZ 0x30000
struct printer{
    #ifndef TEST_DECODE
    Original_State_t& sample_muts;
    #endif
    void operator()(const uint8_t* in) const{
        uint8_t buffer[MAX_SIZ];
        unsigned int item_len=*(int*) in;
        size_t out_len=MAX_SIZ;
        auto uncompress_out=uncompress(buffer, &out_len, in+4, item_len);
        switch (uncompress_out) {
            case Z_MEM_ERROR:fprintf(stderr, "decompress memory error\n");
                break;
            case Z_DATA_ERROR:fprintf(stderr, "data corrupted\n");
                break;
            case Z_BUF_ERROR: fprintf(stderr, "buffer not enough");
                break;
        }
        if (uncompress_out!=Z_OK) {
            fprintf(stderr, "corrupted file size %d\n",item_len);
            assert(false);
        }
        const uint8_t* start=buffer;
        auto end=buffer+out_len;
        while(start!=end){
    #ifdef TEST_DECODE
            std::vector<Pos_Mut> mutations_test;
            std::string sample;
            start=parse_buffer(start,sample,mutations_test);
            std::string out(sample);
            out+=":";
            for(int pos_idx=0;pos_idx<mutations_test.size();pos_idx++){
                out+="\t";
                out+=std::to_string(mutations_test[pos_idx].position);
                out+=" ";
                out.push_back(Mutation_Annotated_Tree::get_nuc(mutations_test[pos_idx].mut));
            }
        {
            std::lock_guard<std::mutex> lk(print_mutex);
            puts(out.c_str());
        }
#else
            start=parse_buffer(start,sample_muts);
#endif
        }
            //fprintf(stderr, "finished file size %d\n",item_len);
    }
};
#ifdef TEST_DECODER
int main(int argc,char** argv){
    char* path=argv[1];
    //tbb::task_scheduler_init init(1);
#else
void load_mutations(char* path, Original_State_t& sample_muts){
#endif
    auto fh=open(path, O_RDONLY);
    struct stat stat_buf;
    fstat(fh, &stat_buf);
    auto size=stat_buf.st_size;
    const uint8_t* file=(uint8_t*) mmap(nullptr, size, PROT_READ, MAP_SHARED, fh, 0);
    auto end=file+size;
    tbb::flow::graph input_graph;
    tbb::flow::source_node<const uint8_t*> source (input_graph,partitioner{file,end});
    tbb::flow::function_node<const uint8_t*> consumer(input_graph,tbb::flow::unlimited,printer());
    tbb::flow::make_edge(source,consumer);
    input_graph.wait_for_all();
}