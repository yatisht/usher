#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include "src/matOptimize/mutation_annotated_tree.hpp"
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
struct Pos_Mut{
    int position;
    uint8_t mut;
    bool operator<(const Pos_Mut& other)const{
        return position<other.position;
    }
};
static const uint8_t* parse_buffer(const uint8_t* in,std::string& sample_out, std::vector<Pos_Mut>& mutations){
    while(*in!=0){
        sample_out.push_back(*in);
        in++;
    }
    in++;
    size_t mut_size=loadVariant(in);
    mutations.resize(mut_size);
    for(size_t idx=0;idx<mut_size;idx+=2){
        mutations[idx].position=loadVariant(in);
        if (idx+1==mut_size) {
            mutations[idx].mut=*in;
            in++;
        }else{
            mutations[idx+1].position=loadVariant(in);
            mutations[idx].mut=*in&0xf;
            mutations[idx+1].mut=*in>>4;
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
#define MAX_SIZ 0xffff
struct printer{
    mutable uint8_t buffer[MAX_SIZ];
    printer(){}
    printer(const printer&){}
    printer(printer&){}
    void operator=(const printer&){}
    void operator()(const uint8_t* in) const{
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
        /*google::protobuf::io::CodedInputStream data(buffer,out_len);
        data.SetTotalBytesLimit(out_len, out_len);
        while (data.CurrentPosition()!=out_len) {
            unsigned int limit;
            data.ReadVarint32(&limit);
            auto internal_limit=data.PushLimit(limit);
            TransposedVCF::sample sample;
            bool ret=sample.ParseFromCodedStream(&data);
            assert(ret);
            auto nPos=sample.position_size();
            Nuc_Packer nucs{*sample.mutable_allele()};*/
        const uint8_t* start=buffer;
        auto end=buffer+out_len;
        while(start!=end){
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
            //data.PopLimit(internal_limit);
        }
            //fprintf(stderr, "finished file size %d\n",item_len);
    }
};
int main(int argc,char** argv){
    tbb::task_scheduler_init init(1);
    auto fh=open(argv[1], O_RDONLY);
    struct stat stat_buf;
    fstat(fh, &stat_buf);
    auto size=stat_buf.st_size;
    const uint8_t* file=(uint8_t*) mmap(nullptr, size, PROT_READ, MAP_SHARED, fh, 0);
    auto end=file+size;
    tbb::flow::graph input_graph;
    tbb::flow::source_node<const uint8_t*> source (input_graph,partitioner{file,end});
    tbb::flow::function_node<const uint8_t*> consumer(input_graph,tbb::flow::serial,printer());
    tbb::flow::make_edge(source,consumer);
    input_graph.wait_for_all();
}