#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "zlib.h"
#include "tbb/flow_graph.h"
#include <cctype>
#include <cstdint>
#include <cstdio>
//#include <google/protobuf/io/coded_stream.h>
#include <string>
#include "transpose_vcf.hpp"
#include <sys/types.h>
#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <vector>
#include <iostream>
//#include "transposed_vcf.pb.h"
//#include <google/protobuf/io/zero_copy_stream_impl_lite.h>
#include "tbb/parallel_for.h"
std::mutex print_mutex;
#define ZLIB_BUFSIZ 0x10000
#include <atomic>
std::atomic<size_t> buffer_left; 
//Decouple parsing (slow) and decompression, segment file into blocks for parallelized parsing
typedef tbb::flow::source_node<char*> decompressor_node_t;
typedef tbb::flow::function_node<char*> line_parser_t;
static void writeVariant(std::vector<uint8_t>& out,unsigned int to_write){
    out.push_back(to_write&0x7f);
    to_write>>=7;
    while(to_write){
        out.back()|=0x80;
        out.push_back(to_write&0x7f);
        to_write>>=7;
    }
}
/*
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
*/
struct Decompressor{
    gzFile* fd;
    size_t init_read_size;
    size_t cont_read_size;
    char* operator()(tbb::flow_control& fc) const{
        char* buf;
        if (gzeof(*fd)) {
            fc.stop();
            return nullptr;
        }
        buf=new char[init_read_size+cont_read_size];
        int read_size=gzread(*fd, buf, init_read_size);
        if (read_size<0) {
            int z_errnum = 0;
            fputs( gzerror(*fd,&z_errnum),stderr);
        }
        if (!read_size) {
            delete [] (buf);
            fc.stop();
            return nullptr;
        }
        //Make sure the last line is complete in the block.
        if(!gzgets(*fd, buf+read_size, cont_read_size)){
            *(buf+read_size)=0;
        }
        buffer_left++;
        fprintf(stderr, "buffer in flight %zu\n",buffer_left.load());
        return buf;
    }
};
struct Pos_Mut{
    int position;
    uint8_t mut;
    bool operator<(const Pos_Mut& other)const{
        return position<other.position;
    }
};
//Parse a block of lines, assuming there is a complete line in the line_in buffer
struct Line_Parser{
    std::vector<tbb::concurrent_vector<Pos_Mut>>& sample_pos_mut;
    void operator()(char* line_in)const{
        char* const start=line_in;
        while (*line_in!=0) {
            std::vector<nuc_one_hot> allele_translated;
            int pos=0;
            //Chromosome
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //Position
            while (*line_in!='\t') {
                pos=pos*10+(*line_in-'0');
                line_in++;
            }
            assert(pos>0);
            line_in++;
            //ID don't care
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //REF
            line_in++;
            //assert(*line_in=='\t');
            line_in++;
            //ALT
            while (*line_in!='\t') {
                allele_translated.push_back(Mutation_Annotated_Tree::get_nuc_id(*line_in));
                line_in++;
                if(*line_in==','){
                    line_in++;
                }else{
                    //assert(*line_in=='\t');
                }
            }
            line_in++;
            unsigned int field_idx=5;
            for (; field_idx < 9; field_idx++) {
              while (*line_in != '\t') {
                line_in++;
              }
              line_in++;
            }
            //samples
            bool is_last=false;
            while (!is_last) {
                unsigned int allele_idx=(*line_in-'0');
                line_in++;
                while (std::isdigit(*line_in)) {
                    allele_idx*=10;
                    allele_idx+=(*line_in-'0');
                    line_in++;
                }
                while (*line_in!='\t') {
                    if (*line_in=='\n'||*line_in==0) {
                        is_last=true;
                        break;
                    }
                    line_in++;
                }
                if (!sample_pos_mut[field_idx].empty()) {
                    assert(sample_pos_mut[field_idx].back().position!=pos);
                }
                //output prototype of mutation, and a map from sample to non-ref allele
                if (allele_idx>=(allele_translated.size()+1)) {
                    sample_pos_mut[field_idx].push_back(Pos_Mut{pos,0xf});
                }else if (allele_idx) {
                    sample_pos_mut[field_idx].push_back(Pos_Mut{pos,allele_translated[allele_idx-1]});
                }
                field_idx++;
                line_in++;
            }
        }
        delete[] (start);
        buffer_left--;
        fprintf(stderr, "buffer in flight %zu\n",buffer_left.load());
    }
};
//tokenize header, get sample name
static int read_header(gzFile* fd,std::vector<std::string>& out){
    int header_len=0;
    char in=gzgetc(*fd);
    in=gzgetc(*fd);
    bool second_char_pong=(in=='#');

    while (second_char_pong) {
        while (in!='\n') {
            in=gzgetc(*fd);
        }
        in=gzgetc(*fd);
        in=gzgetc(*fd);
        second_char_pong=(in=='#');
    }

    bool eol=false;
    while (!eol) {
        std::string field;
        while (in!='\t') {
            if (in=='\n') {
                eol=true;
                break;
            }
            field.push_back(in);
            in=gzgetc(*fd);
            header_len++;
        }
        in=gzgetc(*fd);
        out.push_back(field);
    }
    return header_len;
}
#define MAX_SIZ 0x30000
struct Sample_Mut_Msg{
    std::vector<uint8_t> buffer;
    Sample_Mut_Msg(const std::string& sample, tbb::concurrent_vector<Pos_Mut>& mutations){
        buffer.reserve(sample.size()+5*mutations.size());
        buffer.insert(buffer.end(),(const uint8_t*)sample.c_str(),(const uint8_t*)(sample.c_str()+1+sample.size()));
        assert(buffer.back()==0);
        writeVariant(buffer, mutations.size());
        for(size_t idx=0;idx<mutations.size();idx+=2){
            writeVariant(buffer, mutations[idx].position);
            if (idx+1==mutations.size()) {
                buffer.push_back(mutations[idx].mut);
                break;
            }else{
                writeVariant(buffer, mutations[idx+1].position);
                buffer.push_back((mutations[idx+1].mut<<4)|mutations[idx].mut);
            }
        }
    }
};
/*
const uint8_t* parse_buffer(const uint8_t* in,std::string& sample_out, std::vector<Pos_Mut>& mutations){
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
*/
Sample_Mut_Msg* serialize(const std::string& sample, tbb::concurrent_vector<Pos_Mut>& mutations){
    std::sort(mutations.begin(),mutations.end());
    std::string p_out(sample);
    p_out+=":";
    for (auto const& mut : mutations) {
        p_out+="\t";
        p_out+=std::to_string(mut.position);
        p_out+=" ";
        p_out.push_back(Mutation_Annotated_Tree::get_nuc(mut.mut));
        //widx++;
    }
    {
        //std::lock_guard<std::mutex> lk(print_mutex);
        //std::cout<<p_out<<"\n";
    }
    /*TransposedVCF::sample* out=new TransposedVCF::sample();
    out->set_sample_name(sample);
    out->mutable_position()->Reserve(mutations.size());
    Nuc_Packer muts{*out->mutable_allele()};
    muts.reserve(mutations.size());
    //size_t widx=0;
    */
    /*for (size_t idx=0;idx<mutations.size();idx++) {
        assert(muts[idx]==mutations[idx].mut);
    }*/
    Sample_Mut_Msg* out=new Sample_Mut_Msg(sample,mutations);
    
/*    std::vector<Pos_Mut> mutations_test;
    std::string sample_name_test;
    auto end_ptr=parse_buffer(out->buffer.data(),sample_name_test,mutations_test);
    assert(sample_name_test==sample);
    assert(mutations_test.size()==mutations.size());
    for (size_t idx=0;idx<mutations.size();idx++) {
        assert(mutations_test[idx].mut==mutations[idx].mut);
        assert(mutations_test[idx].position==mutations[idx].position);
    }
    assert(end_ptr==(&(out->buffer.back())+1));
    */
    return out;
}
struct Packed_Msgs:public std::vector<Sample_Mut_Msg*>{
    size_t acc_size;
    Packed_Msgs():acc_size(0){clear();}
    bool push_back(Sample_Mut_Msg* in){
        auto new_size=acc_size+in->buffer.size();
        if (new_size<=MAX_SIZ) {
            std::vector<Sample_Mut_Msg*>::push_back(in);
            acc_size=new_size;
            return true;
        }
        return false;
    }
    void pop_back(){
        acc_size-=(back()->buffer.size());
        std::vector<Sample_Mut_Msg*>::pop_back();
    }
/*    size_t serialize(uint8_t* out){
        google::protobuf::io::ArrayOutputStream out_wrap(out,MAX_SIZ);
        {
        google::protobuf::io::CodedOutputStream write_buffer(&out_wrap);
        for(auto msg:*this){
            write_buffer.WriteVarint32(msg->GetCachedSize());
            msg->SerializeToCodedStream(&write_buffer);
            delete msg;
        }
        }
        return out_wrap.ByteCount();
    }*/
};
typedef tbb::flow::multifunction_node<Packed_Msgs*, tbb::flow::tuple<Packed_Msgs*>> block_serializer_t;
struct Block_Serializer{
    Packed_Msgs* & out_buffer;
    void operator()(Packed_Msgs * in,block_serializer_t::output_ports_type& out ) const{
        if (!in) {
            if (!out_buffer->empty()) {            
                std::get<0>(out).try_put(out_buffer);
            }
            out_buffer=nullptr;
            return;
        }
        assert(out_buffer);
        while(!in->empty()){
            if (out_buffer->push_back(in->back())) {
                in->pop_back();
            }else {
                break;
            }
        }
        if (in->empty()) {
            delete in;
        }else{
            std::get<0>(out).try_put(out_buffer);
            out_buffer=in;
        }
    }
};
size_t compress_len;
typedef tbb::flow::function_node<Packed_Msgs*,std::pair<unsigned char*,size_t>> compressor_t;
struct Compressor{
    std::pair<unsigned char*,size_t> operator()(Packed_Msgs* in){
        uint8_t* out=new uint8_t[compress_len];
        /*auto raw_size=in->serialize(buffer);
        google::protobuf::io::CodedInputStream data(buffer,raw_size);
        while (data.CurrentPosition()!=raw_size) {
            unsigned int limit;
            data.ReadVarint32(&limit);
            auto internal_limit=data.PushLimit(limit);
            TransposedVCF::sample sample;
            bool ret=sample.ParseFromCodedStream(&data);
            assert(ret);
            data.PopLimit(internal_limit);
        }*/
        z_stream stream;
        stream.zalloc = (alloc_func)0;
        stream.zfree = (free_func)0;
        stream.opaque = (voidpf)0;
        auto err = deflateInit(&stream, Z_DEFAULT_COMPRESSION);
        assert(err==Z_OK);
        stream.next_out=out;
        stream.avail_out=compress_len;
        for(size_t idx=0;idx<in->size()-1;idx++){
            stream.next_in=(*in)[idx]->buffer.data();
            stream.avail_in=(*in)[idx]->buffer.size();
            err=deflate(&stream,Z_NO_FLUSH);
            assert(err==Z_OK);
            assert(stream.avail_out);
            assert(!stream.avail_in);
            delete (*in)[idx];
        }
        stream.next_in=(*in)[in->size()-1]->buffer.data();
        stream.avail_in=(*in)[in->size()-1]->buffer.size();
        err=deflate(&stream,Z_FINISH);
        assert(err==Z_STREAM_END);
        assert(!stream.avail_in);
        delete (*in)[in->size()-1];
        deflateEnd(&stream);
        size_t comp_len=stream.total_out;

/*        uint8_t buffer[MAX_SIZ];
        size_t iiii=MAX_SIZ;
        auto uncompress_out=uncompress(buffer, &iiii, out, comp_len);
        if (uncompress_out!=Z_OK) {
            sleep(1);
            assert(false);
        }
        auto msg_count=in->size();
        const uint8_t* start=buffer;
        auto end=buffer+iiii;
        while(start!=end){
            std::vector<Pos_Mut> mutations_test;
            std::string sample_name_test;
            start=parse_buffer(start,sample_name_test,mutations_test);
            msg_count--;
        }
        assert(msg_count==0);
*/
        delete in;
        return std::make_pair(out,comp_len);
    }
};
typedef tbb::flow::function_node<std::pair<unsigned char*,size_t>> write_node_t;
struct Write_Node{
    FILE* file;
    void operator()(std::pair<unsigned char*,size_t> in){
        unsigned int b_length=in.second;
        //fprintf(stderr,"%d\n",b_length);
        std::fwrite(&b_length,4,1,file);
        std::fwrite(in.first,1,b_length,file);
        delete[] in.first;
    }
};

#define CHUNK_SIZ 5
void VCF_input(const char * name,const char * out_name){
    std::vector<std::string> fields;
    //open file set increase buffer size
    gzFile fd=gzopen(name, "r");
    if (!fd) {
        fprintf(stderr, "cannnot open vcf file : %s, exiting.\n",name);
        exit(EXIT_FAILURE);
    }
    gzbuffer(fd,ZLIB_BUFSIZ);

    unsigned int header_size=read_header(&fd, fields);
    std::vector<tbb::concurrent_vector<Pos_Mut>> sample_pos_mut(fields.size());
    /*tbb::flow::graph input_graph;
    decompressor_node_t decompressor(input_graph,Decompressor{&fd,CHUNK_SIZ*header_size,2*header_size});
    line_parser_t parser(input_graph,tbb::flow::unlimited,Line_Parser{sample_pos_mut});
    //feed used buffer back to decompressor
    tbb::flow::make_edge(decompressor,parser);
    input_graph.wait_for_all();*/
    tbb::parallel_pipeline(80,tbb::make_filter<void,char*>(tbb::filter::serial_in_order,Decompressor{&fd,CHUNK_SIZ*header_size,2*header_size})&tbb::make_filter<char*,void>(tbb::filter::parallel,Line_Parser{sample_pos_mut}));
    gzclose(fd);
    compress_len=compressBound(MAX_SIZ);
    tbb::flow::graph output_graph;
    Packed_Msgs* blk_str=new Packed_Msgs;
    block_serializer_t serializer_head(output_graph,tbb::flow::serial,Block_Serializer{blk_str});
    compressor_t compressor(output_graph,4,Compressor{});
    FILE* out_file=fopen(out_name, "a");
    write_node_t writer(output_graph,tbb::flow::serial,Write_Node{out_file});
    tbb::flow::make_edge(serializer_head,compressor);
    tbb::flow::make_edge(compressor,writer);
    tbb::parallel_for(tbb::blocked_range<size_t>(9,fields.size()),[&fields,&sample_pos_mut,&serializer_head,&compressor](tbb::blocked_range<size_t>& range){
        auto packed_out=new Packed_Msgs();
        for(size_t idx=range.begin();idx<range.end();idx++){
            auto out=serialize(fields[idx],sample_pos_mut[idx]);
            if (!packed_out->push_back(out)) {
                compressor.try_put(packed_out);
                packed_out=new Packed_Msgs();
                auto ret=packed_out->push_back(out);
                if (!ret) {
                    fprintf(stderr, "Message size %zu, cannot fit %d, new content %zu\n",out->buffer.size(),MAX_SIZ,packed_out->acc_size);
                }
                assert(ret);
            }
        }
        serializer_head.try_put(packed_out);

    });
    serializer_head.try_put(nullptr);
    output_graph.wait_for_all();
    fclose(out_file);
}
int main(int argc,char** argv){
    /*
    std::vector<uint8_t> test_vect;
    writeVariant(test_vect, -1);
    assert(test_vect.size()==5);
    const uint8_t* ptr=test_vect.data();
    auto res=loadVariant(ptr);
    assert(0xffffffff==res);
    test_vect.clear();
    writeVariant(test_vect, 0xffff);
    assert(test_vect.size()==3);
    ptr=test_vect.data();
    assert(0xffff==loadVariant(ptr));
    test_vect.clear();
    writeVariant(test_vect, 0xfff);
    assert(test_vect.size()==2);
    ptr=test_vect.data();
    assert(0xfff==loadVariant(ptr));
    //tbb::task_scheduler_init init(1);
    */
    VCF_input(argv[1], "vcf_test_out.pb");
    std::flush(std::cout);
}
