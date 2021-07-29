#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "zlib.h"
#include "tbb/flow_graph.h"
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <ios>
#include <string>
#include "transpose_vcf.hpp"
#include <sys/types.h>
#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include "tbb/parallel_for.h"
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
std::mutex print_mutex;
#define ZLIB_BUFSIZ 0x10000
#include <atomic>
#define SAMPLE_START_IDX 9
std::atomic<size_t> buffer_left; 
//Decouple parsing (slow) and decompression, segment file into blocks for parallelized parsing
typedef tbb::flow::source_node<char*> decompressor_node_t;
typedef tbb::flow::function_node<char*> line_parser_t;
static void writeVariant(std::vector<uint8_t>& out,unsigned int to_write){
    assert(to_write);
    out.push_back(to_write&0x7f);
    to_write>>=7;
    while(to_write){
        out.back()|=0x80;
        out.push_back(to_write&0x7f);
        to_write>>=7;
    }
}
struct Pos_Mut{
    int position;
    uint8_t mut;
    bool operator<(const Pos_Mut& other)const{
        return position<other.position;
    }
};
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
        //buffer_left++;
        //fprintf(stderr, "buffer in flight %zu\n",buffer_left.load());
        return buf;
    }
};

struct Pos_Mut_Block{
    std::vector<Pos_Mut> not_N;
    std::vector<std::pair<int,int>> Ns;
    void add_N(int position){
        if(Ns.empty()||(Ns.back().second+1)!=position){
            Ns.emplace_back(position,position);
        }else {
            Ns.back().second++;
        }
    }
    bool empty(){
        return not_N.empty()&&Ns.empty();
    }
    int min_pos;
    int max_pos;
    void finalize(){
        min_pos=std::min(not_N.empty()?INT_MAX:not_N[0].position,Ns.empty()?INT_MAX:Ns[0].first);
        assert(min_pos!=INT_MAX);
        max_pos=std::max(not_N.empty()?0:not_N.back().position,Ns.empty()?0:Ns.back().second);
        assert(max_pos!=0);
    }

    bool operator<(const Pos_Mut_Block& other)const {
        return max_pos<other.min_pos;
    }
};
//Parse a block of lines, assuming there is a complete line in the line_in buffer
struct Line_Parser{
    size_t sample_size;
    std::vector<Pos_Mut_Block>* operator()(char* line_in)const{
        char* const start=line_in;
        std::vector<Pos_Mut_Block>* local_block= new std::vector<Pos_Mut_Block>(sample_size);
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
            for (; field_idx < SAMPLE_START_IDX; field_idx++) {
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
                //output prototype of mutation, and a map from sample to non-ref allele
                if (allele_idx>=(allele_translated.size()+1)) {
                    (*local_block)[field_idx].add_N(pos);
                }else if (allele_idx) {
                    (*local_block)[field_idx].not_N.push_back(Pos_Mut{pos,allele_translated[allele_idx-1]});
                }
                field_idx++;
                line_in++;
            }
        }
        delete[] (start);
        //buffer_left--;
        //fprintf(stderr, "buffer in flight %zu\n",buffer_left.load());
        return local_block;
    }
};
struct Appender{
    std::vector<std::vector<Pos_Mut_Block>>& sample_pos_mut;
    void operator()(std::vector<Pos_Mut_Block>* in) const{
        for (size_t idx=0; idx<sample_pos_mut.size(); idx++) {
            if (!(*in)[idx].empty()) {
                sample_pos_mut[idx].push_back(std::move((*in)[idx]));
            }
        }
        delete in;
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
struct Sample_Mut_Msg{
    std::vector<uint8_t> buffer;
    Sample_Mut_Msg(const std::string& sample, const std::vector<Pos_Mut>& not_Ns,const std::vector<std::pair<int, int>>& Ns){
        buffer.reserve(sample.size()+5*not_Ns.size()+8*Ns.size());
        buffer.insert(buffer.end(),(const uint8_t*)sample.c_str(),(const uint8_t*)(sample.c_str()+1+sample.size()));
        assert(buffer.back()==0);
        auto loop_end=not_Ns.size()&0xfffffffe;
        for(size_t idx=0;idx<loop_end;idx+=2){
            writeVariant(buffer, not_Ns[idx].position);
            assert(idx+1<not_Ns.size());
            writeVariant(buffer, not_Ns[idx+1].position);
            buffer.push_back((not_Ns[idx+1].mut<<4)|not_Ns[idx].mut);
        }
        if (not_Ns.size()&1) {
            writeVariant(buffer, not_Ns[not_Ns.size()-1].position);
            buffer.push_back(not_Ns[not_Ns.size()-1].mut);
        }
        buffer.push_back(0);
        for(size_t idx=0;idx<Ns.size();idx++){
            writeVariant(buffer, Ns[idx].second);
            if (Ns[idx].first!=Ns[idx].second) {
                writeVariant(buffer, Ns[idx].first);
            }
        }
        buffer.push_back(0);
    }
};

Sample_Mut_Msg* serialize(const std::string& sample, std::vector<Pos_Mut_Block>& mutation_blocks){
    for (auto& block : mutation_blocks) {
        block.finalize();
    }
    std::sort(mutation_blocks.begin(),mutation_blocks.end());
    std::vector<std::pair<int,int>> Ns;
    std::vector<Pos_Mut> not_Ns;
    Ns.reserve(mutation_blocks.size());
    not_Ns.reserve(mutation_blocks.size());
    for (const auto& block : mutation_blocks) {
        not_Ns.insert(not_Ns.end(),block.not_N.begin(),block.not_N.end());
        if (block.Ns.empty()) {
            continue;
        }
        auto iter=block.Ns.begin();
        if ((!Ns.empty())&&Ns.back().second==(iter->first-1)) {
            Ns.back().second=iter->second;
            iter++;
        }
        Ns.insert(Ns.end(),iter,block.Ns.end());
    }
    Sample_Mut_Msg* out=new Sample_Mut_Msg(sample,not_Ns,Ns);
 #ifdef TEST_IN_ENCODE   
    std::vector<Pos_Mut> not_Ns_test;
    std::string sample_name_test;
    std::vector<std::pair<int,int>> Ns_test;
    auto end_ptr=parse_buffer(out->buffer.data(),sample_name_test,not_Ns_test,Ns_test);
    assert(sample_name_test==sample);
    assert(not_Ns_test.size()==not_Ns.size());
    for (size_t idx=0;idx<not_Ns.size();idx++) {
        assert(not_Ns_test[idx].mut==not_Ns[idx].mut);
        assert(not_Ns_test[idx].position==not_Ns[idx].position);
    }
    assert(Ns_test.size()==Ns.size());
    for (size_t idx=0;idx<Ns.size();idx++) {
        assert(Ns_test[idx]==Ns[idx]);
    }
    assert(end_ptr==(&(out->buffer.back())+1));
#endif
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
    std::pair<unsigned char*,size_t> operator()(Packed_Msgs* in) const {
        uint8_t* out=new uint8_t[compress_len];
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
    void operator()(std::pair<unsigned char*,size_t> in) const {
        unsigned int b_length=in.second;
        //fprintf(stderr,"%d\n",b_length);
        std::fwrite(&b_length,4,1,file);
        std::fwrite(in.first,1,b_length,file);
        delete[] in.first;
    }
};
void get_samp_names(const std::string &sample_names_fn,const std::vector<std::string>& fields,std::vector<bool>& do_add) {
    std::unordered_set<std::string> sample_set;
    std::ifstream infile(sample_names_fn, std::ios::in | std::ios::binary);
    if (infile.good()) {
        boost::iostreams::filtering_istream instream;
        instream.push(boost::iostreams::gzip_decompressor());
        instream.push(infile);
        std::string sample;
        std::getline(instream, sample);
        sample_set.reserve(std::stoi(sample));
        while (instream.good()) {
            sample.clear();
            std::getline(instream, sample);
            sample_set.insert(sample);
        }
    }
    do_add.resize(fields.size());
    for (size_t idx=SAMPLE_START_IDX; idx<fields.size(); idx++) {
        auto res=sample_set.insert(fields[idx]);
        do_add[idx]=res.second;
    }
    infile.close();
    std::ofstream outfile(sample_names_fn, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
    outbuf.push(boost::iostreams::gzip_compressor());
    outbuf.push(outfile);
    std::ostream  out(&outbuf);
    out<<sample_set.size()<<'\n';
    for (const auto& sample_name : sample_set) {
        out<<sample_name<<'\n';
    }
}
#define CHUNK_SIZ 5
void VCF_input(const char * name,const char * out_name,uint32_t nthreads, const std::string& sample_names_fn){
    std::vector<std::string> fields;
    //open file set increase buffer size
    gzFile fd=gzopen(name, "r");
    if (!fd) {
        fprintf(stderr, "cannnot open vcf file : %s, exiting.\n",name);
        exit(EXIT_FAILURE);
    }
    gzbuffer(fd,ZLIB_BUFSIZ);

    unsigned int header_size=read_header(&fd, fields);
    std::vector<bool> do_add;
    std::thread sample_name_thread(get_samp_names,std::ref(sample_names_fn), std::ref(fields), std::ref(do_add));
    std::vector<std::vector<Pos_Mut_Block>> sample_pos_mut(fields.size());
    for(auto& samp:sample_pos_mut){
        samp.reserve(30);
    }

    tbb::parallel_pipeline(nthreads,
        tbb::make_filter<void,char*>(tbb::filter::serial_in_order,Decompressor{&fd,CHUNK_SIZ*header_size,2*header_size})&
        tbb::make_filter<char*,std::vector<Pos_Mut_Block>*>(tbb::filter::parallel,Line_Parser{fields.size()})&
        tbb::make_filter<std::vector<Pos_Mut_Block>*,void>(tbb::filter::serial_out_of_order,Appender{sample_pos_mut}));
    gzclose(fd);
    sample_name_thread.join();
    compress_len=compressBound(MAX_SIZ);
    tbb::flow::graph output_graph;
    Packed_Msgs* blk_str=new Packed_Msgs;
    block_serializer_t serializer_head(output_graph,tbb::flow::serial,Block_Serializer{blk_str});
    compressor_t compressor(output_graph,4,Compressor{});
    FILE* out_file=fopen(out_name, "a");
    write_node_t writer(output_graph,tbb::flow::serial,Write_Node{out_file});
    tbb::flow::make_edge(serializer_head,compressor);
    tbb::flow::make_edge(compressor,writer);
    tbb::parallel_for(tbb::blocked_range<size_t>(SAMPLE_START_IDX,fields.size()),[&fields,&sample_pos_mut,&serializer_head,&compressor,&do_add](tbb::blocked_range<size_t>& range){
        auto packed_out=new Packed_Msgs();
        for(size_t idx=range.begin();idx<range.end();idx++){
            if (!do_add[idx]) {
                continue;
            }
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
namespace po = boost::program_options;

int main(int argc,char** argv){
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string input_vcf_path;
    std::string output_path;
    uint32_t num_threads;
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
        ("vcf,v", po::value<std::string>(&input_vcf_path)->default_value(""), "Input VCF file (in uncompressed or gzip-compressed .gz format) ")
        ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("output_path,o", po::value<std::string>(&output_path)->default_value(""), "Save transposed VCF");
    po::options_description all_options;
    all_options.add(desc);
    po::variables_map vm;
   try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        // Return with error code 1 unless the user specifies help
        std::cerr << desc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }
    tbb::task_scheduler_init init(num_threads);
    VCF_input(input_vcf_path.c_str(), output_path.c_str(),num_threads,output_path+".sample_names.gz");
    std::flush(std::cout);
}
