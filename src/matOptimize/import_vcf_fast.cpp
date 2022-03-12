#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <mutex>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unordered_map>
#include "isa-l/igzip_lib.h"
#include <signal.h>
#include <unordered_set>
#include <vector>
#include "import_vcf.hpp"
#include "tbb/parallel_for.h"
#define ZLIB_BUFSIZ 0x10000
size_t read_size;
size_t alloc_size;
//Decouple parsing (slow) and decompression, segment file into blocks for parallelized parsing
typedef tbb::flow::source_node<char*> decompressor_node_t;
std::condition_variable progress_bar_cv;
struct raw_input_source {
    FILE* fh;
    raw_input_source(const char* fname) {
        fh=fopen(fname, "r");
    }
    int getc() {
        return fgetc(fh);
    }
    void operator()(tbb::concurrent_bounded_queue<std::pair<char*,uint8_t*>>& out) const {
        while(!feof(fh)) {
            auto line_out=new char[alloc_size];
            auto bytes_read=fread(line_out, 1, read_size, fh);
            if (bytes_read==0) {
                delete[] line_out;
                break;
            }
            out.emplace(line_out,(uint8_t*)line_out+bytes_read);
        }
        out.emplace(nullptr,nullptr);
        out.emplace(nullptr,nullptr);
    }
    void unalloc() {
        fclose(fh);
    }
};
struct line_start_later {
    char* start;
    char* alloc_start;
};
struct line_align {
    mutable line_start_later prev;
    mutable uint8_t* prev_end;
    tbb::concurrent_bounded_queue<std::pair<char*,uint8_t*>>& in;
    line_align(tbb::concurrent_bounded_queue<std::pair<char*,uint8_t*>>& out):in(out) {
        prev.start=nullptr;
    }
    bool operator()(line_start_later& out) const {
        std::pair<char*, unsigned char*> line;
        in.pop(line);
        if (line.first==nullptr) {
            if (prev.start!=nullptr) {
                out=prev;
                prev.start=nullptr;
                prev_end[0]=0;
                return true;
            } else {
                return false;
            }
        }
        if (!prev.start) {
            prev.start=line.first;
            prev.alloc_start=line.first;
            prev_end=line.second;
            in.pop(line);
            if (line.first==nullptr) {
                out=prev;
                prev.start=nullptr;
                return true;
            }
        }
        auto start_ptr=strchr(line.first, '\n');
        if (*start_ptr!='\n') {
            raise(SIGTRAP);
        }
        start_ptr++;
        auto cpy_siz=start_ptr-line.first;
        memcpy(prev_end, line.first, cpy_siz);
        prev_end[cpy_siz]=0;
        out=prev;
        prev.start=start_ptr;
        prev.alloc_start=line.first;
        prev_end=line.second;
        return true;
    }
};
struct gzip_input_source {
    unsigned char* map_start;
    //char* read_curr;
    size_t mapped_size;
    struct isal_gzip_header gz_hdr;
    mutable unsigned char* getc_buf;
    mutable unsigned char* get_c_ptr;
    struct inflate_state* state;
    gzip_input_source(const char* fname) {
        struct stat stat_buf;
        stat(fname, &stat_buf);
        mapped_size=stat_buf.st_size;
        auto fh=open(fname, O_RDONLY);
        map_start=(unsigned char*)mmap(nullptr, mapped_size, PROT_READ, MAP_SHARED, fh, 0);
        //madvise(map_start, mapped_size, MADV_SEQUENTIAL|MADV_WILLNEED);
        close(fh);
        state=new struct inflate_state;
        getc_buf=new unsigned char[BUFSIZ];
        get_c_ptr=0;

        isal_inflate_init(state);
        state->next_in=map_start;
        state->avail_in=mapped_size;
        state->crc_flag=IGZIP_GZIP;
        isal_gzip_header_init(&gz_hdr);
        auto ret = isal_read_gzip_header(state, &gz_hdr);
        fprintf(stderr, "Header ret: %d\n",ret);
        decompress_to_buffer(getc_buf, BUFSIZ);
        get_c_ptr=getc_buf;
    }
    void unalloc() {
        munmap(map_start, mapped_size);
    }
    bool decompress_to_buffer(unsigned char* buffer, size_t buffer_size) const {
        if (!state->avail_in) {
            return false;
        }
        state->next_out=buffer;
        state->avail_out=buffer_size;
        auto out=ISAL_DECOMP_OK;
        //while (out==ISAL_DECOMP_OK&&state->avail_out>0&&state->avail_in>0&&state->block_state!=ISAL_BLOCK_FINISH) {
        out=isal_inflate(state);
        //}
        if (out!=ISAL_DECOMP_OK) {
            fprintf(stderr, "decompress error %d\n",out);
        }
        //copied from https://github.com/intel/isa-l/blob/78f5c31e66fedab78328e592c49eefe4c2a733df/programs/igzip_cli.c#L952
        while (state->avail_out>0&&state->avail_in>0&&state->next_in[0] == 31) {
            //fprintf(stderr,"continuing\n");
            if (state->avail_in > 1 && state->next_in[1] != 139)
                break;
            isal_inflate_reset(state);
            //while (out==ISAL_DECOMP_OK&&state->avail_out>0&&state->avail_in>0&&state->block_state!=ISAL_BLOCK_FINISH) {
            out=isal_inflate(state);
            //}
            if (out!=ISAL_DECOMP_OK) {
                fprintf(stderr, "restart error %d\n",out);
            }
        }
        return state->avail_out<buffer_size;
    }
    int getc() {
        if (get_c_ptr==state->next_out) {
            get_c_ptr=getc_buf;
            if(!decompress_to_buffer(getc_buf, BUFSIZ)) {
                return -1;
            }
        }
        return *(get_c_ptr++);
    }

    void operator()(tbb::concurrent_bounded_queue<std::pair<char*,uint8_t*>>& out) const {
        char* line_out=new char[alloc_size];
        //if (getc_buf) {
        auto load_size=getc_buf+BUFSIZ-get_c_ptr;
        memcpy(line_out, get_c_ptr, load_size);
        decompress_to_buffer((unsigned char*)line_out+load_size, read_size-load_size);
        delete[](getc_buf);
        getc_buf=nullptr;
        fprintf(stderr,"getc_buff_deallocated\n");
        out.emplace(line_out,state->next_out);
        //}
        line_out=new char[alloc_size];
        while (decompress_to_buffer((unsigned char*)line_out, read_size)) {
            out.emplace(line_out,state->next_out);
            line_out=new char[alloc_size];
        }
        out.emplace(nullptr,nullptr);
        out.emplace(nullptr,nullptr);
        delete [] line_out;
    }
};
//Parse a block of lines, assuming there is a complete line in the line_in buffer
typedef tbb::flow::multifunction_node<line_start_later,tbb::flow::tuple<Parsed_VCF_Line*>> line_parser_t;
struct line_parser {
    const std::vector<long>& header;
    void operator()(line_start_later line_in_struct, line_parser_t::output_ports_type& out)const {
        char*  line_in=line_in_struct.start;
        char*  to_free=line_in_struct.alloc_start;
        Parsed_VCF_Line* parsed_line;
        while (*line_in!=0) {
            std::vector<nuc_one_hot> allele_translated;
            std::string chromosome;
            int pos=0;
            //Chromosome
            while (*line_in!='\t') {
                chromosome.push_back(*line_in);
                line_in++;
            }
            line_in++;
            //Position
            while (*line_in!='\t') {
                pos=pos*10+(*line_in-'0');
                line_in++;
            }
            if(pos<=0) {
                raise(SIGTRAP);
            }
            line_in++;
            //ID don't care
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //REF
            parsed_line=new Parsed_VCF_Line{MAT::Mutation(chromosome,pos, 0, 0, 1,
                                            MAT::get_nuc_id(*line_in))};
            auto& non_ref_muts_out=parsed_line->mutated;
            line_in++;
            //assert(*line_in=='\t');
            line_in++;
            //ALT
            while (*line_in!='\t') {
                allele_translated.push_back(MAT::get_nuc_id(*line_in));
                line_in++;
                if(*line_in==',') {
                    line_in++;
                } else {
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
                if(header[field_idx]>0) {
                    //output prototype of mutation, and a map from sample to non-ref allele
                    if (allele_idx>=(allele_translated.size()+1)) {
                        non_ref_muts_out.emplace_back(header[field_idx],0xf);
                    } else if (allele_idx&&header[field_idx]>=0) {
                        non_ref_muts_out.emplace_back(header[field_idx],allele_translated[allele_idx-1]);
                    }
                }
                field_idx++;
                line_in++;
            }
            //assert(field_idx==header.size());
            std::sort(non_ref_muts_out.begin(),non_ref_muts_out.end(),mutated_t_comparator());
            non_ref_muts_out.emplace_back(0,0);
            std::get<0>(out).try_put(parsed_line);
        }
        delete[] (to_free);
    }
};
//tokenize header, get sample name
template<typename infile_t>
static void read_header(infile_t& fd,std::vector<std::string>& out) {
    char in=fd.getc();
    in=fd.getc();
    bool second_char_pong=(in=='#');

    while (second_char_pong) {
        while (in!='\n') {
            in=fd.getc();
        }
        in=fd.getc();
        in=fd.getc();
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
            in=fd.getc();
        }
        if(!eol) {
            in=fd.getc();
        }
        out.push_back(field);
    }
}
std::atomic<size_t> assigned_count;
struct Assign_State {
    const std::vector<backward_pass_range>& child_idx_range;
    const std::vector<forward_pass_range>& parent_idx;
    FS_result_per_thread_t &output;
    void operator()(const Parsed_VCF_Line* vcf_line)const {
        auto& this_out=output.local();
        this_out.init(child_idx_range.size());
        assert(vcf_line->mutation.get_position()>0);
        Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx,vcf_line->mutation,vcf_line->mutated,this_out);
        assigned_count.fetch_add(1,std::memory_order_relaxed);
        delete vcf_line;
    }
};
void print_progress(std::atomic<bool>* done,std::mutex* done_mutex) {
    while (true) {
        {
            std::unique_lock<std::mutex> lk(*done_mutex);
            progress_bar_cv.wait_for(lk,std::chrono::seconds(1));
            if (done->load()) {
                return;
            }
        }
        fprintf(stderr,"\rAssigned %zu locus",assigned_count.load(std::memory_order_relaxed));
    }
}
template<typename infile_t>
static line_start_later try_get_first_line(infile_t& f,size_t& size ) {
    std::string temp;
    char c;
    while ((c=f.getc())!='\n') {
        temp.push_back(c);
    }
    temp.push_back('\n');
    size=temp.size()+1;
    char* buf=new char[size];
    strcpy(buf, temp.data());
    return line_start_later{buf,buf};
}
#define CHUNK_SIZ 200ul
#define ONE_GB 0x4ffffffful
template<typename infile_t>
static void process(MAT::Tree& tree,infile_t& fd) {
    std::vector<std::string> fields;
    read_header(fd, fields);
    std::unordered_set<std::string> inserted_samples;
    std::vector<MAT::Node*> bfs_ordered_nodes=tree.breadth_first_expansion();
    std::vector<long> idx_map(9);
    for (size_t idx=9; idx<fields.size(); idx++) {
        auto iter=tree.all_nodes.find(fields[idx]);
        if (iter==tree.all_nodes.end()) {
            fprintf(stderr, "sample %s cannot be found in tree\n",fields[idx].c_str());
            idx_map.push_back(-1);
            //exit(EXIT_FAILURE);
        } else {
            auto res=inserted_samples.insert(fields[idx]);
            if (res.second) {
                idx_map.push_back(iter->second->bfs_index);
            } else {
                idx_map.push_back(-1);
            }

        }
    }
    tbb::flow::graph input_graph;
    line_parser_t parser(input_graph,tbb::flow::unlimited,line_parser{idx_map});
    //feed used buffer back to decompressor

    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    FS_result_per_thread_t output;
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    tbb::flow::function_node<Parsed_VCF_Line*> assign_state(input_graph,tbb::flow::unlimited,Assign_State{child_idx_range,parent_idx,output});
    tbb::flow::make_edge(tbb::flow::output_port<0>(parser),assign_state);
    size_t single_line_size;
    parser.try_put(try_get_first_line(fd, single_line_size));
    size_t first_approx_size=std::min(CHUNK_SIZ,ONE_GB/single_line_size)-2;
    read_size=first_approx_size*single_line_size;
    alloc_size=(first_approx_size+2)*single_line_size;
    tbb::concurrent_bounded_queue<std::pair<char*,uint8_t*>> queue;
    queue.set_capacity(10);
    tbb::flow::source_node<line_start_later> line(input_graph,line_align(queue));
    tbb::flow::make_edge(line,parser);
    fd(queue);
    input_graph.wait_for_all();
    deallocate_FS_cache(output);
    fill_muts(output, bfs_ordered_nodes);
    size_t total_mutation_size=0;
    for(const auto node:bfs_ordered_nodes) {
        total_mutation_size+=node->mutations.size();
    }
    fprintf(stderr,"Total mutation size %zu \n", total_mutation_size);
    fd.unalloc();
}
void VCF_input(const char * name,MAT::Tree& tree) {
    assigned_count=0;
    std::atomic<bool> done(false);
    std::mutex done_mutex;
    std::thread progress_meter(print_progress,&done,&done_mutex);
    //open file set increase buffer size
    std::string vcf_filename(name);
    if (vcf_filename.find(".gz\0") != std::string::npos) {
        gzip_input_source fd(name);
        process(tree, fd);
    } else {
        raw_input_source fd(name);
        process(tree, fd);
    }
    done=true;
    progress_bar_cv.notify_all();
    progress_meter.join();
}
