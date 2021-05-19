#include "mutation_annotated_tree.hpp"
#include "zlib.h"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include <string>
#include <vector>
#include "tbb/parallel_for.h"
#define ZLIB_BUFSIZ 0x10000
#define PARSE_THREADS 10
namespace MAT = Mutation_Annotated_Tree;
struct Parsed_VCF_Line{
    Mutation_Annotated_Tree::Mutation mutation;
    std::unordered_map<std::string, uint8_t> mutated;
};
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<char*>> decompressor_node_t;
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<Parsed_VCF_Line*,char*>> line_parser_t;
struct Decompressor{
    gzFile* fd;
    size_t cont_read_size;
    size_t init_read_size;
    void operator()(char* buf,decompressor_node_t::output_ports_type& out) const{
        if (gzeof(*fd)) {
            free(buf);
            return;
        }
        int read_size=gzread(*fd, buf, init_read_size);
        if (!read_size) {
            free(buf);
            return;
        }
        if(!gzgets(*fd, buf+read_size, cont_read_size)){
            *(buf+read_size)=0;
        }
        std::get<0>(out).try_put(buf);
    }
};

struct line_parser{
    const std::vector<std::string>& header;
    //First output port is parsed line, and second output port returns the buffer to decompressor
    void operator()(char* line_in, line_parser_t::output_ports_type& out)const{
        char* start=line_in;
        Parsed_VCF_Line* parsed_line;
        while (*line_in!=0) {
            std::vector<uint8_t> allele_translated;
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
            line_in++;
            //ID don't care
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //REF
            MAT::Mutation temp;
            temp.chrom=chromosome;
            temp.position=pos;
            temp.ref_nuc=MAT::get_nuc_id(*line_in);
            parsed_line=new Parsed_VCF_Line{temp};
            line_in++;
            assert(*line_in=='\t');
            line_in++;
            //ALT
            while (*line_in!='\t') {
                allele_translated.push_back(MAT::get_nuc_id(*line_in));
                line_in++;
                if(*line_in==','){
                    line_in++;
                }else{
                    assert(*line_in=='\t');
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
                int allele_idx=(*line_in-'0');
                while (*line_in!='\t') {
                    if (*line_in=='\n'||*line_in==0) {
                        is_last=true;
                        break;
                    }
                    line_in++;
                }
                if (allele_idx>=(allele_translated.size()+1)||allele_idx<0) {
                    parsed_line->mutated.emplace(header[field_idx],0xf);
                }else if (allele_idx) {
                    parsed_line->mutated.emplace(header[field_idx],allele_translated[allele_idx-1]);
                }
                field_idx++;
                line_in++;
            }
            assert(field_idx==header.size());
            std::get<0>(out).try_put(parsed_line);
        }
        std::get<1>(out).try_put(start);
    }
};

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

void VCF_input(const char * name,MAT::Tree& tree){
    std::vector<std::string> fields;
    gzFile fd=gzopen(name, "r");
    gzbuffer(fd,ZLIB_BUFSIZ);
    unsigned int header_size=read_header(&fd, fields);
    tbb::flow::graph input_graph;
    decompressor_node_t decompressor(input_graph,1,Decompressor{&fd,50*header_size,2*header_size});
    line_parser_t parser(input_graph,tbb::flow::unlimited,line_parser{fields});
    tbb::flow::make_edge(tbb::flow::output_port<1>(parser),decompressor);
    tbb::flow::make_edge(tbb::flow::output_port<0>(decompressor),parser);
//TODO Connect tbb::flow::output_port<0>(parser) to consumer of parsed line

    for (int i=0; i<PARSE_THREADS; i++) {
        decompressor.try_put((char*)malloc(52*header_size));
    }
    input_graph.wait_for_all();

}
