#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "zlib.h"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include <string>
#include <vector>
#include "import_vcf.hpp"
#include "Fitch_Sankoff.hpp"
#include "tbb/parallel_for.h"
#define ZLIB_BUFSIZ 0x10000
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
    void operator()(char* line_in, line_parser_t::output_ports_type& out)const{
        char* start=line_in;
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
            line_in++;
            //ID don't care
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //REF
            parsed_line=new Parsed_VCF_Line{MAT::Mutation(chromosome,pos, 0, 0, 0, 1,
                      MAT::get_nuc_id(*line_in))};
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

struct Assign_State{
    const std::vector<Mutation_Annotated_Tree::Node *>& bfs_ordered_nodes;
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> &output;
    void operator()(const Parsed_VCF_Line* vcf_line)const{
        Fitch_Sankoff_Whole_Tree(bfs_ordered_nodes,vcf_line->mutation,vcf_line->mutated,output);
        delete vcf_line;
    }
};

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

    std::vector<MAT::Node*> bfs_ordered_nodes=tree.breadth_first_expansion();
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> output(bfs_ordered_nodes.size());
    tbb::flow::function_node<Parsed_VCF_Line*> assign_state(input_graph,tbb::flow::unlimited,Assign_State{bfs_ordered_nodes,output});
    tbb::flow::make_edge(tbb::flow::output_port<0>(parser),assign_state);
    for (int i=0; i<1; i++) {
        decompressor.try_put((char*)malloc(52*header_size));
    }
    input_graph.wait_for_all();

    tbb::affinity_partitioner ap;
        tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto &to_refill = output[i];
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size(),false);
            }
        },
    ap);
}
