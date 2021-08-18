#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "zlib.h"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>

#include <unordered_set>
#include <vector>
#include "import_vcf.hpp"
#include "tbb/parallel_for.h"
#define ZLIB_BUFSIZ 0x10000
//Decouple parsing (slow) and decompression, segment file into blocks for parallelized parsing
typedef tbb::flow::source_node<char*> decompressor_node_t;
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<Parsed_VCF_Line*>> line_parser_t;

struct Decompressor {
    gzFile* fd;
    size_t init_read_size;
    size_t cont_read_size;
    bool operator()(char*& buf) const {
        if (gzeof(*fd)) {
            return false;
        }
        buf=new char[init_read_size+cont_read_size];
        int read_size=gzread(*fd, buf, init_read_size);
        if (read_size<0) {
            int z_errnum = 0;
            fputs( gzerror(*fd,&z_errnum),stderr);
        }
        if (!read_size) {
            delete [] (buf);
            return false;
        }
        //Make sure the last line is complete in the block.
        if(!gzgets(*fd, buf+read_size, cont_read_size)) {
            *(buf+read_size)=0;
        }
        return true;
    }
};
//Parse a block of lines, assuming there is a complete line in the line_in buffer
struct line_parser {
    const std::vector<long>& header;
    void operator()(char* line_in, line_parser_t::output_ports_type& out)const {
        char* const start=line_in;
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
            assert(pos>0);
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
                    } else if (allele_idx) {
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
        delete[] (start);
    }
};
//tokenize header, get sample name
static int read_header(gzFile* fd,std::vector<std::string>& out) {
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
std::atomic<size_t> assigned_count;
struct Assign_State {
    const std::vector<backward_pass_range>& child_idx_range;
    const std::vector<forward_pass_range>& parent_idx;
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> &output;
    void operator()(const Parsed_VCF_Line* vcf_line)const {
        assert(vcf_line->mutation.get_position()>0);
        Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx,vcf_line->mutation,vcf_line->mutated,output);
        assigned_count.fetch_add(1,std::memory_order_relaxed);
        delete vcf_line;
    }
};
void print_progress(std::atomic<bool>* done,std::mutex* done_mutex) {
    while (true) {
        {
            std::unique_lock<std::mutex> lk(*done_mutex);
            if (timed_print_progress) {
                progress_bar_cv.wait_for(lk,std::chrono::seconds(1));

            } else {
                progress_bar_cv.wait(lk);
            }
            if (done->load()) {
                return;
            }
        }
        fprintf(stderr,"\rAssigned %zu locus",assigned_count.load(std::memory_order_relaxed));
    }
}
#define CHUNK_SIZ 10
void VCF_input(const char * name,MAT::Tree& tree) {
    assigned_count=0;
    std::vector<std::string> fields;
    //open file set increase buffer size
    gzFile fd=gzopen(name, "r");
    if (!fd) {
        fprintf(stderr, "cannnot open vcf file : %s, exiting.\n",name);
        exit(EXIT_FAILURE);
    }
    gzbuffer(fd,ZLIB_BUFSIZ);

    unsigned int header_size=read_header(&fd, fields);
    tbb::flow::graph input_graph;
    std::atomic<bool> done(false);
    std::mutex done_mutex;
    std::thread progress_meter(print_progress,&done,&done_mutex);
    decompressor_node_t decompressor(input_graph,Decompressor{&fd,CHUNK_SIZ*header_size,2*header_size});
    std::vector<long> idx_map(9);
    std::vector<MAT::Node*> bfs_ordered_nodes=tree.breadth_first_expansion();
    std::unordered_set<std::string> inserted_samples;
    for (size_t idx=9; idx<fields.size(); idx++) {
        auto iter=tree.all_nodes.find(fields[idx]);
        if (iter==tree.all_nodes.end()) {
            fprintf(stderr, "sample %s cannot be found\n",fields[idx].c_str());
            exit(EXIT_FAILURE);
        } else {
            auto res=inserted_samples.insert(fields[idx]);
            if (res.second) {
                idx_map.push_back(iter->second->bfs_index);
            } else {
                idx_map.push_back(-1);
            }

        }
    }
    line_parser_t parser(input_graph,tbb::flow::unlimited,line_parser{idx_map});
    //feed used buffer back to decompressor
    tbb::flow::make_edge(decompressor,parser);

    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> output(bfs_ordered_nodes.size());
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    tbb::flow::function_node<Parsed_VCF_Line*> assign_state(input_graph,tbb::flow::unlimited,Assign_State{child_idx_range,parent_idx,output});
    tbb::flow::make_edge(tbb::flow::output_port<0>(parser),assign_state);
    input_graph.wait_for_all();
    gzclose(fd);
    done=true;
    progress_bar_cv.notify_all();
    //Filling mutation vector
    tbb::affinity_partitioner ap;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
    [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            const auto &to_refill = output[i];
            bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                         to_refill.size());
        }
    },
    ap);
    size_t total_mutation_size=0;
    for(const auto node:bfs_ordered_nodes) {
        total_mutation_size+=node->mutations.size();
    }
    fprintf(stderr,"Total mutation size %zu \n", total_mutation_size);
    progress_meter.join();
}
