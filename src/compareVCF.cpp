#include "mutation_annotated_tree.hpp"
#include "zlib.h"
#include <array>
#include <cstdio>
#include <mutex>
#include <sstream>
#include <string>
#include <tbb/concurrent_hash_map.h>
#include <tbb/flow_graph.h>
#include <unordered_map>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
static int read_header(gzFile *fd, std::vector<std::string> &out) {
    int header_len = 0;
    char in = gzgetc(*fd);
    in = gzgetc(*fd);
    bool second_char_pong = (in == '#');

    while (second_char_pong) {
        while (in != '\n') {
            in = gzgetc(*fd);
        }
        in = gzgetc(*fd);
        in = gzgetc(*fd);
        second_char_pong = (in == '#');
    }

    bool eol = false;
    while (!eol) {
        std::string field;
        while (in != '\t') {
            if (in == '\n') {
                eol = true;
                break;
            }
            field.push_back(in);
            in = gzgetc(*fd);
            header_len++;
        }
        in = gzgetc(*fd);
        out.push_back(field);
    }
    return header_len;
}
#define ZLIB_BUFSIZ 0x40000
typedef tbb::flow::source_node<char *>
    decompressor_node_t;
struct Decompressor {
    gzFile *fd;
    size_t cont_read_size;
    size_t init_read_size;
    bool operator()(char *&buf) const {
        if (gzeof(*fd)) {
            return false;
        }
        buf=new char[init_read_size+cont_read_size];
        int read_size = gzread(*fd, buf, init_read_size);
        if (!read_size) {
            free(buf);
            return false;
        }
        if (!gzgets(*fd, buf + read_size, cont_read_size)) {
            *(buf + read_size) = 0;
        }
        return true;
    }
};
static std::array<std::string, 2> filenames;
static std::vector<std::string> samples;
struct VCFline {
    int file_idx;
    std::vector<int8_t> alleles;
};
static tbb::concurrent_hash_map<int, VCFline> parsed_lines;
std::mutex print_lock;
static void insert_line(int position,std::vector<int8_t>& variants,int file_idx){
    std::vector<int8_t> other;
    tbb::concurrent_hash_map<int, VCFline>::accessor accessor;
    parsed_lines.insert(accessor,position);
    if (accessor->second.alleles.empty()) {
        accessor->second.alleles.swap(variants);
        accessor->second.file_idx=file_idx;
        accessor.release();
    }else {
        assert(accessor->second.file_idx!=file_idx);
        other.swap(accessor->second.alleles);
        parsed_lines.erase(accessor);
        for (size_t i=0; i<other.size(); i++) {
            if (!(other[i]&variants[i])) {
                std::lock_guard<std::mutex> lock(print_lock);
                printf("At %d , sample %s , %c in file %s, %c in file %s\n",position,samples[i].c_str(),MAT::get_nuc(other[i]),filenames[!file_idx].c_str(),MAT::get_nuc(variants[i]),filenames[file_idx].c_str());
            }
        }
    }
}
struct line_parser {
    static const int CHROM_IDX = 0;
    static const int POS_IDX = 1;
    static const int REF_IDX = 3;
    static const int ALT_IDX = 4;
    static const int SAMPLE_START_IDX = 9;
    const std::vector<unsigned int> index_translate;
    int file_idx;
    void operator()(char *line_in) const {
        auto start=line_in;
        while (*line_in!=0) {
            std::vector<uint8_t> allele_translated;
            //std::string chromosome;
            int pos=0;
            //Chromosome
            while (*line_in!='\t') {
                //chromosome.push_back(*line_in);
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
            allele_translated.push_back(MAT::get_nuc_id(*line_in));
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
            std::vector<int8_t> out(index_translate.size());
            while (!is_last) {
                int allele_idx=(*line_in-'0');
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
                if (allele_idx>=(int)(allele_translated.size()+1)||allele_idx<0) {
                    out[index_translate[field_idx-SAMPLE_START_IDX]]=0xf;
                }else {
                    out[index_translate[field_idx-SAMPLE_START_IDX]]=allele_translated[allele_idx];
                }
                field_idx++;
                line_in++;
            }
            assert(field_idx==index_translate.size()+9);
            insert_line(pos, out, file_idx);
        }
        /*std::stringstream instream(line_in);
        while (instream) {
            std::string s;
            std::getline(instream, s);
            std::vector<std::string> words;
            MAT::string_split(s, words);

            if (words.size() != index_translate.size()+9) {
                fprintf(stderr, "ERROR! Incorrect VCF format.\n");
                exit(1);
            }

            std::vector<std::string> alleles;
            MAT::string_split(words[ALT_IDX], ',', alleles);
            std::vector<int8_t> allele_translated{MAT::get_nuc_id(words[REF_IDX][0])};
            allele_translated.reserve(alleles.size()+1);
            for (const auto &allele : alleles) {
                allele_translated.push_back(MAT::get_nuc_id(allele[0]));
            }
            std::vector<int8_t> out(index_translate.size());
            for (size_t j = SAMPLE_START_IDX; j < words.size(); j++) {
                if (isdigit(words[j][0])) {
                    int allele_id = std::stoi(words[j]);
                    out[index_translate[j-SAMPLE_START_IDX]]=allele_translated[allele_id];
                } else {
                    out[index_translate[j-SAMPLE_START_IDX]]=0xf;
                }
            }
            auto position=std::stoi(words[POS_IDX]);
            insert_line(position, out, file_idx);
        }*/
        delete[] start;
    }
};
static int  map_index(gzFile* fd1,gzFile* fd2,std::vector<unsigned int>& index_map1,std::vector<unsigned int>& index_map2){
    std::vector<std::string> header1;
    int header_size=read_header(fd1, header1);
    std::vector<std::string> header2;
    read_header(fd2, header2);

    std::unordered_map<std::string, unsigned int> idx_map;
    idx_map.reserve(2*header1.size());
    for (size_t i=0; i<header1.size()-line_parser::SAMPLE_START_IDX; i++) {
        samples.push_back(header1[i+line_parser::SAMPLE_START_IDX]);
        auto emplace_result=idx_map.emplace(header1[i+line_parser::SAMPLE_START_IDX],i);
        if (!emplace_result.second) {
            printf( "%s\t %d repeated\n",emplace_result.first->first.c_str(),emplace_result.first->second);
        }
        assert(emplace_result.second);
        index_map1.push_back(i);
    }
    for (size_t i=0; i<header2.size()-line_parser::SAMPLE_START_IDX; i++) {
        auto iter=idx_map.find(header2[i+line_parser::SAMPLE_START_IDX]);
        if(iter!=idx_map.end()){
            index_map2.push_back(iter->second);
            idx_map.erase(iter);
        }else {
            printf("sample %s missing\n",header2[i+line_parser::SAMPLE_START_IDX].c_str());
        }
    }
    for(auto s:idx_map){
            printf("sample %s missing\n",s.first.c_str());
    }
    return header_size;
}
#define LINE_PER_BLOCK 5
int main(int argc, char**argv){
    gzFile fd1=gzopen(argv[1], "r");
    filenames[0]=argv[1];
    gzbuffer(fd1,ZLIB_BUFSIZ);

    gzFile fd2=gzopen(argv[2], "r");
    filenames[1]=argv[2];
    gzbuffer(fd2,ZLIB_BUFSIZ);

    std::vector<unsigned int> index_map1;
    std::vector<unsigned int> index_map2;
    size_t header_size=map_index(&fd1, &fd2, index_map1, index_map2);

    tbb::flow::graph input_graph;
    decompressor_node_t decompressor1(input_graph,Decompressor{&fd1,LINE_PER_BLOCK*header_size,1*header_size});
    decompressor_node_t decompressor2(input_graph,Decompressor{&fd2,LINE_PER_BLOCK*header_size,1*header_size});

    tbb::flow::function_node<char*,tbb::flow::continue_msg> parser1(input_graph,tbb::flow::unlimited,line_parser{index_map1,0});
    tbb::flow::function_node<char*,tbb::flow::continue_msg> parser2(input_graph,tbb::flow::unlimited,line_parser{index_map2,1});
    tbb::flow::make_edge(decompressor1,parser1);

    tbb::flow::make_edge(decompressor2,parser2);
    input_graph.wait_for_all();
}

