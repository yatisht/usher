#include "parsimony.pb.h"
#include <bits/types/FILE.h>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <google/protobuf/io/coded_stream.h>
#include <unordered_map>
#include <unordered_set>
#include "../mutation_annotated_tree.hpp"

int main(int argc, char** argv){
    Parsimony::data data;
    struct stat stat_buf;
    stat(argv[1],&stat_buf);
    size_t file_size=stat_buf.st_size;
    auto fd=open(argv[1], O_RDONLY);
    uint8_t* maped_file=(uint8_t*)mmap(nullptr, file_size, PROT_READ, MAP_SHARED,fd , 0);
    close(fd);
    google::protobuf::io::CodedInputStream input(maped_file,file_size);
    input.SetTotalBytesLimit(file_size*4, file_size*4);
    data.ParseFromCodedStream(&input);
    munmap(maped_file, file_size);
    std::ofstream tree_out(argv[2]);
    tree_out<<data.newick()<<"\n";
    tree_out.close();
    int num_condensed_nodes = data.condensed_nodes_size();
    struct condensed_remap_t {
        std::string mapped_name;
        bool is_set;
        condensed_remap_t(const std::string& name):mapped_name(name),is_set(false){}
        condensed_remap_t():is_set(false){}
    } ;
    std::unordered_map<std::string, std::string> condensed_map;
    std::unordered_set<std::string> written_condensed;
    condensed_map.reserve(2*num_condensed_nodes);
    written_condensed.reserve(2*num_condensed_nodes);
    for (int i=0; i<num_condensed_nodes; i++) {
        auto cn = data.condensed_nodes(i);
        for (int k = 0; k < cn.condensed_leaves_size(); k++) {
            if (cn.node_name()=="node_1938_condensed_3_leaves") {
                puts(cn.condensed_leaves(k).c_str());
            }
           condensed_map.emplace(cn.condensed_leaves(k),cn.node_name());
        }
    }
    FILE* fasta_in=fopen(argv[3], "r");
    FILE* fasta_out=fopen(argv[4], "w");
    perror("");
    int last_got_char=fgetc(fasta_in);
    bool write=false;
    char* seqname=(char*)malloc(BUFSIZ);
    size_t seq_len=BUFSIZ;
    fprintf(stdout, "%zu \n",condensed_map.size());
    puts(condensed_map.begin()->first.c_str());
    while (last_got_char!=EOF) {
        while (last_got_char!='>'&&last_got_char!=EOF) {
            if (write) {
                fputc(last_got_char, fasta_out);
            }
            last_got_char=fgetc(fasta_in);
        }
        if (last_got_char=='>') {
            auto len=getline(&seqname, &seq_len, fasta_in);
            fprintf(stdout, "%zu \n",len);
            fputs(seqname,stdout);
            std::string to_find(seqname);
            to_find.pop_back();
            auto iter=condensed_map.find(to_find);
            if (iter==condensed_map.end()) {
                fputc('>', fasta_out);
                fputs(seqname, fasta_out);
                write=true;
            }else {
                if (written_condensed.count(iter->second)) {
                    puts("skiped\n");
                    write=false;
                }else {
                    written_condensed.emplace(iter->second);
                    fputc('>', fasta_out);
                    fputs(iter->second.c_str(), fasta_out);
                    fputs("\n",fasta_out);
                    fputs(iter->second.c_str(), stdout);
                    write=true;
                }
            }
            puts("");
            last_got_char=fgetc(fasta_in);
        }
    }
    fclose(fasta_out);
}