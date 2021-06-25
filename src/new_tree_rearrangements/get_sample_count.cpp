#include "parsimony.pb.h"
#include <cstddef>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <google/protobuf/io/coded_stream.h>
#include "../mutation_annotated_tree.hpp"
int count_leaves(Mutation_Annotated_Tree::Node* node){
    if(node->children.size()==0){
        return 1;
    }
    size_t count=0;
    for(auto n:node->children){
        count+=count_leaves(n);
    }
    return count;
}
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
    Mutation_Annotated_Tree::Tree tree=Mutation_Annotated_Tree::create_tree_from_newick_string(data.newick());
    int tree_leaves=count_leaves(tree.root);
    int num_condensed_nodes = data.condensed_nodes_size();
    for(int idx=0;idx<num_condensed_nodes;idx++){
        tree_leaves+=(data.condensed_nodes(idx).condensed_leaves_size()-1);
    }
    printf("%s : %d\n",argv[1],tree_leaves);
}