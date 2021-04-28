#include "mutation_detailed.pb.h"
#include <google/protobuf/generated_message_util.h>
#include <sys/mman.h>
#define LOAD
#include "mutation_annotated_tree.hpp"
#include <fstream>
#include "google/protobuf/io/coded_stream.h"
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
namespace MAT = Mutation_Annotated_Tree;
void load_chrom_vector(const Mutation_Detailed::data &to_load) {
    size_t chrom_size = to_load.chromosomes_size();
    MAT::Mutation::chromosomes.reserve(chrom_size);
    for (size_t idx = 0; idx < chrom_size; idx++) {
        MAT::Mutation::chromosomes.push_back(to_load.chromosomes(idx));
    }
    size_t ref_nuc_size=to_load.ref_nuc_size();
    MAT::Mutation::refs.reserve(ref_nuc_size);
    for (size_t i=0; i<ref_nuc_size; i++) {
        MAT::Mutation::refs.push_back(nuc_one_hot(to_load.ref_nuc(i),true));
    }
}
void save_chrom_vector(Mutation_Detailed::data &to_save) {
    for (const std::string &chrom : MAT::Mutation::chromosomes) {
        to_save.add_chromosomes(chrom);
    }
    for(nuc_one_hot nuc:MAT::Mutation::refs){
        to_save.add_ref_nuc(nuc.get_nuc_no_check());
    }
}
void set_protobuf(const MAT::Mutation &m,
                  Mutation_Detailed::detailed_mutation *out) {
    out->set_position((uint32_t)m.position);
    out->set_chrom_idx((uint32_t)m.chrom_idx);
    out->set_par_mut_nuc((uint32_t)m.par_mut_nuc);
    out->set_boundary1_tie((uint32_t)m.boundary1_tie);
    out->set_boundary2_flag((uint32_t)m.boundary2_flag);
}
static MAT::Mutation load_mutation_from_protobuf(
    const Mutation_Detailed::detailed_mutation &to_load) {
    MAT::Mutation m;
    m.position = to_load.position();
    m.chrom_idx = to_load.chrom_idx();
    m.par_mut_nuc = to_load.par_mut_nuc();
    m.boundary1_tie = to_load.boundary1_tie();
    m.boundary2_flag = to_load.boundary2_flag();
    return m;
}

void dump_to_protobuf(const MAT::Mutations_Collection &to_dump,
                      Mutation_Detailed::mutation_list *out) {
    for (const MAT::Mutation &m : to_dump.mutations) {
        auto t = out->add_mutation();
        set_protobuf(m, t);
    }
}
void load_from_protobuf(MAT::Mutations_Collection &to_load,
                        const Mutation_Detailed::mutation_list *in) {
    std::vector<MAT::Mutation> &mutations = to_load.mutations;
    mutations.clear();
    auto mut_count = in->mutation_size();
    mutations.reserve(mut_count);
    for (int i = 0; i < mut_count; i++) {
        mutations.push_back(load_mutation_from_protobuf(in->mutation(i)));
    }
}

void Mutation_Annotated_Tree::Tree::save_detailed_mutations(
    const std::string &path) const {
    Mutation_Detailed::data data;
    data.set_newick(get_newick_string(*this, true, true));
    save_chrom_vector(data);
    data.set_internal_nodes_count(this->curr_internal_node);
    std::vector<Node *> dfs_ordered_nodes = depth_first_expansion();
    for (const Node *node : dfs_ordered_nodes) {
        auto mut_list = data.add_node_mutations();
        dump_to_protobuf(node->mutations, mut_list);
    }
    std::ofstream outfile(path, std::ios::out | std::ios::binary);
    data.SerializeToOstream(&outfile);
}
void Mutation_Annotated_Tree::Tree::load_detatiled_mutations(
    const std::string &path) {
    Mutation_Detailed::data data;
    int fd=open(path.c_str(),O_RDONLY);
    struct stat stat_buf;
    fstat(fd,&stat_buf);
    void* file=mmap(nullptr, stat_buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
    google::protobuf::io::CodedInputStream inputi((const uint8_t*) file,stat_buf.st_size);
    inputi.SetTotalBytesLimit(stat_buf.st_size*4, stat_buf.st_size);
    data.ParseFromCodedStream(&inputi);
    this->curr_internal_node = data.internal_nodes_count();
    load_from_newick(data.newick(), true);
    load_chrom_vector(data);
    std::vector<Node *> dfs_ordered_nodes = depth_first_expansion();
    for (int i = 0; i < (int)dfs_ordered_nodes.size(); i++) {
        const auto &mut_list = data.node_mutations(i);
        load_from_protobuf(dfs_ordered_nodes[i]->mutations, &mut_list);
    }
}