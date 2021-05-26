#include "mutation_detailed.pb.h"
#include <cstddef>
#include <cstdint>
#include <google/protobuf/generated_message_util.h>
#include <ios>
#include <sys/mman.h>
#include <tbb/task.h>
#include <utility>
#define LOAD
#include "google/protobuf/io/coded_stream.h"
#include "mutation_annotated_tree.hpp"
#include <fcntl.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
void load_chrom_vector(const Mutation_Detailed::meta &to_load) {
    size_t chrom_size = to_load.chromosomes_size();
    MAT::Mutation::chromosomes.reserve(chrom_size);
    for (size_t idx = 0; idx < chrom_size; idx++) {
        MAT::Mutation::chromosomes.push_back(to_load.chromosomes(idx));
    }
    size_t ref_nuc_size = to_load.ref_nuc_size();
    MAT::Mutation::refs.reserve(ref_nuc_size);
    for (size_t i = 0; i < ref_nuc_size; i++) {
        MAT::Mutation::refs.push_back(nuc_one_hot(to_load.ref_nuc(i), true));
    }
}
void save_chrom_vector(Mutation_Detailed::meta &to_save) {
    for (const std::string &chrom : MAT::Mutation::chromosomes) {
        to_save.add_chromosomes(chrom);
    }
    for (nuc_one_hot nuc : MAT::Mutation::refs) {
        to_save.add_ref_nuc(nuc.get_nuc_no_check());
    }
}
void set_protobuf(const MAT::Mutation &m,
                  Mutation_Detailed::detailed_mutation *out) {
    out->set_position((uint32_t)m.position);
    uint32_t temp = *((uint32_t *)(&m) + 1);
    out->set_other_fields(temp);
}
static MAT::Mutation load_mutation_from_protobuf(
    const Mutation_Detailed::detailed_mutation &to_load) {
    MAT::Mutation m;
    m.position = to_load.position();
    int temp = to_load.other_fields();
    *((uint32_t *)(&m) + 1) = temp;
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
u_int64_t save_meta(const MAT::Tree &tree, u_int64_t root_offset,
                    u_int64_t root_length, std::ofstream &outfile) {
    Mutation_Detailed::meta meta;
    meta.set_internal_nodes_count(tree.curr_internal_node);
    save_chrom_vector(meta);
    meta.set_root_offset(root_offset);
    meta.set_root_length(root_length);
    u_int64_t meta_offset = outfile.tellp();
    meta.SerializePartialToOstream(&outfile);
    return meta_offset;
}

static std::pair<u_int64_t, u_int64_t> write_subtree(MAT::Node *root,
                                                     std::ofstream &outfile) {
    std::vector<u_int64_t> children_offset;
    std::vector<u_int64_t> children_length;
    Mutation_Detailed::node this_node;
    for (auto child : root->children) {
        auto temp = write_subtree(child, outfile);
        this_node.add_children_offsets(temp.first);
        this_node.add_children_lengths(temp.second);
    }
    this_node.set_identifier(root->identifier);
    Mutation_Detailed::mutation_list *temp =
        new Mutation_Detailed::mutation_list;
    dump_to_protobuf(root->mutations, temp);
    this_node.set_allocated_node_mutations(temp);
    u_int64_t start = outfile.tellp();
    this_node.SerializeToOstream(&outfile);
    u_int64_t end = outfile.tellp();
    return std::make_pair(start, end - start);
}

void Mutation_Annotated_Tree::Tree::save_detailed_mutations(
    const std::string &path) const {
    std::ofstream outfile(path, std::ios::out | std::ios::binary);
    outfile.seekp(8, std::ios::cur);
    auto ret = write_subtree(this->root, outfile);
    u_int64_t meta_offset = save_meta(*this, ret.first, ret.second, outfile);
    outfile.seekp(0, std::ios::beg);
    outfile.write((char *)&meta_offset, 8);
}
static std::pair<uint64_t, uint64_t>
load_meta(MAT::Tree *tree, const uint8_t *start, int length) {
    google::protobuf::io::CodedInputStream inputi(start, length);
    Mutation_Detailed::meta meta;
    meta.ParseFromCodedStream(&inputi);
    load_chrom_vector(meta);
    tree->curr_internal_node = meta.internal_nodes_count();
    return std::make_pair(meta.root_offset(), meta.root_length());
}
struct Load_Subtree_pararllel : public tbb::task {
    MAT::Node *parent;
    const uint8_t *file_start;
    int64_t start_offset;
    int length;
    MAT::Node *&out;
    Load_Subtree_pararllel(MAT::Node *parent, const uint8_t *file_start,
                           int64_t start_offset, int length, MAT::Node *&out)
        : parent(parent), file_start(file_start), start_offset(start_offset),
          length(length), out(out) {}
    tbb::task *execute() override {
        out = new MAT::Node();
        out->parent = parent;
        google::protobuf::io::CodedInputStream inputi(file_start + start_offset,
                                                      length);
        Mutation_Detailed::node node;
        node.ParseFromCodedStream(&inputi);
        out->identifier = node.identifier();
        load_from_protobuf(out->mutations, &node.node_mutations());
        size_t child_size = node.children_offsets_size();
        out->children = std::vector<MAT::Node *>(child_size);
        tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
        empty->set_ref_count(child_size);
        for (size_t child_idx = 0; child_idx < child_size; child_idx++) {
            empty->spawn(*new (
                empty->allocate_child()) Load_Subtree_pararllel{
                out, file_start, node.children_offsets(child_idx),
                node.children_lengths(child_idx), out->children[child_idx]});
        }
        return child_size?nullptr:empty;
    }
};
static void load_subtree(MAT::Node *parent, const uint8_t *file_start,
                         int start_offset, int length, MAT::Node *&out) {
    out = new MAT::Node();
    out->parent = out;
    google::protobuf::io::CodedInputStream inputi(file_start + start_offset,
                                                  length);
    Mutation_Detailed::node node;
    node.ParseFromCodedStream(&inputi);
    out->identifier = node.identifier();
    load_from_protobuf(out->mutations, &node.node_mutations());
    size_t child_size = node.children_offsets_size();
    out->children = std::vector<MAT::Node *>(child_size);
    for (size_t child_idx = 0; child_idx < child_size; child_idx++) {
        load_subtree(out, file_start, node.children_offsets(child_idx),
                     node.children_lengths(child_idx),
                     out->children[child_idx]);
    }
}
void Mutation_Annotated_Tree::Tree::load_detatiled_mutations(
    const std::string &path) {
    int fd = open(path.c_str(), O_RDONLY);
    struct stat stat_buf;
    fstat(fd, &stat_buf);
    uint64_t file_size = stat_buf.st_size;
    void *file = mmap(nullptr, file_size, PROT_READ, MAP_SHARED, fd, 0);
    uint64_t meta_offset = *((u_int64_t *)file);
    auto temp =
        load_meta(this, (uint8_t *)file + meta_offset, file_size - meta_offset);
    //load_subtree(nullptr, (uint8_t *)file, temp.first, temp.second, this->root);
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Load_Subtree_pararllel(nullptr, (uint8_t *)file, temp.first, temp.second, this->root));
    std::vector<MAT::Node*> bfs_nodes=breadth_first_expansion();
    for (auto node : bfs_nodes) {
        all_nodes.emplace(node->identifier,node);
    }
}