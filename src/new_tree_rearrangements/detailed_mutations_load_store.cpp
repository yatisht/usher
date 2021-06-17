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
/*
File structure:
offset of metadata message (8 byte)
repeated node messages (node in mutation_detailed.proto), one message for each node in the tree, each contain the offset and size of its children
metadata message (meta in mutation_detailed.proto), contain offset and size of root node
*/
namespace MAT = Mutation_Annotated_Tree;
//mapping between chromosome string and chromosome index
static void load_chrom_vector(const Mutation_Detailed::meta &to_load) {
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
static void save_chrom_vector(Mutation_Detailed::meta &to_save) {
    for (const std::string &chrom : MAT::Mutation::chromosomes) {
        to_save.add_chromosomes(chrom);
    }
    for (nuc_one_hot nuc : MAT::Mutation::refs) {
        to_save.add_ref_nuc(nuc.get_nuc_no_check());
    }
}
//each mutation, serialization format is dependent on in-memory representation, as it write all the state information
//toghether as a int, to save space, as protobuf don't have char data type
static void serialize_mutation(const MAT::Mutation &m,
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
//meta mesaage
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
//serialize each node recursively, post-order traversl, so parent node know the offset of its children
static std::pair<u_int64_t, u_int64_t> write_subtree(MAT::Node *root,const MAT::Tree& tree,
                                                     std::ofstream &outfile) {
    std::vector<u_int64_t> children_offset;
    std::vector<u_int64_t> children_length;
    Mutation_Detailed::node this_node;
    for (auto child : root->children) {
        auto temp = write_subtree(child, tree,outfile);
        this_node.add_children_offsets(temp.first);
        this_node.add_children_lengths(temp.second);
    }
    this_node.set_identifier(root->identifier);
    for(const auto& mut:root->mutations){
        auto serialized_mut=this_node.add_node_mutations();
        serialize_mutation(mut, serialized_mut);
    }
    auto condensed_node_iter=tree.condensed_nodes.find(root->identifier);
    if(condensed_node_iter!=tree.condensed_nodes.end()){
        for(const auto& child_id:condensed_node_iter->second){
            this_node.add_condensed_nodes(child_id);
        }
    }
    u_int64_t start = outfile.tellp();
    this_node.SerializeToOstream(&outfile);
    u_int64_t end = outfile.tellp();
    return std::make_pair(start, end - start);
}
// main save function
void Mutation_Annotated_Tree::Tree::save_detailed_mutations(
    const std::string &path) const {
    std::ofstream outfile(path, std::ios::out | std::ios::binary);
    outfile.seekp(8, std::ios::cur);
    auto ret = write_subtree(this->root,*this, outfile);
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
    MAT::Tree::condensed_node_t& condensed_nodes;
    Load_Subtree_pararllel(MAT::Node *parent, const uint8_t *file_start,
                           int64_t start_offset, int length, MAT::Node *&out,MAT::Tree::condensed_node_t& condensed_nodes)
        : parent(parent), file_start(file_start), start_offset(start_offset),
          length(length), out(out),condensed_nodes(condensed_nodes) {}
    tbb::task *execute() override {
        out = new MAT::Node();
        out->parent = parent;
        google::protobuf::io::CodedInputStream inputi(file_start + start_offset,
                                                      length);
        Mutation_Detailed::node node;
        node.ParseFromCodedStream(&inputi);
        out->identifier = node.identifier();
        //deserialize mutations
        auto mut_size=node.node_mutations_size();
        out->mutations.reserve(mut_size);
        for (int i=0; i<mut_size; i++) {
            out->mutations.push_back(load_mutation_from_protobuf(node.node_mutations(i)));
        }
        //deserialize condensed nodes
        auto condensed_node_size=node.condensed_nodes_size();
        if (condensed_node_size) {
            std::vector<std::string> this_node_condensed;
            this_node_condensed.reserve(condensed_node_size);
            for (int i=0; i<condensed_node_size; i++) {
                this_node_condensed.push_back(node.condensed_nodes(i));
            }
            condensed_nodes.emplace(out->identifier,std::move(this_node_condensed));
        }
        //spawn a new task for deserializing each children
        size_t child_size = node.children_offsets_size();
        out->children = std::vector<MAT::Node *>(child_size);
        tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
        empty->set_ref_count(child_size);
        for (size_t child_idx = 0; child_idx < child_size; child_idx++) {
            empty->spawn(*new (
                empty->allocate_child()) Load_Subtree_pararllel{
                out, file_start, node.children_offsets(child_idx),
                node.children_lengths(child_idx), out->children[child_idx],condensed_nodes});
        }
        return child_size?nullptr:empty;
    }
};
//main load function
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
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Load_Subtree_pararllel(nullptr, (uint8_t *)file, temp.first, temp.second, this->root,this->condensed_nodes));
    std::vector<MAT::Node*> dfs_nodes=depth_first_expansion();
    for (auto node : dfs_nodes) {
        all_nodes.emplace(node->identifier,node);
    }
}