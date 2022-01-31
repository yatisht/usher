#include "detailed_mutation_load_store.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tbb/flow_graph.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
/*
File structure:
offset of metadata message (8 byte)
repeated node messages (node in mutation_detailed.proto), one message for each
node in the tree, each contain the offset and size of its children metadata
message (meta in mutation_detailed.proto), contain offset and size of root node
*/
namespace MAT = Mutation_Annotated_Tree;
size_t get_memory() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}
static void save_chrom_vector(Mutation_Detailed::meta &to_save) {
    for (const std::string &chrom : MAT::Mutation::chromosomes) {
        to_save.add_chromosomes(chrom);
    }
    for (nuc_one_hot nuc : MAT::Mutation::refs) {
        to_save.add_ref_nuc(nuc.get_nuc_no_check());
    }
}
struct to_compress_container {
    std::vector<std::string *> to_compress;
    size_t curr_size;
    const size_t offset;
    to_compress_container(size_t offset) : offset(offset) {
        to_compress.reserve(0x10000);
        curr_size = 0;
    }
    void push_back(std::string *in) {
        curr_size += in->size();
        to_compress.push_back(in);
    }
    bool exceed_limit() {
        return curr_size >= BLOCK_SIZE;
    }
    operator bool() {
        return curr_size != 0;
    }
};
struct compressor_out {
    uint8_t *content;
    size_t size;
    compressor_out() {}
    compressor_out(size_t reserve_size) {
        content = (uint8_t *)malloc(reserve_size);
    }
    operator uint8_t *() {
        return content;
    }
};
struct compressor_node {
    compressor_out operator()(to_compress_container *to_compress) {
        size_t out_size = compressBound(to_compress->curr_size)+16;
        compressor_out output(out_size);
        uint8_t *input = (uint8_t *)malloc(to_compress->curr_size);
        int offset = 0;
        for (auto ele : to_compress->to_compress) {
            memcpy(input + offset, ele->data(), ele->size());
            offset += ele->size();
            delete ele;
        }
        compress(output+16, &out_size, input, offset);
        *(uint64_t*)output.content=to_compress->offset;
        *((uint64_t*)(output.content+8))=out_size;
        free(input);
        output.content=(uint8_t*)realloc(output.content,out_size+16);
        output.size=out_size+16;
        delete to_compress;
        return output;
    }
};

typedef tbb::flow::function_node<to_compress_container *, compressor_out>
compressor_node_t;
struct serializer_t {
    std::mutex out_mutex;
    size_t curr_offset;
    to_compress_container *to_compress;
    compressor_node_t &out_stream;
    serializer_t(compressor_node_t &out_stream) : out_stream(out_stream) {
        curr_offset = 0;
        to_compress = new to_compress_container(curr_offset);
    }
    void operator()(std::string *in, u_int64_t &offset_out) {
        std::lock_guard<std::mutex> lk(out_mutex);
        offset_out = curr_offset;
        curr_offset += in->size();
        to_compress->push_back(in);
        if (to_compress->exceed_limit()) {
            out_stream.try_put(to_compress);
            to_compress = new to_compress_container(curr_offset);
        }
    }
    size_t finalize() {
        if (to_compress) {
            out_stream.try_put(to_compress);
        } else {
            delete to_compress;
        }
        return curr_offset;
    }
};
struct serialize_condensed_node {
    void operator()(const MAT::Node *root, const MAT::Tree &tree,
                    Mutation_Detailed::node &this_node) {
        auto condensed_node_iter = tree.condensed_nodes.find(root->node_id);
        if (condensed_node_iter != tree.condensed_nodes.end()) {
            for (const auto &child_id : condensed_node_iter->second) {
                this_node.add_condensed_nodes(child_id);
            }
        }
    }
};
struct no_serialize_condensed_node {
    void operator()(const MAT::Node *root, const MAT::Tree &tree,
                    Mutation_Detailed::node &this_node) {}
};
static void serialize_node(Mutation_Detailed::node &this_node,
                           const MAT::Node *root, const MAT::Tree &tree) {
    this_node.set_node_id(root->node_id);
    this_node.set_changed(root->changed);
    this_node.mutable_mutation_positions()->Reserve(root->mutations.size());
    this_node.mutable_mutation_other_fields()->Reserve(root->mutations.size());
    for (const auto &mut : root->mutations) {
        if (mut.get_all_major_allele() != 0xf) {
            this_node.add_mutation_positions(mut.get_position());
            this_node.add_mutation_other_fields(*((uint32_t *)(&mut) + 1));
        }
    }
    if (root->ignore.size()) {
        this_node.mutable_ignored_range_start()->Reserve(root->ignore.size());
        this_node.mutable_ignored_range_end()->Reserve(root->ignore.size());
        for (const auto &ignored_one : root->ignore) {
            if (ignored_one.first==INT_MAX) {
                break;
            }
            this_node.add_ignored_range_start(ignored_one.first);
            this_node.add_ignored_range_end(ignored_one.second);
        }
    }
}
static void serialize_node(Mutation_Detailed::node &this_node,
                           u_int64_t &offset_out, u_int64_t &length_out,
                           serializer_t &out_stream) {
    std::string * out_string = new std::string;
    this_node.SerializeToString(out_string);
    length_out = out_string->size();
    out_stream(out_string,offset_out);
}
template <typename condensed_node_serializer>
static void serialize_node(const MAT::Node *root, const MAT::Tree &tree,
                           u_int64_t &offset_out, u_int64_t &length_out,
                           serializer_t &out_stream,
                           std::vector<u_int64_t> &children_offset,
                           std::vector<u_int64_t> &children_length) {
    Mutation_Detailed::node this_node;
    serialize_node(this_node, root, tree);
    this_node.mutable_children_lengths()->Reserve(children_offset.size());
    this_node.mutable_children_offsets()->Reserve(children_offset.size());
    for (size_t child_idx = 0; child_idx < children_offset.size();
            child_idx++) {
        this_node.add_children_lengths(children_length[child_idx]);
        this_node.add_children_offsets(children_offset[child_idx]);
    }
    condensed_node_serializer()(root, tree, this_node);
    serialize_node(this_node, offset_out, length_out, out_stream);
}
template <typename condensed_node_serializer>
static void serialize_node(const MAT::Node *root, const MAT::Tree &tree,
                           u_int64_t &offset_out, u_int64_t &length_out,
                           serializer_t &out_stream) {
    Mutation_Detailed::node this_node;
    serialize_node(this_node, root, tree);
    condensed_node_serializer()(root, tree, this_node);
    serialize_node(this_node, offset_out, length_out, out_stream);
}
template <typename condensed_node_serializer>
struct subtree_serializer_continuation : public tbb::task {
    std::vector<u_int64_t> children_offset;
    std::vector<u_int64_t> children_length;
    u_int64_t &offset_out;
    u_int64_t &length_out;
    const MAT::Node *root;
    const MAT::Tree &tree;
    serializer_t &out_stream;
    subtree_serializer_continuation(u_int64_t &offset_out,
                                    u_int64_t &length_out,
                                    const MAT::Node *root,
                                    const MAT::Tree &tree,
                                    serializer_t &out_stream,size_t children_size)
        : children_offset(children_size),children_length(children_size),offset_out(offset_out), length_out(length_out), root(root),
          tree(tree), out_stream(out_stream) {}
    tbb::task *execute() {
        serialize_node<condensed_node_serializer>(
            root, tree, offset_out, length_out, out_stream, children_offset,
            children_length);
        return nullptr;
    }
};
template <typename condensed_node_serializer>
struct subtree_serializer : public tbb::task {
    u_int64_t &offset_out;
    u_int64_t &length_out;
    const MAT::Node *root;
    const MAT::Tree &tree;
    serializer_t &out_stream;
    subtree_serializer(u_int64_t &offset_out, u_int64_t &length_out,
                       const MAT::Node *root, const MAT::Tree &tree,
                       serializer_t &out_stream)
        : offset_out(offset_out), length_out(length_out), root(root),
          tree(tree), out_stream(out_stream) {}
    tbb::task *execute() {
        auto child_size = root->children.size();
        if (child_size) {
            auto continuation = new (allocate_continuation())
            subtree_serializer_continuation<condensed_node_serializer>(
                offset_out, length_out, root, tree, out_stream,child_size);
            std::vector<std::pair<u_int64_t, u_int64_t>> child_addr(child_size);
            continuation->set_ref_count(child_size);
            for (size_t child_idx = 0; child_idx < child_size; child_idx++) {
                continuation->spawn(
                    *(new (continuation->allocate_child())
                      subtree_serializer<condensed_node_serializer>(
                          continuation->children_offset[child_idx],
                          continuation->children_length[child_idx],
                          root->children[child_idx], tree, out_stream)));
            }
            return nullptr;
        } else {
            serialize_node<condensed_node_serializer>(root, tree, offset_out,
                    length_out, out_stream);
            return new (allocate_continuation()) tbb::empty_task();
        }
    }
};
// meta mesaage
u_int64_t save_meta(const MAT::Tree &tree, u_int64_t root_offset,
                    u_int64_t root_length, serializer_t &serializer) {
    Mutation_Detailed::meta meta;
    meta.set_nodes_idx_next(tree.node_idx);
    auto node_names=meta.mutable_node_idx_map();
    node_names->Reserve(tree.node_names.size());
    for (const auto& name_map: tree.node_names) {
        auto next_idx=meta.add_node_idx_map();
        next_idx->set_node_id(name_map.first);
        next_idx->set_node_name(name_map.second);
    }
    save_chrom_vector(meta);
    meta.set_root_offset(root_offset);
    meta.set_root_length(root_length);
    auto serialized = new std::string;
    meta.SerializeToString(serialized);
    u_int64_t meta_offset;
    serializer(serialized, meta_offset);
    u_int64_t ignored;
    serializer(new std::string((char *)&meta_offset, 8), ignored);
    return serializer.finalize();
}
struct file_writer_t {
    const int fd;
    file_writer_t(int fd) : fd(fd) {}
    void operator()(compressor_out to_write) {
        auto size_written=write(fd, to_write.content, to_write.size);
        if (size_written!=(long) to_write.size) {
            puts("Failed to write intermediate protobuf");
        }
        free(to_write.content);
    }
};
struct mpi_writer {
    void operator()(compressor_out to_write) {
        MPI_Bcast( &to_write.size,1,MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        //fprintf(stderr,":::::Sending segmane tof length %zu\n", to_write.size );
        MPI_Bcast(to_write.content, to_write.size, MPI_BYTE, 0, MPI_COMM_WORLD);
        free(to_write.content);
    }
};
template<typename T>
size_t serialize_tree_general(const MAT::Tree* tree,serializer_t& serializer) {
    u_int64_t root_offset;
    u_int64_t root_length;
    tbb::task::spawn_root_and_wait(
        *(new (tbb::task::allocate_root())
          subtree_serializer<T>(
              root_offset, root_length, tree->root, *tree, serializer)));
    return save_meta(*tree, root_offset, root_length, serializer);
}
// main save function
void Mutation_Annotated_Tree::Tree::save_detailed_mutations(
    const std::string &path) const {
    tbb::flow::graph g;
    compressor_node_t compressor(g, tbb::flow::unlimited, compressor_node());
    auto fd = open(path.c_str(), O_CREAT | O_WRONLY | O_TRUNC,
                   S_IRUSR | S_IRGRP | S_IROTH);
    perror("");
    tbb::flow::function_node<compressor_out,tbb::flow::continue_msg,tbb::flow::queueing> file_writer(
        g, 1, file_writer_t(fd));
    tbb::flow::make_edge(compressor, file_writer);
    serializer_t serializer(compressor);
    auto uncompressed_length=serialize_tree_general<serialize_condensed_node>(this,serializer);
    g.wait_for_all();
    if(write(fd, &uncompressed_length, 8)!=8) {
        puts("Failed to write intermediate protobuf");
    }
    close(fd);
}

void Mutation_Annotated_Tree::Tree::MPI_send_tree() const {
    size_t max_memory=get_memory();
    MPI_Bcast(&max_memory, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    //fprintf(stderr, "Sending memory requirement of %zu \n",max_memory);
    tbb::flow::graph g;
    compressor_node_t compressor(g, tbb::flow::unlimited, compressor_node());
    tbb::flow::function_node<compressor_out,tbb::flow::continue_msg,tbb::flow::queueing> sender(
        g, 1, mpi_writer());
    tbb::flow::make_edge(compressor, sender);
    serializer_t serializer(compressor);
    auto uncompressed_length=serialize_tree_general<no_serialize_condensed_node>(this,serializer);
    g.wait_for_all();
    size_t size_t_max=UINT64_MAX;
    MPI_Bcast(&size_t_max, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&uncompressed_length, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
}
