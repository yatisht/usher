#include "detailed_mutation_load_store.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <mpi.h>
#include <sys/mman.h>
#include <tbb/flow_graph.h>
#include <thread>
#include <unistd.h>
#include <utility>
namespace MAT = Mutation_Annotated_Tree;
// mapping between chromosome string and chromosome index
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
static std::pair<uint64_t, uint64_t>
load_meta(MAT::Tree *tree, const uint8_t *start, int length) {
    google::protobuf::io::CodedInputStream inputi(start, length);
    Mutation_Detailed::meta meta;
    meta.ParseFromCodedStream(&inputi);
    load_chrom_vector(meta);
    tree->curr_internal_node = meta.internal_nodes_count();
    return std::make_pair(meta.root_offset(), meta.root_length());
}
static MAT::Mutation get_mutation(Mutation_Detailed::node &to_load,
                                  size_t mut_idx) {
    MAT::Mutation mut;
    mut.position = to_load.mutation_positions(mut_idx);
    *((uint32_t *)(&mut) + 1) = to_load.mutation_other_fields(mut_idx);
    return mut;
}
struct output_mutation {
    MAT::Mutations_Collection &mut_out;
    output_mutation(MAT::Mutations_Collection &mut_out) : mut_out(mut_out) {}
    void operator()(MAT::Mutations_Collection::const_iterator &par_iter,
                    const MAT::Mutation &last_mut) {
        while (par_iter->position < last_mut.position) {
            mut_out.push_back(*par_iter);
            par_iter++;
        }
        if (par_iter->position == last_mut.position) {
            par_iter++;
        }
        mut_out.push_back(last_mut);
    }
    void push_back(const MAT::Mutation &mut) { mut_out.push_back(mut); }
    void reserve(size_t in) { mut_out.reserve(in); }
    void exhaust(MAT::Mutations_Collection::const_iterator &par_iter) {
        while (par_iter->position < INT_MAX) {
            mut_out.push_back(*par_iter);
            par_iter++;
        }
        mut_out.push_back(*par_iter);
    }
};
struct no_output_mutation {
    void operator()(MAT::Mutations_Collection::const_iterator &par_iter,
                    const MAT::Mutation &last_mut) {}
    void push_back(const MAT::Mutation &mut) {}
    void exhaust(MAT::Mutations_Collection::const_iterator &par_iter) {}
    void reserve(size_t in) {}
};
template <typename outputer_type>
void fill_ignored_mutations(
    MAT::Node *&node, MAT::Mutations_Collection::const_iterator &par_iter,
    Mutation_Annotated_Tree::ignored_t::const_iterator &ignored_iter,
    outputer_type next_level_muts) {
    for (int position = ignored_iter->first; position <= ignored_iter->second;
         position++) {
        while (par_iter->position < position) {
            next_level_muts.push_back(*par_iter);
            par_iter++;
        }
        if (par_iter->position == position) {
            node->mutations.mutations.emplace_back(
                0, position, par_iter->get_mut_one_hot(),
                par_iter->get_mut_one_hot(), MAT::Mutation::ignored());
        } else {
            node->mutations.mutations.emplace_back(
                0, position, MAT::Mutation::refs[position],
                MAT::Mutation::refs[position], MAT::Mutation::ignored());
        }
    }
}

template <typename outputer_type>
static void load_mutations(MAT::Node *node, Mutation_Detailed::node &to_load,
                           size_t ignored_size,
                           const MAT::Mutations_Collection &parent_mutations,
                           outputer_type mut_out) {
    size_t valid_mut_size = to_load.mutation_positions_size();
    size_t mut_size = valid_mut_size + ignored_size;
    node->mutations.reserve(mut_size);
    mut_out.reserve(valid_mut_size + parent_mutations.size());
    auto par_iter = parent_mutations.begin();
    if (ignored_size == 0) {
        for (size_t mut_idx = 0; mut_idx < valid_mut_size; mut_idx++) {
            node->mutations.push_back(get_mutation(to_load, mut_idx));
            mut_out(par_iter, node->mutations.back());
        }
        mut_out.exhaust(par_iter);
        return;
    }
    MAT::ignored_t::const_iterator ignored_iter = node->ignore.begin();
    for (size_t mut_idx = 0; mut_idx < valid_mut_size; mut_idx++) {
        auto mut = get_mutation(to_load, mut_idx);
        while (ignored_iter->first < mut.position) {
            assert(ignored_iter->second < mut.position);
            fill_ignored_mutations(node, par_iter, ignored_iter, mut_out);
            ignored_iter++;
        }
        assert(ignored_iter->first > mut.position);
        node->mutations.push_back(mut);
        mut_out(par_iter, node->mutations.back());
    }
    while (ignored_iter->first < INT_MAX) {
        fill_ignored_mutations(node, par_iter, ignored_iter,
                               mut_out);
        ignored_iter++;
    }
    mut_out.exhaust(par_iter);
}
static void load_mutations(MAT::Node *node, Mutation_Detailed::node &to_load,
                           size_t ignored_size,
                           const MAT::Mutations_Collection &parent_mutations) {
    no_output_mutation foo;
    load_mutations(node, to_load, ignored_size, parent_mutations, foo);
}
static void load_mutations(MAT::Node *node, Mutation_Detailed::node &to_load,
                           size_t ignored_size,
                           const MAT::Mutations_Collection &parent_mutations,
                           MAT::Mutations_Collection &mutation_out) {
    output_mutation foo(mutation_out);
    load_mutations(node, to_load, ignored_size, parent_mutations, foo);
}
struct deserialize_condensed_nodes {
    void operator()(Mutation_Detailed::node node,
                    MAT::Tree::condensed_node_t &condensed_nodes,
                    const MAT::Node *out) {
        auto condensed_node_size = node.condensed_nodes_size();
        if (condensed_node_size) {
            std::vector<std::string> this_node_condensed;
            this_node_condensed.reserve(condensed_node_size);
            for (int i = 0; i < condensed_node_size; i++) {
                this_node_condensed.push_back(node.condensed_nodes(i));
            }
            condensed_nodes.emplace(out->identifier,
                                    std::move(this_node_condensed));
        }
    }
};
struct no_deserialize_condensed_nodes {
    void operator()(Mutation_Detailed::node node,
                    MAT::Tree::condensed_node_t &condensed_nodes,
                    const MAT::Node *out) {}
};
struct Load_Subtree_pararllel_Continuation : public tbb::task {
    MAT::Mutations_Collection mutation_so_far;
    tbb::task *execute() { return nullptr; }
};
template <typename do_serialize_condensed>
struct Load_Subtree_pararllel : public tbb::task {
    MAT::Node *parent;
    const uint8_t *file_start;
    int64_t start_offset;
    int length;
    MAT::Node *&out;
    MAT::Tree::condensed_node_t &condensed_nodes;
    const MAT::Mutations_Collection &parent_mutations;
    Load_Subtree_pararllel(MAT::Node *parent, const uint8_t *file_start,
                           int64_t start_offset, int length, MAT::Node *&out,
                           MAT::Tree::condensed_node_t &condensed_nodes,
                           const MAT::Mutations_Collection &parent_mutations)
        : parent(parent), file_start(file_start), start_offset(start_offset),
          length(length), out(out), condensed_nodes(condensed_nodes),
          parent_mutations(parent_mutations) {}
    tbb::task *execute() override {
        out = new MAT::Node();
        out->parent = parent;
        google::protobuf::io::CodedInputStream inputi(file_start + start_offset,
                                                      length);
        Mutation_Detailed::node node;
        node.ParseFromCodedStream(&inputi);
        out->identifier = node.identifier();
        // deserialize ignored range
        size_t ignore_range_size = node.ignored_range_end_size();
        int ignored_size = 0;
        if (ignore_range_size) {
        out->ignore.reserve(ignore_range_size+1);
        for (size_t ignore_idx = 0; ignore_idx < ignore_range_size;
             ignore_idx++) {
            auto start = node.ignored_range_start(ignore_idx);
            auto end = node.ignored_range_end(ignore_idx);
            out->ignore.emplace_back(start, end);
            ignored_size += (end - start)+1;
        }
        out->ignore.emplace_back(INT_MAX, INT_MAX);
        }
        // deserialize condensed nodes
        do_serialize_condensed()(node, condensed_nodes, out);
        size_t child_size = node.children_offsets_size();
        if (child_size) {
            out->children.resize(child_size);
            auto continuation = new (allocate_continuation())
                Load_Subtree_pararllel_Continuation;
            continuation->set_ref_count(child_size);
            load_mutations(out, node, ignored_size, parent_mutations,
                           continuation->mutation_so_far);
            for (size_t child_idx = 0; child_idx < child_size; child_idx++) {
                continuation->spawn(
                    *new (continuation->allocate_child())
                        Load_Subtree_pararllel(
                            out, file_start, node.children_offsets(child_idx),
                            node.children_lengths(child_idx),
                            out->children[child_idx], condensed_nodes,
                            continuation->mutation_so_far));
            }
            return nullptr;
        }
        load_mutations(out, node, ignored_size, parent_mutations);
        return new (allocate_continuation()) tbb::empty_task;
    }
};
struct no_free_input {
    void operator()(uint8_t *in) {}
};
struct free_input {
    void operator()(uint8_t *in) { free(in); }
};
template <typename do_free_input> struct decompressor_t {
    uint8_t *start;
    void operator()(uint8_t *to_load) {
        size_t offset = *(uint64_t*)to_load;
        size_t length = *(uint64_t*)(to_load + 8);
        unsigned long uncompressed_out = INT_MAX;
        auto result =
            uncompress(start + offset, &uncompressed_out, to_load + 16, length);
        if (result != Z_OK) {
            fprintf(stderr, "corrupted\n");
        }
        do_free_input()(to_load);
    }
};
struct file_loader {
    uint8_t *&start;
    uint8_t *end;
    bool operator()(uint8_t *&curr) {
        if (start==end) {
            return false;
        }
        curr=start;
        size_t length = *(uint64_t*)(start + 8);
        assert(start < end);
        start += (16 + length);
        return true;
    }
};
typedef tbb::flow::function_node<uint8_t *> decompressor_node_t;
static std::pair<uint8_t *, uint8_t *>
uncompress_file(const std::string &path) {
    int fd = open(path.c_str(), O_RDONLY);
    struct stat stat_buf;
    fstat(fd, &stat_buf);
    uint64_t file_size = stat_buf.st_size;
    uint8_t *file =
        (uint8_t *)mmap(nullptr, file_size, PROT_READ, MAP_SHARED, fd, 0);
    uint8_t *end_of_file = file + file_size;
    uint64_t uncompressed_size = *(uint64_t*)(end_of_file-8);
    uint8_t *uncompressed_buffer = (uint8_t *)malloc(uncompressed_size);
    tbb::flow::graph g;
    uint8_t *curr = file;
    decompressor_node_t decompressor(
        g, tbb::flow::unlimited,
        decompressor_t<no_free_input>{uncompressed_buffer});
    tbb::flow::source_node<uint8_t *> src(g, file_loader{curr, end_of_file-8});
    tbb::flow::make_edge(src, decompressor);
    g.wait_for_all();
    munmap(file, file_size);
    close(fd);
    return std::make_pair(uncompressed_buffer,
                          uncompressed_buffer + uncompressed_size);
}
static void receive_MPI(decompressor_node_t& decompressor,size_t& uncompressed_size){
    while (true) {
        size_t msg_size;
        //fprintf(stderr,"====Waiting segment length\n");
        MPI_Bcast(&msg_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        //fprintf(stderr,"====Receiving segmane tof length %zu \n", msg_size );
        if (msg_size==UINT64_MAX) {
            MPI_Bcast(&uncompressed_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
            return;
        }
        uint8_t* buffer=(uint8_t*) malloc(msg_size);
        MPI_Bcast(buffer, msg_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        //fprintf(stderr,"====Finished Receiving segmane tof length %zu \n", msg_size );
        decompressor.try_put(buffer);
    }
    //fprintf(stderr,"====Finished recieving \n" );
}
static std::pair<uint8_t *, uint8_t *>
receive_mpi_uncompress() {
    size_t max_memory;
    MPI_Bcast(&max_memory, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    //fprintf(stderr, "received memory requirement of %zu \n",max_memory);
    uint8_t *uncompressed_buffer = (uint8_t *)malloc(max_memory<<10);
    size_t uncompressed_size;
    tbb::flow::graph g;
    decompressor_node_t decompressor(
        g, tbb::flow::unlimited,
        decompressor_t<free_input>{uncompressed_buffer});
    receive_MPI(decompressor,uncompressed_size);
    g.wait_for_all();
    return std::make_pair(uncompressed_buffer,
                          uncompressed_buffer + uncompressed_size);
}

template<typename T>
static void deserialize_common(std::pair<uint8_t*,uint8_t*> uncompressed,MAT::Tree* tree){
    auto file = uncompressed.first;
    auto file_end = uncompressed.second;
    uint64_t meta_offset = *(uint64_t*)(file_end - 8);
    auto temp = load_meta(tree, file + meta_offset,
                          file_end - 8 - (file + meta_offset));
    MAT::Mutations_Collection root_muts;
    root_muts.mutations.emplace_back(INT_MAX);
    tbb::task::spawn_root_and_wait(
        *new (tbb::task::allocate_root())
            Load_Subtree_pararllel<T>(
                nullptr, (uint8_t *)file, temp.first, temp.second, tree->root,
                tree->condensed_nodes, root_muts));
    free(uncompressed.first);
    std::vector<MAT::Node *> dfs_nodes = tree->depth_first_expansion();
    for (auto node : dfs_nodes) {
        tree->all_nodes.emplace(node->identifier, node);
    }
}
// main load function
void Mutation_Annotated_Tree::Tree::load_detatiled_mutations(
    const std::string &path) {
    fputs("Loading intermediate protobuf\n", stderr);
    auto uncompressed = uncompress_file(path);
    deserialize_common<deserialize_condensed_nodes>(uncompressed, this);
    fputs("Finished loading intermediate protobuf\n", stderr);
}
// main load function
void Mutation_Annotated_Tree::Tree::MPI_receive_tree() {
    fputs("Loading intermediate protobuf\n", stderr);
    auto uncompressed = receive_mpi_uncompress();
    deserialize_common<no_deserialize_condensed_nodes>(uncompressed, this);
    fputs("Finished loading intermediate protobuf\n", stderr);
}