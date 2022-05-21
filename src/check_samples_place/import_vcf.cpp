#include "isa-l/igzip_lib.h"
#include "oneTBB-2019_U9/include/tbb/blocked_range.h"
#include "oneTBB-2019_U9/include/tbb/enumerable_thread_specific.h"
#include "oneTBB-2019_U9/include/tbb/flow_graph.h"
#include "oneTBB-2019_U9/include/tbb/parallel_for.h"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include "tbb/parallel_for.h"
#include "src/matOptimize/check_samples.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <mutex>
#include <signal.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
#define ZLIB_BUFSIZ 0x10000
size_t read_size;
size_t alloc_size;
std::mutex ref_lock;
// Decouple parsing (slow) and decompression, segment file into blocks for
// parallelized parsing
typedef tbb::flow::source_node<char *> decompressor_node_t;
std::condition_variable progress_bar_cv;
struct raw_input_source {
    FILE *fh;
    raw_input_source(const char *fname) { fh = fopen(fname, "r"); }
    int getc() { return fgetc(fh); }
    void operator()(tbb::concurrent_bounded_queue<std::pair<char *, uint8_t *>>
                        &out) const {
        while (!feof(fh)) {
            auto line_out = new char[alloc_size];
            auto bytes_read = fread(line_out, 1, read_size, fh);
            if (bytes_read == 0) {
                delete[] line_out;
                break;
            }
            out.emplace(line_out, (uint8_t *)line_out + bytes_read);
        }
        out.emplace(nullptr, nullptr);
        out.emplace(nullptr, nullptr);
    }
    void unalloc() { fclose(fh); }
};
struct line_start_later {
    char *start;
    char *alloc_start;
};
struct line_align {
    mutable line_start_later prev;
    mutable uint8_t *prev_end;
    tbb::concurrent_bounded_queue<std::pair<char *, uint8_t *>> &in;
    line_align(tbb::concurrent_bounded_queue<std::pair<char *, uint8_t *>> &out)
        : in(out) {
        prev.start = nullptr;
    }
    bool operator()(line_start_later &out) const {
        std::pair<char *, unsigned char *> line;
        in.pop(line);
        if (line.first == nullptr) {
            if (prev.start != nullptr) {
                out = prev;
                prev.start = nullptr;
                prev_end[0] = 0;
                return true;
            } else {
                return false;
            }
        }
        if (!prev.start) {
            prev.start = line.first;
            prev.alloc_start = line.first;
            prev_end = line.second;
            in.pop(line);
            if (line.first==nullptr) {
                out=prev;
                prev.start=nullptr;
                *prev_end=0;
                return true;
            }
        }
        auto start_ptr = strchr(line.first, '\n');
        if (*start_ptr != '\n') {
            raise(SIGTRAP);
        }
        start_ptr++;
        auto cpy_siz = start_ptr - line.first;
        memcpy(prev_end, line.first, cpy_siz);
        prev_end[cpy_siz] = 0;
        out = prev;
        prev.start = start_ptr;
        prev.alloc_start = line.first;
        prev_end = line.second;
        return true;
    }
};
struct gzip_input_source {
    unsigned char *map_start;
    // char* read_curr;
    size_t mapped_size;
    struct isal_gzip_header gz_hdr;
    mutable unsigned char *getc_buf;
    mutable unsigned char *get_c_ptr;
    struct inflate_state *state;
    gzip_input_source(const char *fname) {
        struct stat stat_buf;
        stat(fname, &stat_buf);
        mapped_size = stat_buf.st_size;
        auto fh = open(fname, O_RDONLY);
        map_start = (unsigned char *)mmap(nullptr, mapped_size, PROT_READ,
                                          MAP_SHARED, fh, 0);
        // madvise(map_start, mapped_size, MADV_SEQUENTIAL|MADV_WILLNEED);
        close(fh);
        state = new struct inflate_state;
        getc_buf = new unsigned char[BUFSIZ];
        get_c_ptr = 0;

        isal_inflate_init(state);
        state->next_in = map_start;
        state->avail_in = mapped_size;
        state->crc_flag = IGZIP_GZIP;
        isal_gzip_header_init(&gz_hdr);
        auto ret = isal_read_gzip_header(state, &gz_hdr);
        fprintf(stderr, "Header ret: %d\n", ret);
        decompress_to_buffer(getc_buf, BUFSIZ);
        get_c_ptr = getc_buf;
    }
    void unalloc() {
        munmap(map_start, mapped_size);
    }
    bool decompress_to_buffer(unsigned char *buffer, size_t buffer_size) const {
        if (!state->avail_in) {
            return false;
        }
        state->next_out = buffer;
        state->avail_out = buffer_size;
        auto out = ISAL_DECOMP_OK;
        // while
        // (out==ISAL_DECOMP_OK&&state->avail_out>0&&state->avail_in>0&&state->block_state!=ISAL_BLOCK_FINISH)
        // {
        out = isal_inflate(state);
        //}
        if (out != ISAL_DECOMP_OK) {
            fprintf(stderr, "decompress error %d\n", out);
        }
        // copied from
        // https://github.com/intel/isa-l/blob/78f5c31e66fedab78328e592c49eefe4c2a733df/programs/igzip_cli.c#L952
        while (state->avail_out > 0 && state->avail_in > 0 &&
               state->next_in[0] == 31) {
            // fprintf(stderr,"continuing\n");
            if (state->avail_in > 1 && state->next_in[1] != 139)
                break;
            isal_inflate_reset(state);
            // while
            // (out==ISAL_DECOMP_OK&&state->avail_out>0&&state->avail_in>0&&state->block_state!=ISAL_BLOCK_FINISH)
            // {
            out = isal_inflate(state);
            //}
            if (out != ISAL_DECOMP_OK) {
                fprintf(stderr, "restart error %d\n", out);
            }
        }
        return state->avail_out < buffer_size;
    }
    int getc() {
        if (get_c_ptr == state->next_out) {
            get_c_ptr = getc_buf;
            if (!decompress_to_buffer(getc_buf, BUFSIZ)) {
                return -1;
            }
        }
        return *(get_c_ptr++);
    }

    void operator()(tbb::concurrent_bounded_queue<std::pair<char *, uint8_t *>>
                        &out) const {
        char *line_out = new char[alloc_size];
        // if (getc_buf) {
        auto load_size = getc_buf + BUFSIZ - get_c_ptr;
        memcpy(line_out, get_c_ptr, load_size);
        decompress_to_buffer((unsigned char *)line_out + load_size,
                             read_size - load_size);
        delete[](getc_buf);
        getc_buf = nullptr;
        fprintf(stderr, "getc_buff_deallocated\n");
        out.emplace(line_out, state->next_out);
        //}
        line_out = new char[alloc_size];
        while (decompress_to_buffer((unsigned char *)line_out, read_size)) {
            out.emplace(line_out, state->next_out);
            line_out = new char[alloc_size];
        }
        out.emplace(nullptr, nullptr);
        out.emplace(nullptr, nullptr);
        delete[] line_out;
    }
};
typedef tbb::enumerable_thread_specific<
    std::vector<std::vector<MAT::Mutation>>>
    Sampled_Tree_Mutations_t;
// Parse a block of lines, assuming there is a complete line in the line_in
// buffer
typedef tbb::flow::function_node<line_start_later> line_parser_t;
struct line_parser {
    Sampled_Tree_Mutations_t &header;
    std::vector<size_t> &sample_idx;
    size_t sample_size;
    void operator()(line_start_later line_in_struct) const {
        char *line_in = line_in_struct.start;
        char *to_free = line_in_struct.alloc_start;
        auto &this_blk = header.local();
        this_blk.resize(sample_size);
        while (*line_in != 0) {
            std::vector<nuc_one_hot> allele_translated;
            std::string chromosome;
            int pos = 0;
            // Chromosome
            while (*line_in != '\t') {
                chromosome.push_back(*line_in);
                line_in++;
            }
            line_in++;
            // Position
            while (*line_in != '\t') {
                pos = pos * 10 + (*line_in - '0');
                line_in++;
            }
            if (pos <= 0) {
                raise(SIGTRAP);
            }
            line_in++;
            // ID don't care
            while (*line_in != '\t') {
                line_in++;
            }
            line_in++;
            // REF
            auto ref_nuc=MAT::get_nuc_id(*line_in);
            MAT::Mutation mut_template(chromosome,pos, 0, ref_nuc, 0,ref_nuc);
            line_in++;
            // assert(*line_in=='\t');
            line_in++;
            // ALT
            while (*line_in != '\t') {
                allele_translated.push_back(MAT::get_nuc_id(*line_in));
                line_in++;
                if (*line_in == ',') {
                    line_in++;
                } else {
                    // assert(*line_in=='\t');
                }
            }
            line_in++;
            unsigned int field_idx = 5;
            for (; field_idx < 9; field_idx++) {
                while (*line_in != '\t') {
                    line_in++;
                }
                line_in++;
            }
            // samples
            bool is_last = false;
            while (!is_last) {
                unsigned int allele_idx = (*line_in - '0');
                line_in++;
                while (std::isdigit(*line_in)) {
                    allele_idx *= 10;
                    allele_idx += (*line_in - '0');
                    line_in++;
                }
                while (*line_in != '\t') {
                    if (*line_in == '\n' || *line_in == 0) {
                        is_last = true;
                        break;
                    }
                    line_in++;
                }
                auto output_idx = sample_idx[field_idx];
                if (output_idx != SIZE_MAX) {
                    // output prototype of mutation, and a map from sample to
                    // non-ref allele
                    MAT::Mutation this_mut(mut_template);
                    if (allele_idx >= (allele_translated.size() + 1)) {
                        this_mut.set_mut_one_hot(0xf);
                        this_mut.set_auxillary(0xf, 0);
                        this_mut.set_descendant_mut(0xf);
                        this_blk[output_idx].push_back(this_mut);
                    } else if (allele_idx) {
                        this_mut.set_mut_one_hot(allele_translated[allele_idx - 1]);
                        this_mut.set_descendant_mut(allele_translated[allele_idx - 1]);
                        this_mut.set_auxillary(allele_translated[allele_idx - 1], 0);
                        assert(this_mut.get_par_one_hot()!=this_mut.get_mut_one_hot());
                        this_blk[output_idx].push_back(this_mut);
                    }
                }

                field_idx++;
                line_in++;
            }
        }
        delete[](to_free);
    }
};
// tokenize header, get sample name
template <typename infile_t>
static void read_header(infile_t &fd, std::vector<std::string> &out) {
    char in = fd.getc();
    in = fd.getc();
    bool second_char_pong = (in == '#');

    while (second_char_pong) {
        while (in != '\n') {
            in = fd.getc();
        }
        in = fd.getc();
        in = fd.getc();
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
            in = fd.getc();
        }
        if (!eol) {
            in = fd.getc();
        }
        out.push_back(field);
    }
}
std::atomic<size_t> assigned_count;
void print_progress(std::atomic<bool> *done, std::mutex *done_mutex) {
    while (true) {
        {
            std::unique_lock<std::mutex> lk(*done_mutex);
            progress_bar_cv.wait_for(lk, std::chrono::seconds(1));
            if (done->load()) {
                return;
            }
        }
        fprintf(stderr, "\rAssigned %zu locus",
                assigned_count.load(std::memory_order_relaxed));
    }
}
template <typename infile_t>
static line_start_later try_get_first_line(infile_t &f, size_t &size) {
    std::string temp;
    char c;
    while ((c = f.getc()) != '\n') {
        temp.push_back(c);
    }
    temp.push_back('\n');
    size = temp.size() + 1;
    char *buf = new char[size];
    strcpy(buf, temp.data());
    return line_start_later{buf, buf};
}
#define CHUNK_SIZ 200ul
#define ONE_GB 0x4ffffffful
template <typename infile_t>
static void process(infile_t &fd, Original_State_t& ori_state,
                    MAT::Tree &tree) {
    std::vector<std::string> fields;
    read_header(fd, fields);
    tbb::flow::graph input_graph;
    Sampled_Tree_Mutations_t tree_mutations;
    std::vector<size_t> sample_idx(fields.size(), SIZE_MAX);
    std::vector<std::string> sample_names;
    sample_names.reserve(fields.size());
    int curr_sample_idx = 0;
    for (size_t field_idx = 9; field_idx < fields.size(); field_idx++) {
        if (tree.get_node(fields[field_idx]) == nullptr) {
            sample_names.push_back(fields[field_idx]);
            sample_idx[field_idx] = curr_sample_idx;
            curr_sample_idx++;
        }else {
            //fprintf(stderr, "%s already present\n",fields[field_idx].c_str());
        }
    }
    line_parser_t parser(
        input_graph, tbb::flow::unlimited,
        line_parser{tree_mutations, sample_idx, sample_names.size()});
    // feed used buffer back to decompressor

    size_t single_line_size;
    parser.try_put(try_get_first_line(fd, single_line_size));
    size_t first_approx_size =
        std::min(CHUNK_SIZ, ONE_GB / single_line_size) - 2;
    read_size = first_approx_size * single_line_size;
    alloc_size = (first_approx_size + 2) * single_line_size;
    tbb::concurrent_bounded_queue<std::pair<char *, uint8_t *>> queue;
    tbb::flow::source_node<line_start_later> line(input_graph,
                                                  line_align(queue));
    tbb::flow::make_edge(line, parser);
    fd(queue);
    input_graph.wait_for_all();
    fd.unalloc();
    auto start_idx=tree.map_samp_name_only(sample_names[0]);
    for (size_t idx=1; idx<sample_names.size(); idx++) {
        tree.map_samp_name_only(sample_names[idx]);
    }
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, sample_names.size()),
        [&](tbb::blocked_range<size_t> range) {
            for (size_t idx = range.begin(); idx < range.end(); idx++) {
                size_t mut_count = 2;
                for (const auto &one_thread : tree_mutations) {
                    mut_count += one_thread[idx].size();
                }
                Mutation_Set this_out;
                //sample_mutations[idx].sample_name = ;
                //auto &this_out = sample_mutations[idx].muts;
                this_out.reserve(mut_count);
                for (auto &one_thread : tree_mutations) {
                    this_out.insert(one_thread[idx].begin(),
                                    one_thread[idx].end());
                }
                ori_state.emplace(start_idx+idx,std::move(this_out));
            }
        });
}
void Sample_Input(const char *name, Original_State_t& ori_state,
                  MAT::Tree &tree) {
    assigned_count = 0;
    std::atomic<bool> done(false);
    std::mutex done_mutex;
    std::thread progress_meter(print_progress, &done, &done_mutex);
    // open file set increase buffer size
    std::string vcf_filename(name);
    if (vcf_filename.find(".gz\0") != std::string::npos) {
        gzip_input_source fd(name);
        process(fd, ori_state, tree);
    } else {
        raw_input_source fd(name);
        process(fd, ori_state, tree);
    }
    done = true;
    progress_bar_cv.notify_all();
    progress_meter.join();
}
