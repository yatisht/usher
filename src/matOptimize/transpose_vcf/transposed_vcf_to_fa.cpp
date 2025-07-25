#include "../mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "transpose_vcf.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_pipeline.h>
#include <tbb/global_control.h>
#include <tbb/info.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>
std::string chrom;
std::string ref;
struct Pos_Mut {
    int position;
    uint8_t mut;
    Pos_Mut(int position, uint8_t mut) : position(position), mut(mut) {}
    bool operator<(const Pos_Mut &other) const {
        return position < other.position;
    }
};

struct Sample_Pos_Mut {
    std::string name;
    std::vector<Pos_Mut> not_Ns;
    std::vector<std::pair<int, int>> Ns;
    Sample_Pos_Mut(std::string &&name) : name(std::move(name)) {}
};
struct Sample_Pos_Mut_Printer {
    std::vector<Pos_Mut>::const_iterator not_Ns_iter;
    const std::vector<Pos_Mut>::const_iterator not_Ns_end;
    std::vector<std::pair<int, int>>::const_iterator Ns_iter;
    const std::vector<std::pair<int, int>>::const_iterator Ns_end;
    Sample_Pos_Mut_Printer(const Sample_Pos_Mut &in)
        : not_Ns_iter(in.not_Ns.begin()), not_Ns_end(in.not_Ns.end()),
          Ns_iter(in.Ns.begin()), Ns_end(in.Ns.end()) {}
    int print(uint8_t *out) {
        int end;
        if (not_Ns_iter == not_Ns_end ||
                (Ns_iter != Ns_end && Ns_iter->first < not_Ns_iter->position)) {
            memset(out, 'N', Ns_iter->second - Ns_iter->first + 1);
            end = Ns_iter->second;
            Ns_iter++;
        } else {
            *out = Mutation_Annotated_Tree::get_nuc(not_Ns_iter->mut);
            end = not_Ns_iter->position;
            not_Ns_iter++;
        }
        assert(not_Ns_iter<=not_Ns_end&&Ns_iter<=Ns_end);
        return end;
    }
    int next_pos() const {
        auto next_pos=std::min(not_Ns_end == not_Ns_iter ? INT_MAX
                               : not_Ns_iter->position,
                               Ns_end == Ns_iter ? INT_MAX : Ns_iter->first);
        assert(next_pos!=INT_MAX);
        return (next_pos);
    }
    operator bool() const {
        return (Ns_iter != Ns_end) || (not_Ns_iter != not_Ns_end);
    }
};

struct Sample_Pos_Mut_Wrap {
    Sample_Pos_Mut &content;
    void add_Not_N(int position, uint8_t allele) {
        content.not_Ns.emplace_back(position-1, allele);
    }
    void add_N(int first, int second) {
        content.Ns.emplace_back(first-1, second-1);
    }
};

typedef tbb::enumerable_thread_specific<std::vector<Sample_Pos_Mut>>
        sample_pos_mut_local_t;
sample_pos_mut_local_t sample_pos_mut_local;
struct All_Sample_Appender {
    Sample_Pos_Mut_Wrap set_name(std::string &&name) {
        sample_pos_mut_local_t::reference my_sample_pos_mut_local =
            sample_pos_mut_local.local();
        my_sample_pos_mut_local.emplace_back(std::move(name));
        return Sample_Pos_Mut_Wrap{my_sample_pos_mut_local.back()};
    }
};
void load_reference(const char *ref_path, std::string &seq_name,
                    std::string &reference) {
    std::ifstream fasta_f(ref_path);
    std::string temp;
    std::getline(fasta_f, temp);
    seq_name = temp.substr(1, 50);
    while (fasta_f) {
        std::getline(fasta_f, temp);
        reference += temp;
    }
}
namespace po = boost::program_options;
size_t write_sample(uint8_t *out, const Sample_Pos_Mut &in,
                    const std::string &seq) {
    out[0] = '>';
    memcpy(out + 1, in.name.c_str(), in.name.size());
    out[1 + in.name.size()] = '\n';
    auto start_offset = 2 + in.name.size();
    int last_print = 0;
    Sample_Pos_Mut_Printer printer(in);
    while (printer) {
        auto next = printer.next_pos();
        assert(next<seq.size());
        auto ref_print_siz = next - last_print;
        memcpy(out + last_print + start_offset, seq.data() + last_print,
               ref_print_siz);
        last_print = printer.print(out + next + start_offset) + 1;
        assert(last_print<seq.size());
    }
    memcpy(out + last_print + start_offset, seq.data() + last_print,
           seq.size() - last_print);
    out[start_offset + seq.size()] = '\n';
    return start_offset + seq.size() + 1;
}
#define OUT_BUF_SIZ 0x100000
struct Batch_Printer {
    const std::vector<Sample_Pos_Mut> &in;
    const std::string &seq;
    const size_t inbuf_siz;
    const int fd;
    std::mutex& out_mutex;
    mutable z_stream stream;
    mutable uint8_t *inbuf;
    mutable uint8_t *outbuf;
    Batch_Printer(const std::vector<Sample_Pos_Mut> &in, const std::string &seq,
                  size_t inbuf_siz,int fd,std::mutex& out_mutex)
        : in(in), seq(seq), inbuf_siz(inbuf_siz),fd(fd),out_mutex(out_mutex), inbuf(nullptr),
          outbuf(nullptr) {
        stream.opaque = nullptr;
        stream.zalloc = nullptr;
        stream.zfree = nullptr;
    }
    Batch_Printer(const Batch_Printer &other)
        : in(other.in), seq(other.seq), inbuf_siz(other.inbuf_siz),fd(other.fd), out_mutex(other.out_mutex),inbuf(nullptr),
          outbuf(nullptr) {
        stream.opaque = nullptr;
        stream.zalloc = nullptr;
        stream.zfree = nullptr;
    }
    void init() const {
        if (!inbuf) {
            inbuf=(uint8_t*) malloc(inbuf_siz);
            outbuf=(uint8_t*) malloc(OUT_BUF_SIZ);
            //auto err =
            deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                         15 | 16, 9, Z_DEFAULT_STRATEGY);
            //assert(err==Z_OK);
        } else {
            //auto err=
            deflateReset(&stream);
            //assert(err==Z_OK);
        }
        stream.avail_out=OUT_BUF_SIZ;
        stream.next_out=outbuf;
    }
    void output()const {
        if (stream.avail_out<OUT_BUF_SIZ/2) {
            //auto err =
            deflate(&stream, Z_FINISH);
            //assert(err==Z_STREAM_END);
            out_mutex.lock();
            auto out=write(fd,outbuf,stream.next_out-outbuf);
            if (out!=stream.next_out-outbuf) {
                fprintf(stderr, "Couldn't output\n");
                exit(EXIT_FAILURE);
            }
            out_mutex.unlock();
            deflateReset(&stream);
            stream.avail_out=OUT_BUF_SIZ;
            stream.next_out=outbuf;
        }
    }
    void operator()(tbb::blocked_range<size_t> range) const {
        //fprintf(stderr, "%zu ",range.size());
        init();
        stream.avail_out = OUT_BUF_SIZ;
        stream.next_out = outbuf;
        for (size_t samp_idx = range.begin(); samp_idx < range.end();
                samp_idx++) {
            auto bytes_written=write_sample(inbuf, in[samp_idx], seq);
            stream.next_in=inbuf;
            stream.avail_in=bytes_written;
            auto err = deflate(&stream, Z_NO_FLUSH);
            assert(stream.avail_out);
            if (err!=Z_OK) {
                fprintf(stderr, "%d\n",err);
            }
            assert(err==Z_OK);
            output();
        }
        //auto err =
        deflate(&stream, Z_FINISH);
        //assert(err==Z_STREAM_END);
        out_mutex.lock();
        auto written_bytes=write(fd, outbuf, stream.next_out-outbuf);
        if (written_bytes!=stream.next_out-outbuf) {
            fprintf(stderr, "Couldn't output\n");
            exit(EXIT_FAILURE);
        }
        out_mutex.unlock();
    }
    ~Batch_Printer() {
        if (inbuf) {
            free(inbuf);
            free(outbuf);
            deflateEnd(&stream);
        }
    }
};
int main(int argc, char **argv) {
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::info::default_concurrency();
    std::string output_fa_file;
    std::string input_path;
    std::string rename_file;
    uint32_t num_threads;
    std::string reference;
    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_cores) +
                                      " detected on this machine]";
    desc.add_options()(
        "fa,o", po::value<std::string>(&output_fa_file)->required(),
        "Output FASTA file (in uncompressed or gzip-compressed .gz format) ")(
            "threads,T",
            po::value<uint32_t>(&num_threads)->default_value(num_cores),
            num_threads_message
            .c_str())("reference,r",
                      po::value<std::string>(&reference)->required(),
                      "Reference file")("samples,s",
                                        po::value<std::string>(&rename_file),
                                        "rename file")("output_path,i",
                                                po::value<std::string>(
                                                        &input_path)
                                                ->required(),
                                                "Load transposed VCF");
    po::options_description all_options;
    all_options.add(desc);
    po::variables_map vm;
    try {
        po::store(
            po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            return 0;
        else
            return 1;
    }
    std::unordered_map<std::string, std::string> rename_mapping;
    if (rename_file != "") {
        parse_rename_file(rename_file, rename_mapping);
    }
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads);
    load_reference(reference.c_str(), chrom, ref);
    All_Sample_Appender appender;
    load_mutations(input_path.c_str(), 80, appender);
    std::vector<Sample_Pos_Mut> all_samples;
    if (rename_mapping.empty()) {
        for (auto &sample_block : sample_pos_mut_local) {
            all_samples.insert(all_samples.end(),
                               std::make_move_iterator(sample_block.begin()),
                               std::make_move_iterator(sample_block.end()));
        }
    } else {
        all_samples.reserve(rename_mapping.size());
        for (auto &sample_block : sample_pos_mut_local) {
            for (const auto &sample : sample_block) {
                auto iter = rename_mapping.find(sample.name);
                if (iter != rename_mapping.end()) {
                    all_samples.push_back(std::move(sample));
                    all_samples.back().name = iter->second;
                    rename_mapping.erase(iter);
                }
            }
            sample_block.clear();
        }
    }
    for (const auto &name : rename_mapping) {
        fprintf(stderr, "sample %s not found \n", name.first.c_str());
    }
    size_t max_name_len=0;
    for(const auto& samp:all_samples) {
        max_name_len=std::max(max_name_len,samp.name.size());
    }
    auto inbuf_size=max_name_len+ref.size()+8;
    auto out_fd=open(output_fa_file.c_str(),O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
    std::mutex f_mutex;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,all_samples.size(),all_samples.size()/num_threads), Batch_Printer{all_samples,ref,inbuf_size,out_fd,f_mutex});
    close(out_fd);
}
