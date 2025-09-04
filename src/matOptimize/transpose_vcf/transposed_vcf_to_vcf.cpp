#include "../mutation_annotated_tree.hpp"
#include "transpose_vcf.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_pipeline.h>
#include <tbb/global_control.h>
#include <utility>
#include <unordered_map>
#include <vector>
#include <tbb/info.h>
#include <zlib.h>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options.hpp>
namespace MAT = Mutation_Annotated_Tree;
std::string chrom;
std::string ref;
struct Pos_Mut {
    int position;
    uint8_t mut;
    Pos_Mut(int position,uint8_t mut):position(position),mut(mut) {}
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

struct Sample_Pos_Mut_Wrap {
    Sample_Pos_Mut& content;
    void add_Not_N(int position, uint8_t allele) {
        content.not_Ns.emplace_back(position, allele);
    }
    void add_N(int first, int second) {
        content.Ns.emplace_back(first, second);
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
std::string write_vcf_header(std::vector<Sample_Pos_Mut> &samples) {
    std::string out("##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (const auto &samp : samples) {
        out.push_back('\t');
        out += samp.name;
    }
    out += '\n';
    return out;
}
struct mut_iterator {
#ifndef NDEBUG
    std::vector<Pos_Mut>::const_iterator not_N_begin;
#endif
    std::vector<Pos_Mut>::const_iterator not_N_iter;
    std::vector<Pos_Mut>::const_iterator not_N_end;
#ifndef NDEBUG
    std::vector<std::pair<int, int>>::const_iterator N_begin;
#endif
    std::vector<std::pair<int, int>>::const_iterator N_iter;
    std::vector<std::pair<int, int>>::const_iterator N_end;
    mut_iterator() {}
    mut_iterator(const Sample_Pos_Mut &mut_ele, int position) {
        not_N_end = mut_ele.not_Ns.end();
        N_end = mut_ele.Ns.end();
#ifndef NDEBUG
        not_N_begin = mut_ele.not_Ns.begin();
#endif
        not_N_iter = std::lower_bound(mut_ele.not_Ns.begin(), not_N_end,
                                      Pos_Mut{position, 0});
        assert((mut_ele.not_Ns.begin() == not_N_iter) ||
               (not_N_iter - 1)->position < position);
        assert((mut_ele.not_Ns.end() == not_N_iter) ||
               not_N_iter->position >= position);
        N_iter = std::lower_bound(mut_ele.Ns.begin(), N_end,
                                  std::make_pair(position, position),
                                  [](const std::pair<int, int> &first,
        const std::pair<int, int> &second) {
            return first.first < second.first;
        });
#ifndef NDEBUG
        N_begin = mut_ele.Ns.begin();
#endif
        assert((mut_ele.Ns.begin() == N_iter) ||
               (N_iter - 1)->first < position);
        assert((mut_ele.Ns.end() == N_iter) || N_iter->first >= position);
        if (N_iter != N_end && N_iter->second < position) {
            N_iter++;
        }
        if (N_iter != mut_ele.Ns.begin() && (N_iter - 1)->second >= position) {
            N_iter--;
        }
    }
    uint8_t get_allele(int position) {
        uint8_t ret_val = 0;
        if (not_N_iter != not_N_end) {
            if (not_N_iter->position == position) {
                ret_val = not_N_iter->mut;
            }
            if (not_N_iter->position <= position) {
                not_N_iter++;
            }
            assert(not_N_iter == not_N_end||not_N_iter->position > position);
        }
        assert(not_N_iter == not_N_begin ||
               (not_N_iter - 1)->position <= position);
        if (N_iter != N_end && N_iter->second < position) {
            N_iter++;
        }
        assert(N_iter<=N_end);
        assert(N_iter == N_begin || (N_iter - 1)->second < position);
        assert(N_iter == N_end || N_iter->first >= position ||
               N_iter->second >= position);
        if (N_iter != N_end && N_iter->second >= position&&position>=N_iter->first) {
            ret_val = 0xf;
        }
        assert(!(ret_val&0xf0));
        return ret_val;
    }
};
struct vcf_rows {
    std::vector<uint8_t> alleles;
    uint allele_count[16];
    vcf_rows(size_t count) {
        alleles.reserve(count);
        std::memset(allele_count, 0, 16 * 4);
    }
    void add_allele(uint8_t allele) {
        alleles.push_back(allele);
        allele_count[allele]++;
    }
};
std::vector<std::pair<int8_t, uint>> count_alleles(const uint *allele_count) {
    std::vector<std::pair<int8_t, uint>> out;
    for (uint8_t allele = 1; allele < 16; allele++) {
        if (allele_count[allele]) {
            out.emplace_back(allele, allele_count[allele]);
        }
    }
    std::sort(out.begin(), out.end(),
              [](const std::pair<int8_t, uint> &first,
    const std::pair<int8_t, uint> &second) {
        return first.second > second.second;
    });
    return out;
}
void make_id(int8_t ref, uint pos, std::vector<std::pair<int8_t, uint>> &alts,
             std::string &id) {
    assert(!alts.empty());
    id.push_back(ref);
    id += std::to_string(pos) + MAT::get_nuc(alts[0].first);
    for (size_t idx=1; idx<alts.size(); idx++) {
        id += ",";
        id.push_back(ref);
        id += std::to_string(pos) + MAT::get_nuc(alts[idx].first);
    }
}

void make_info(std::vector<std::pair<int8_t, uint>> &alts, uint leaf_count,
               std::string &alt_count_str) {
    assert(!alts.empty());
    alt_count_str+="AC=";
    alt_count_str += std::to_string(alts[0].second);
    for (size_t idx=1; idx<alts.size(); idx++) {
        alt_count_str += ",";
        alt_count_str += std::to_string(alts[idx].second);
    }
    alt_count_str += ";AN=" + std::to_string(leaf_count);
}
void make_alt_str(std::vector<std::pair<int8_t, uint>> &alts,
                  std::string &alt_str) {
    assert(!alts.empty());
    alt_str += MAT::get_nuc(alts[0].first);

    for (size_t idx=1; idx<alts.size(); idx++) {
        alt_str += ",";
        alt_str += MAT::get_nuc(alts[idx].first);
    }

}

void make_allele_codes(std::vector<std::pair<int8_t, uint>> &alts,
                       int *allele_code) {
    memset(allele_code, 0, 16 * 4);
    for (size_t idx = 0; idx < alts.size(); idx++) {
        allele_code[alts[idx].first] = idx + 1;
    }
}

void write_individual_row(std::string *out, int position, const vcf_rows &row) {
    int8_t ref_nuc = ref[position-1];
    out->append(chrom);
    out->push_back('\t');
    out->append(std::to_string(position));
    out->push_back('\t');
    std::vector<std::pair<int8_t, uint>> alts =
                                          count_alleles(row.allele_count);
    make_id(ref_nuc, position, alts, *out);
    out->push_back('\t');
    out->push_back(ref_nuc);
    out->push_back('\t');
    make_alt_str(alts, *out);
    out->append("\t.\t.\t");
    auto sample_size = row.alleles.size();
    make_info(alts, sample_size, *out);
    int allele_codes[16];
    make_allele_codes(alts, allele_codes);
    // fprintf(vcf_file, "\tGT");
    out->append("\tGT");
    for (uint i = 0; i < sample_size; i++) {
        int8_t allele = row.alleles[i];
        out->append("\t");
        out->append(std::to_string(allele_codes[allele]));
    }
    out->append("\n");
}
std::pair<uint8_t *, size_t> compress(std::string *in) {
    auto malloc_size = compressBound(in->size()) + 32;
    uint8_t *buffer = (uint8_t *)malloc(malloc_size);
    z_stream stream;
    stream.opaque = nullptr;
    stream.zalloc = nullptr;
    stream.zfree = nullptr;
    stream.next_in = (uint8_t *)const_cast<char *>(in->data());
    stream.avail_in = in->size();
    stream.avail_out=malloc_size;
    stream.next_out = buffer;
    //auto err =
    deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 | 16,
                 9, Z_DEFAULT_STRATEGY);
    //assert(err == Z_OK);
    //err =
    deflate(&stream, Z_FINISH);
    //assert(err == Z_STREAM_END);
    size_t writen_size = stream.total_out;
    buffer = (uint8_t *)realloc((void *)buffer, writen_size);
    deflateEnd(&stream);
    return std::make_pair(buffer, writen_size);
}

struct File_Writer {
    FILE *fh;
    void operator()(std::pair<uint8_t *, size_t> in) const {
        fwrite((void *)in.first, 1, in.second, fh);
        free(in.first);
    }
};
typedef std::pair<std::vector<int>::const_iterator,
        std::vector<int>::const_iterator>
        iter_range;
struct Output_Genotypes {
    const std::vector<Sample_Pos_Mut> &all_samples;
    std::pair<uint8_t *, size_t> operator()(const iter_range &range) const {
        auto nPos = range.second - range.first;
        auto nSamples = all_samples.size();
        std::vector<vcf_rows> rows(nPos, vcf_rows(nSamples));
        std::vector<mut_iterator> iters;
        iters.reserve(all_samples.size());
        for (const auto &samp : all_samples) {
            iters.emplace_back(samp, *range.first);
        }
        size_t samp_idx_start = 0;
        while (samp_idx_start < nSamples) {
            auto samp_idx_end = std::min(64 + samp_idx_start, nSamples);
            for (auto pos_iter = range.first; pos_iter < range.second;
                    pos_iter++) {
                for (size_t samp_idx = samp_idx_start; samp_idx < samp_idx_end;
                        samp_idx++) {
                    rows[pos_iter-range.first].add_allele(
                        iters[samp_idx].get_allele(*pos_iter));
                }
            }
            samp_idx_start = samp_idx_end;
        }
        std::string in;
        in.reserve((range.second - range.first) *
                   (100 + all_samples.size() * 2));
        auto iter = range.first;
        for (auto &row : rows) {
            write_individual_row(&in, *iter, row);
            iter++;
        }
        return compress(&in);
    }
};

struct Pos_Iter_Gen {
    std::vector<int>::const_iterator &last_pos_iter;
    std::vector<int>::const_iterator end_iter;
    iter_range operator()(tbb::flow_control &fc) const {
        if (last_pos_iter == end_iter) {
            fc.stop();
        }
        auto first = last_pos_iter;
        auto second = std::min(end_iter, first + 64);
        last_pos_iter = second;
        return std::make_pair(first, second);
    }
};

void write_header(std::vector<Sample_Pos_Mut>& all_samples,FILE* fh) {
    auto hdr_str=write_vcf_header(all_samples);
    auto compressed_hdr=compress(&hdr_str);
    fwrite(compressed_hdr.first,1,compressed_hdr.second,fh);
    free(compressed_hdr.first);
}
namespace po = boost::program_options;
int main(int argc, char **argv) {
    po::options_description desc{"Options"};
    uint32_t num_cores = tbb::info::default_concurrency();
    std::string output_vcf_path;
    std::string input_path;
    std::string rename_file;
    uint32_t num_threads;
    std::string reference;
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
    ("vcf,v", po::value<std::string>(&output_vcf_path)->required(), "Output VCF file (in uncompressed or gzip-compressed .gz format) ")
    ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
    ("reference,r", po::value<std::string>(&reference)->required(),"Reference file")
    ("samples,s", po::value<std::string>(&rename_file),"rename file")
    ("output_path,i", po::value<std::string>(&input_path)->required(), "Load transposed VCF");
    po::options_description all_options;
    all_options.add(desc);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }
    std::unordered_map<std::string,std::string> rename_mapping;
    if(rename_file!="") {
        parse_rename_file(rename_file,rename_mapping);
    }
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads);
    load_reference(reference.c_str(), chrom, ref);
    All_Sample_Appender appender;
    load_mutations(input_path.c_str(), 80, appender);
    std::vector<Sample_Pos_Mut> all_samples;
    if(rename_mapping.empty()) {
        for (auto &sample_block : sample_pos_mut_local) {
            all_samples.insert(all_samples.end(),
                               std::make_move_iterator(sample_block.begin()),
                               std::make_move_iterator(sample_block.end()));
        }
    } else {
        all_samples.reserve(rename_mapping.size());
        for (auto &sample_block : sample_pos_mut_local) {
            for(const auto& sample:sample_block) {
                auto iter=rename_mapping.find(sample.name);
                if(iter!=rename_mapping.end()) {
                    all_samples.push_back(std::move(sample));
                    all_samples.back().name=iter->second;
                    rename_mapping.erase(iter);
                }
            }
            sample_block.clear();
        }
    }
    for(const auto& name:rename_mapping) {
        fprintf(stderr,"sample %s not found \n",name.first.c_str());
    }

    /*for (const auto& sample : all_samples) {
        fprintf(stdout, "%s:",sample.name.c_str());
        for (const auto& mut : sample.not_Ns) {
            fprintf(stdout,"\t%d%c",mut.position,MAT::get_nuc(mut.mut));
        }
        for (const auto& N : sample.Ns) {
            fprintf(stdout, "\t%d-%dN",N.first,N.second);
        }
        putc('\n',stdout);
    }*/
    std::vector<bool> pos_with_mut(ref.size(),false);
    for (const auto& samp : all_samples) {
        for (const auto& mut : samp.not_Ns) {
            pos_with_mut[mut.position]=1;
        }
        for (const auto& N_range : samp.Ns) {
            for(int pos=N_range.first; pos<=N_range.second; pos++) {
                pos_with_mut[pos]=1;
            }
        }
    }
    std::vector<int> pos_mut_idx;
    for (size_t idx = 0; idx < ref.size(); idx++) {
        if (pos_with_mut[idx]) {
            pos_mut_idx.push_back(idx);
        }
    }
    std::vector<int>::const_iterator pos_mut_iter = pos_mut_idx.begin();
    std::vector<int>::const_iterator pos_mut_end = pos_mut_idx.end();
    FILE *vcf_out = fopen(output_vcf_path.c_str(), "w");
    write_header(all_samples, vcf_out);
    tbb::parallel_pipeline(
        1000, tbb::make_filter<void, iter_range>(
            tbb::filter_mode::serial_in_order,
            Pos_Iter_Gen{pos_mut_iter, pos_mut_end}) &
        tbb::make_filter<iter_range, std::pair<uint8_t *, size_t>>(
            tbb::filter_mode::parallel, Output_Genotypes{all_samples}) &
        tbb::make_filter<std::pair<uint8_t *, size_t>, void>(
            tbb::filter_mode::serial_in_order, File_Writer{vcf_out}));
    fclose(vcf_out);
}
