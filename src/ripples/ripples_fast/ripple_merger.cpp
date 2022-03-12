#include "ripples.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cassert>
#include <stdio.h>
#include <emmintrin.h>
#include <vector>
typedef unsigned short __v8hu __attribute__((__vector_size__(16)));
typedef  short __v8h __attribute__((__vector_size__(16)));
typedef  int __v4i __attribute__((__vector_size__(16)));
typedef char __v16b __attribute__((__vector_size__(16)));
typedef unsigned short __v8hu_u __attribute__((__vector_size__(16), __aligned__(1)));
static int acceptor (const Mut_Count_Out_t &counts, size_t i, size_t j,
                     size_t curr_node_idx, size_t node_size,
                     size_t num_mutations) {
    auto first_half_mut =
        counts[i * node_size + curr_node_idx].count_before_exclusive();
    auto second_half_mut = counts[num_mutations * node_size + curr_node_idx]
                           .count_before_exclusive() -
                           counts[(j - 1) * node_size + curr_node_idx]
                           .count_before_inclusive();
    return first_half_mut + second_half_mut;
}
static int donor(const Mut_Count_Out_t &counts, size_t i, size_t j,
                 size_t curr_node_idx, size_t node_size,
                 size_t num_mutations) {
    return counts[(j - 1) * node_size + curr_node_idx]
           .count_before_inclusive() -
           (counts[i * node_size + curr_node_idx].count_before_exclusive());
}

static void push_val(std::vector<int> &filtered_idx,
                     std::vector<unsigned short>& filtered_par_score,
                     int start_idx, int mask, __v8hu par_score) {
    for (int indi_idx = 0; indi_idx < 8; indi_idx++) {
        if (mask & (1 << (2 * indi_idx))) {
            filtered_idx.push_back(start_idx + indi_idx);
            filtered_par_score.push_back(par_score[indi_idx]);
        }
    }
}
static unsigned short min_8(__v8h in) {
    __v8h temp1=(__v8h)__builtin_ia32_pshufd((__v4i)in,0x4e);
    auto min1=__builtin_ia32_pminsw128(temp1,in);
    temp1=__builtin_ia32_pshuflw(min1,0x4e);
    min1=__builtin_ia32_pminsw128(temp1,min1);
    temp1=__builtin_ia32_pshuflw(min1,0xb1);
    min1=__builtin_ia32_pminsw128(temp1,min1);
    auto out=min1[0];
#ifndef NDEBUG
    unsigned short test=in[0];
    for (int i=0; i<8; i++) {
        test=std::min((unsigned short)in[i],test);
    }
    assert(test==out);
#endif
    return out;
}
static std::pair<int,int> filter(const Ripples_Mapper_Output_Interface &out_ifc, size_t i,
                                 size_t j, size_t node_size, size_t num_mutations,
                                 int idx_start,int idx_end,
                                 std::vector<int> &donor_filtered_idx,
                                 std::vector<unsigned short> &donor_filtered_par_score,
                                 std::vector<int> &acceptor_filtered_idx,
                                 std::vector<unsigned short> &acceptor_filtered_par_score,
                                 int pasimony_threshold) {
    /*int donor_min_par = pasimony_threshold + 1;
    int acceptor_min_par = pasimony_threshold + 1;
    const auto& counts = out_ifc.mut_count_out;
    for (; idx_start < idx_end;
         idx_start++) {
        auto donor_par =
            donor(counts, i, j, idx_start, node_size, num_mutations);
        auto acceptor_par =
            acceptor(counts, i, j, idx_start, node_size, num_mutations);

        if (donor_par <= pasimony_threshold) {
            donor_min_par = std::min(donor_min_par, donor_par);
            donor_filtered_idx.push_back(idx_start);
            donor_filtered_par_score.push_back(donor_par);
        }
        if (acceptor_par <= pasimony_threshold) {
            acceptor_min_par = std::min(acceptor_min_par, acceptor_par);
            acceptor_filtered_idx.push_back(idx_start);
            acceptor_filtered_par_score.push_back(acceptor_par);
        }
    }*/
#ifndef NDEBUG
    int donor_par_min_debug= pasimony_threshold + 1;
    int acceptor_par_min_debug= pasimony_threshold + 1;
    std::vector<int> donor_idx_debug(donor_filtered_idx);
    std::vector<int> acceptor_idx_debug(acceptor_filtered_idx);
#endif
    unsigned short threshold_par = pasimony_threshold + 1;
    __v8hu load_cmp_mask= {1,2,3,4,5,6,7,8};
    __v8hu threshold_par_vec=threshold_par-(__v8hu) {};
    __v8h donor_min_par_vec=(__v8h)threshold_par_vec;
    __v8h acceptor_min_par_vec=(__v8h)threshold_par_vec;

    const auto& counts = out_ifc.mut_count_out;
    const auto* base_i_offset=counts.data()+i*node_size;
    const auto* base_j_offset=counts.data()+(j - 1) * node_size;
    const auto* base_all_mut_offset=counts.data()+num_mutations*node_size;
    __v8hu exlusive_count_extract_mask=0x7fff-(__v8hu) {};
    //__m128i_u inclusive_extract_mask=_mm_set1_epi16(0x8000);
    for (; idx_start < (idx_end-8);
            idx_start+=8) {
        __v8hu first_half_raw=*(__v8hu_u*)(base_i_offset+idx_start);
        __v8hu first_half_exclusive=exlusive_count_extract_mask&first_half_raw;
        __v8hu all_muts_raw=*((__v8hu_u*)(base_all_mut_offset+idx_start));
        __v8hu all_muts=exlusive_count_extract_mask& all_muts_raw;
        __v8hu second_half_before_raw=*((__v8hu_u*)(base_j_offset+idx_start));
        __v8hu second_half_before_included_flag=(second_half_before_raw>>15);
        __v8hu second_half_before_inclusive=(second_half_before_raw&exlusive_count_extract_mask)+ second_half_before_included_flag;
        __v8hu acceptor_par=all_muts-second_half_before_inclusive+first_half_exclusive;
        __v8hu donor_par=second_half_before_inclusive-first_half_exclusive;
        int donor_pass=__builtin_ia32_pmovmskb128((__v16b)(donor_par<threshold_par_vec));
        int acceptor_pass=__builtin_ia32_pmovmskb128((__v16b)(acceptor_par<threshold_par_vec));
#ifndef NDEBUG
        auto end_idx=std::min(8,idx_end-idx_start);
        for (int indi_idx=0; indi_idx<end_idx; indi_idx++) {
            auto donor_par_i =
                donor(counts, i, j, idx_start+indi_idx, node_size, num_mutations);
            auto acceptor_par_i =
                acceptor(counts, i, j, idx_start+indi_idx, node_size, num_mutations);
            if (donor_par_i<=pasimony_threshold) {
                donor_idx_debug.push_back(idx_start+indi_idx);
            }
            if (acceptor_par_i<=pasimony_threshold) {
                acceptor_idx_debug.push_back(idx_start+indi_idx);
            }
            assert(donor_par_i==donor_par[indi_idx]);
            assert(acceptor_par_i==acceptor_par[indi_idx]);
            assert((donor_par_i<=pasimony_threshold)==((donor_pass&(1<<(2*indi_idx)))!=0));
            assert((acceptor_par_i<=pasimony_threshold)==((acceptor_pass&(1<<(2*indi_idx)))!=0));
            donor_par_min_debug=std::min(donor_par_min_debug,(int)donor_par_i);
            acceptor_par_min_debug=std::min(acceptor_par_min_debug,(int)acceptor_par_i);
        }
        for (; end_idx<8; end_idx++) {
            assert(donor_par[end_idx]==0x7fff);
            assert(acceptor_par[end_idx]==0x7fff);
        }
#endif
        if (donor_pass) {
            push_val(donor_filtered_idx, donor_filtered_par_score, idx_start, donor_pass, donor_par);
            donor_min_par_vec=__builtin_ia32_pminsw128(donor_min_par_vec,(__v8h)donor_par);
        }
        if (acceptor_pass) {
            push_val(acceptor_filtered_idx, acceptor_filtered_par_score, idx_start, acceptor_pass, acceptor_par);
            acceptor_min_par_vec=__builtin_ia32_pminsw128(acceptor_min_par_vec,(__v8h)acceptor_par);
        }
        assert(donor_idx_debug.size()==donor_filtered_idx.size());
        assert(acceptor_idx_debug.size()==acceptor_filtered_idx.size());
    }
    __v8hu first_half_raw=*(__v8hu_u*)(base_i_offset+idx_start);
    __v8hu first_half_exclusive=exlusive_count_extract_mask&first_half_raw;
    __v8hu all_muts_raw=*((__v8hu_u*)(base_all_mut_offset+idx_start));
    __v8hu all_muts=exlusive_count_extract_mask& all_muts_raw;
    __v8hu second_half_before_raw=*((__v8hu_u*)(base_j_offset+idx_start));
    __v8hu second_half_before_included_flag=(second_half_before_raw>>15);
    __v8hu second_half_before_inclusive=(second_half_before_raw&exlusive_count_extract_mask)+ second_half_before_included_flag;
    __v8hu acceptor_par=all_muts-second_half_before_inclusive+first_half_exclusive;
    __v8hu donor_par=second_half_before_inclusive-first_half_exclusive;
    __v8hu load_mask=load_cmp_mask>((unsigned short)(idx_end-idx_start))-(__v8hu) {};
    load_mask>>=1;
    acceptor_par|=load_mask;
    acceptor_par&=exlusive_count_extract_mask;
    donor_par|=load_mask;
    donor_par&=exlusive_count_extract_mask;
    int donor_pass=__builtin_ia32_pmovmskb128((__v16b)(donor_par<threshold_par_vec));
    int acceptor_pass=__builtin_ia32_pmovmskb128((__v16b)(acceptor_par<threshold_par_vec));
    if (donor_pass) {
        push_val(donor_filtered_idx, donor_filtered_par_score, idx_start, donor_pass, donor_par);
        donor_min_par_vec=__builtin_ia32_pminsw128(donor_min_par_vec,(__v8h)donor_par);
    }
    if (acceptor_pass) {
        push_val(acceptor_filtered_idx, acceptor_filtered_par_score, idx_start, acceptor_pass, acceptor_par);
        acceptor_min_par_vec=__builtin_ia32_pminsw128(acceptor_min_par_vec,(__v8h)acceptor_par);
    }
    auto donor_min_par=min_8(donor_min_par_vec);
    auto acceptor_min_par=min_8(acceptor_min_par_vec);
    assert(donor_min_par==donor_par_min_debug);
    assert(acceptor_min_par==acceptor_par_min_debug);

    return std::make_pair(donor_min_par,acceptor_min_par);
}
static void threshold_parsimony(const Ripples_Mapper_Output_Interface &out_ifc,
                                size_t node_size, size_t num_mutations,
                                std::vector<int> &idx,
                                std::vector<unsigned short> &par_score,
                                std::vector<Recomb_Node> &filtered,
                                int pasimony_threshold,
                                const std::vector<MAT::Node*>& nodes_to_search) {
    const auto& counts = out_ifc.mut_count_out;
    for (size_t filtered_idx = 0; filtered_idx < idx.size(); filtered_idx++) {
        if (par_score[filtered_idx] <= pasimony_threshold) {
            auto idx_start = idx[filtered_idx];
            filtered.emplace_back(nodes_to_search[idx_start],
                                  counts[node_size * num_mutations + idx_start]
                                  .count_before_exclusive(),
                                  par_score[filtered_idx],
                                  out_ifc.is_sibling[idx_start]);
        }
    }
}

static void find_pairs(
    const std::vector<Recomb_Node> &donor_nodes,
    const std::vector<Recomb_Node> &acceptor_nodes,
    const std::vector<MAT::Mutation> &pruned_sample_mutations, int i, int j,
    int parsimony_threshold,
    const MAT::Tree &T, tbb::concurrent_vector<Recomb_Interval> &valid_pairs) {
    bool has_printed = false;

    for (auto d : donor_nodes) {
        /*if (T.is_ancestor(nid_to_consider, d.name)) {
            //raise(SIGTRAP);
            continue;
        }*/
        for (auto a : acceptor_nodes) {
            /*if (T.is_ancestor(nid_to_consider, a.name)) {
                //raise(SIGTRAP);
                continue;
            }*/
            // Ensure donor and acceptor are not the same and
            // neither of them is a descendant of the recombinant
            // node total parsimony is less than the maximum allowed
            if (parsimony_threshold < d.parsimony + a.parsimony) {
                break;
            }
            if ((d.node != a.node)
                    /*&& (d.name != nid_to_consider) &&
                        (a.name != nid_to_consider) &&
                        (orig_parsimony >= d.parsimony + a.parsimony +
                                               parsimony_improvement)
                                               */) {
                int start_range_high = pruned_sample_mutations[i].position;
                int start_range_low =
                    (i >= 1) ? pruned_sample_mutations[i - 1].position : 0;

                // int end_range_high = pruned_sample_mutations[j].position;
                int end_range_high = 1e9;
                int end_range_low =
                    (j >= 1) ? pruned_sample_mutations[j - 1].position : 0;
                Pruned_Sample donor;
                donor.sample_mutations.clear();
                Pruned_Sample acceptor;
                acceptor.sample_mutations.clear();

                auto donor_node = d.node;

                for (auto anc : T.rsearch(donor_node->identifier, true)) {
                    for (auto mut : anc->mutations) {
                        donor.add_mutation(mut);
                    }
                }

                for (auto mut : donor.sample_mutations) {
                    if ((mut.position > start_range_low) &&
                            (mut.position <= start_range_high)) {
                        bool in_pruned_sample = false;
                        for (auto mut2 : pruned_sample_mutations) {
                            if (mut.position == mut2.position) {
                                in_pruned_sample = true;
                            }
                        }
                        if (!in_pruned_sample) {
                            start_range_low = mut.position;
                        }
                    }
                    if ((mut.position > end_range_low) &&
                            (mut.position <= end_range_high)) {
                        bool in_pruned_sample = false;
                        for (auto mut2 : pruned_sample_mutations) {
                            if (mut.position == mut2.position) {
                                in_pruned_sample = true;
                            }
                        }
                        if (!in_pruned_sample) {
                            end_range_high = mut.position;
                        }
                    }
                }

                for (auto mut : pruned_sample_mutations) {
                    if ((mut.position > start_range_low) &&
                            (mut.position <= start_range_high)) {
                        bool in_pruned_sample = false;
                        for (auto mut2 : donor.sample_mutations) {
                            if (mut.position == mut2.position) {
                                in_pruned_sample = true;
                            }
                        }
                        if (!in_pruned_sample) {
                            start_range_low = mut.position;
                        }
                    }
                    if ((mut.position > end_range_low) &&
                            (mut.position <= end_range_high)) {
                        bool in_pruned_sample = false;
                        for (auto mut2 : donor.sample_mutations) {
                            if (mut.position == mut2.position) {
                                in_pruned_sample = true;
                            }
                        }
                        if (!in_pruned_sample) {
                            end_range_high = mut.position;
                        }
                    }
                }

                // tbb_lock.lock();
                valid_pairs.push_back(
                    Recomb_Interval(d, a, start_range_low, start_range_high,
                                    end_range_low, end_range_high));
                // tbb_lock.unlock();

                has_printed = true;
                break;
            }
        }
        if (has_printed) {
            break;
        }
    }
}
struct check_breakpoint {
    const Ripples_Mapper_Output_Interface &out_ifc;
    const std::vector<MAT::Mutation> &pruned_sample_mutations;
    int skip_start_idx;
    int skip_end_idx;
    size_t node_size;
    int pasimony_threshold;
    const std::vector<MAT::Node *> &nodes_to_search;
    const MAT::Tree &T;
    tbb::concurrent_vector<Recomb_Interval> &valid_pairs;
    void operator()(std::pair<int,int> in) const {
        int i=in.first;
        int j=in.second;
        if (i==2&&j==14) {
            //raise(SIGTRAP);
        }

        size_t num_mutations = pruned_sample_mutations.size();
        std::vector<int> donor_filtered_idx;
        std::vector<unsigned short> donor_filtered_par_score;
        std::vector<int> acceptor_filtered_idx;
        std::vector<unsigned short> acceptor_filtered_par_score;
        donor_filtered_idx.reserve(nodes_to_search.size());
        donor_filtered_par_score.reserve(nodes_to_search.size());
        acceptor_filtered_idx.reserve(nodes_to_search.size());
        acceptor_filtered_par_score.reserve(nodes_to_search.size());
        auto min_first = filter(
                             out_ifc, i, j, node_size, num_mutations, 0, skip_start_idx,
                             donor_filtered_idx, donor_filtered_par_score, acceptor_filtered_idx,
                             acceptor_filtered_par_score, pasimony_threshold);
        auto min_second = filter(
                              out_ifc, i, j, node_size, num_mutations, skip_end_idx,
                              nodes_to_search.size(), donor_filtered_idx,
                              donor_filtered_par_score, acceptor_filtered_idx,
                              acceptor_filtered_par_score, pasimony_threshold);
        auto donor_min = std::min(min_first.first, min_second.first);
        auto acceptor_min=std::min(min_first.second,min_second.second);
        if (acceptor_min+donor_min>pasimony_threshold) {
            return;
        }
        std::vector<Recomb_Node> donor_filtered;
        std::vector<Recomb_Node> acceptor_filtered;
        donor_filtered.reserve(donor_filtered_idx.size());
        acceptor_filtered.reserve(donor_filtered_idx.size());
        threshold_parsimony(out_ifc, node_size, num_mutations,
                            donor_filtered_idx, donor_filtered_par_score,
                            donor_filtered, pasimony_threshold - acceptor_min,
                            nodes_to_search);
        threshold_parsimony(out_ifc, node_size, num_mutations,
                            acceptor_filtered_idx, acceptor_filtered_par_score,
                            acceptor_filtered, pasimony_threshold - donor_min,
                            nodes_to_search);
        std::sort(acceptor_filtered.begin(), acceptor_filtered.end());
        std::sort(donor_filtered.begin(), donor_filtered.end());
        find_pairs(donor_filtered, acceptor_filtered, pruned_sample_mutations,
                   i, j, pasimony_threshold, T, valid_pairs);
    }
};

struct search_position {
    const std::vector<MAT::Mutation> &pruned_sample_mutations;
    int &i;
    int &j;
    int branch_len;
    int min_range;
    int max_range;
    int last_i;
    bool is_j_end_of_range(int start_range_high, int total_size) const {
        return j >= total_size || total_size - (j - i) < branch_len ||
               pruned_sample_mutations[j-1].position - start_range_high >
               max_range;
    }
    std::pair<int, int> operator()(tbb::flow_control &fc) const {
        int start_range_high=0;// = pruned_sample_mutations[i].position;
        int total_size = pruned_sample_mutations.size();
        // i end
        if (i > last_i) {
            fc.stop();
            return std::make_pair(0, 0);
        }
        j++;
        // j end
        if (i!=-1) {
            start_range_high = pruned_sample_mutations[i].position;
        }

        while (i==-1||is_j_end_of_range(start_range_high, total_size)) {
            i++;
            if (i > last_i) {
                fc.stop();
                return std::make_pair(0, 0);
            }
            start_range_high = pruned_sample_mutations[i].position;
            j = i + branch_len;
            while (j < total_size && pruned_sample_mutations[j - 1].position <
                    start_range_high + min_range) {
                j++;
            }
        }
        return std::make_pair(i, j);
    }
};

void ripplrs_merger(const Pruned_Sample &pruned_sample,
                    const std::vector<int> & idx_map,
                    const std::vector<MAT::Node *> &nodes_to_search,
                    size_t node_size, int pasimony_threshold,
                    const MAT::Tree &T,
                    tbb::concurrent_vector<Recomb_Interval> &valid_pairs,
                    const Ripples_Mapper_Output_Interface &out_ifc,
                    int nthreads, int branch_len, int min_range,
                    int max_range) {
    auto pruned_node = pruned_sample.sample_name;
    int skip_start_idx=std::abs(idx_map[pruned_node->dfs_idx]);
    //The next index after the one corresponding to dfs_idx -1
    int skip_end_idx=std::abs(idx_map[pruned_node->dfs_end_idx]);

    const auto &sample_mutations = pruned_sample.sample_mutations;
    int i = -1;
    int j = 0;
    const auto last_pos = sample_mutations.back().position - min_range;
    int last_i = sample_mutations.size() - branch_len;
    while (last_i > 0 && sample_mutations[last_i].position > last_pos) {
        last_i--;
    }

    tbb::parallel_pipeline(
        nthreads + 1,
        tbb::make_filter<void, std::pair<int, int>>(
            tbb::filter::serial_in_order,
            search_position{sample_mutations, i, j, branch_len, min_range,
                            max_range, last_i}) &
        tbb::make_filter<std::pair<int, int>, void>(
            tbb::filter::parallel,
            check_breakpoint{out_ifc, sample_mutations,
                             skip_start_idx,skip_end_idx,
                             node_size, pasimony_threshold,
                             nodes_to_search, T, valid_pairs})

    );


}