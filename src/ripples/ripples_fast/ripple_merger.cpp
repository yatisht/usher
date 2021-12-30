#include "ripples.hpp"
#include <stdio.h>
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

static std::pair<int,int> filter(const Ripples_Mapper_Output_Interface &out_ifc, size_t i,
                  size_t j, size_t node_size, size_t num_mutations,
                  int idx_start,int idx_end,
                  std::vector<Recomb_Node> &donor_filtered,
                  std::vector<Recomb_Node> &acceptor_filtered,
                  int pasimony_threshold,
                  const std::vector<MAT::Node*>& nodes_to_search) {
    int donor_min_par = pasimony_threshold + 1;
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
            donor_filtered.emplace_back(nodes_to_search[idx_start],
                                  counts[node_size * num_mutations + idx_start]
                                      .count_before_exclusive(),
                                  donor_par, out_ifc.is_sibling[idx_start]);
        }
        if (acceptor_par <= pasimony_threshold) {
            acceptor_min_par = std::min(acceptor_min_par, acceptor_par);
            acceptor_filtered.emplace_back(nodes_to_search[idx_start],
                                  counts[node_size * num_mutations + idx_start]
                                      .count_before_exclusive(),
                                  acceptor_par, out_ifc.is_sibling[idx_start]);
        }
    }
    return std::make_pair(donor_min_par,acceptor_min_par);
}
static void threshold_parsimony(std::vector<Recomb_Node> &donor_filtered,
                                int pasimony_threshold) {
    auto end_iter = donor_filtered.end();
    auto start_iter = donor_filtered.begin();
    while (start_iter != end_iter) {
        if (start_iter->parsimony > pasimony_threshold) {
            *start_iter = *(end_iter - 1);
            end_iter--;
        } else {
            start_iter++;
        }
    }
    donor_filtered.erase(end_iter, donor_filtered.end());
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
        int i=in.first; int j=in.second;
        if (i==2&&j==14)
        {
            //raise(SIGTRAP);
        }
        
        size_t num_mutations = pruned_sample_mutations.size();
        std::vector<Recomb_Node> donor_filtered;
        std::vector<Recomb_Node> acceptor_filtered;
        auto min_first=filter(out_ifc,i,j,node_size,num_mutations,0,skip_start_idx,donor_filtered,acceptor_filtered,pasimony_threshold,nodes_to_search);
        auto min_second=filter(out_ifc,i,j,node_size,num_mutations,skip_end_idx,nodes_to_search.size(),donor_filtered,acceptor_filtered,pasimony_threshold,nodes_to_search);
        auto donor_min=std::min(min_first.first,min_second.first);
        auto acceptor_min=std::min(min_first.second,min_second.second);
        if (acceptor_min+donor_min>pasimony_threshold) {
            return;
        }
        threshold_parsimony(donor_filtered, pasimony_threshold - acceptor_min);
        threshold_parsimony(acceptor_filtered, pasimony_threshold - donor_min);
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
        int start_range_high;// = pruned_sample_mutations[i].position;
        int total_size = pruned_sample_mutations.size();
        // i end
        if (i > last_i) {
            fc.stop();
            return std::make_pair(0, 0);
        }
        j++;
        // j end
        if (i!=-1)
        {
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