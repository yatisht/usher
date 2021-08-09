#include <cstddef>
#include <functional>
#include <thread>
#define LOAD
#include <cstdio>
#include <deque>
#include <tbb/concurrent_vector.h>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for.h>
#include <tbb/pipeline.h>
#include "../Fitch_Sankoff.hpp"
#include "../mutation_annotated_tree.hpp"
#include "../tree_rearrangement_internal.hpp"
#include "transpose_vcf.hpp"
#include <algorithm>
#include <iterator>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <utility>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
struct Pos_Mut {
    int position;
    uint8_t mut;
    Pos_Mut(int position, uint8_t mut) : position(position), mut(mut) {}
    bool operator<(const Pos_Mut &other) const {
        return position < other.position;
    }
};

struct Sample_Pos_Mut {
    size_t bfs_idx;
    std::vector<Pos_Mut> not_Ns;
    std::vector<std::pair<int, int>> Ns;
    Sample_Pos_Mut(std::string &&name, MAT::Tree &tree) {
        bfs_idx = tree.get_node(name)->bfs_index;
    }
};

struct Sample_Pos_Mut_Wrap {
    Sample_Pos_Mut &content;
    void add_Not_N(int position, uint8_t allele) {
        assert(!(allele&0xf0));
        content.not_Ns.emplace_back(position, allele);
    }
    void add_N(int first, int second) {
        content.Ns.emplace_back(first, second);
    }
};

typedef tbb::enumerable_thread_specific<std::vector<Sample_Pos_Mut>>
    sample_pos_mut_local_t;
struct All_Sample_Appender {
    MAT::Tree &tree;
    sample_pos_mut_local_t& sample_pos_mut_local;
    Sample_Pos_Mut_Wrap set_name(std::string &&name) {
        sample_pos_mut_local_t::reference my_sample_pos_mut_local =
            sample_pos_mut_local.local();
        my_sample_pos_mut_local.emplace_back(std::move(name), tree);
        return Sample_Pos_Mut_Wrap{my_sample_pos_mut_local.back()};
    }
};
size_t p_idx;
struct mut_iterator {
#ifndef NDEBUG
    std::vector<Pos_Mut>::const_iterator not_N_begin;
#endif
    size_t idx;
    std::vector<Pos_Mut>::const_iterator not_N_iter;
    std::vector<Pos_Mut>::const_iterator not_N_end;
#ifndef NDEBUG
    std::vector<std::pair<int, int>>::const_iterator N_begin;
#endif
    std::vector<std::pair<int, int>>::const_iterator N_iter;
    std::vector<std::pair<int, int>>::const_iterator N_end;
    mut_iterator() {}
    mut_iterator(const Sample_Pos_Mut &mut_ele, int position) {
        idx = mut_ele.bfs_idx;
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
            while (not_N_iter != not_N_end&&not_N_iter->position < position) {
                not_N_iter++;
                if (idx==p_idx) {
                    fprintf(stderr, "At %d, not_N incremented to %d, nuc %d\n",position,not_N_iter==not_N_end?INT_MAX:not_N_iter->position,ret_val);
                }
            }
            if (not_N_iter != not_N_end&&not_N_iter->position == position) {
                ret_val = not_N_iter->mut;
                if (ret_val&0xf0) {
                    fprintf(stderr, "At %d, not_N match, nuc %d\n",position,ret_val);
                }
            }
            //assert(not_N_iter == not_N_end || not_N_iter->position > position);
        assert(not_N_iter == not_N_begin ||
               (not_N_iter - 1)->position <= position);
        while (N_iter != N_end && N_iter->second < position) {
            N_iter++;
            if (idx==p_idx) {
                fprintf(stderr, "At %d, N incremented to %d-%d, nuc %d\n",position,N_iter==N_end?INT_MAX:N_iter->first,position,N_iter==N_end?INT_MAX:N_iter->second,ret_val);
            }
        }
        assert(N_iter <= N_end);
        assert(N_iter == N_begin || (N_iter - 1)->second < position);
        if(!(N_iter == N_end || N_iter->first >= position ||
               N_iter->second >= position)){
                   fprintf(stderr, "%d;%d;%d",position,N_iter->first,N_iter->second);
                   assert(false);
               }
        if (N_iter != N_end && N_iter->second >= position &&
            position >= N_iter->first) {
            ret_val = 0xf;
            if (idx==p_idx) {
                fprintf(stderr, "At %d, N match\n",position);
            }
        }
        assert(!(ret_val & 0xf0));
        return ret_val;
    }
};

typedef std::pair<std::vector<int>::const_iterator,
                  std::vector<int>::const_iterator>
    iter_range;
struct row_t {
    MAT::Mutation mut;
    std::vector<std::pair<long, nuc_one_hot>> alleles;
    row_t() {}
    row_t(int pos) : mut(pos) {}
};
typedef tbb::flow::function_node<std::vector<row_t> *> assigner_t;
void output_vcf_rows(size_t start_idx,size_t end_idx,const std::vector<int>& positions, assigner_t& out,const std::vector<Sample_Pos_Mut>& all_samples){
    size_t idx=start_idx;
    std::vector<mut_iterator> iters;
    for (const auto &samp : all_samples) {
        iters.emplace_back(samp, positions[idx]);
    }
    while (idx<end_idx) {
        std::vector<row_t> *rows = new std::vector<row_t>;
        rows->reserve(128);
        for (int count = 0; count < 128; count++) {
            if (idx == end_idx) {
                break;
            }
            rows->emplace_back(positions[idx]);
            idx++;
        }
        for (size_t samp_idx = 0; samp_idx < iters.size(); samp_idx++) {
            for (auto &row : *rows) {
                nuc_one_hot allele = iters[samp_idx].get_allele(row.mut.get_position());
                if (allele) {
                    row.alleles.emplace_back(iters[samp_idx].idx, allele);
                }
            }
        }
        out.try_put(rows);
    }
}
/*struct Transposer {
    std::vector<mut_iterator> &iters;
    const std::vector<int> &positions;
    size_t &idx;
    std::vector<row_t> *operator()(tbb::flow_control &fc) const {
        if (idx == positions.size()) {
            fc.stop();
            return nullptr;
        }
        std::vector<row_t> *rows = new std::vector<row_t>;
        rows->reserve(128);
        for (int count = 0; count < 128; count++) {
            if (idx == positions.size()) {
                break;
            }
            rows->emplace_back(positions[idx]);
            idx++;
        }
        for (size_t samp_idx = 0; samp_idx < iters.size(); samp_idx++) {
            for (auto &row : *rows) {
                nuc_one_hot allele = iters[samp_idx].get_allele(row.mut.get_position());
                if (allele) {
                    row.alleles.emplace_back(iters[samp_idx].idx, allele);
                }
            }
        }
        return rows;
    }
};*/
struct Assigner {
    const std::vector<backward_pass_range> &child_idx_range;
    const std::vector<forward_pass_range> &parent_idx;
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output;
    const std::vector<MAT::Node*>& bfs_ordered_nodes;
    void operator()(std::vector<row_t> *rows) const {
        for (auto &row : *rows) {
            /*if (row.mut.position==3037) {
                for (const auto& samp : row.alleles) {
                    fprintf(stdout, "%s:%d\n",bfs_ordered_nodes[samp.first]->identifier.c_str(),samp.second);
                }
            }*/
            if (row.alleles.empty()) {
                continue;
            }
            std::sort(row.alleles.begin(), row.alleles.end(),
                      mutated_t_comparator());
            row.alleles.emplace_back(0, 0xf);
            Fitch_Sankoff_Whole_Tree(child_idx_range, parent_idx, row.mut,
                                     row.alleles, output);
        }
        delete rows;
    }
};
/*
struct Output_Genotypes {
    const std::vector<int>& positions;
    const std::vector<Sample_Pos_Mut> &all_samples;
    const std::vector<backward_pass_range>& child_idx_range;
    const std::vector<forward_pass_range>& parent_idx;
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>&
            output;
    void operator()(const tbb::blocked_range<size_t> &range) const {
        fprintf(stderr, "%zu ",range.size());
        auto nSamples = all_samples.size();
        std::vector<row_t> rows;
        rows.clear();
        rows.reserve(range.size());
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
                rows.emplace_back(positions[idx]);
        }
        std::vector<mut_iterator> iters;
        iters.clear();
        iters.reserve(all_samples.size());
        for (const auto &samp : all_samples) {
            iters.emplace_back(samp, rows[0].mut.get_position());
        }
        size_t samp_idx_start = 0;
        while (samp_idx_start < nSamples) {
            auto samp_idx_end = std::min(64 + samp_idx_start, nSamples);
            for (auto &row:rows) {
                for (size_t samp_idx = samp_idx_start; samp_idx < samp_idx_end;
                        samp_idx++) {
                    nuc_one_hot
allele=iters[samp_idx].get_allele(row.mut.get_position()); if (allele) {
                        row.alleles.emplace_back(iters[samp_idx].idx,allele);
                    }
                }
            }
            samp_idx_start = samp_idx_end;
        }
        for (auto &row:rows) {
            assert (!row.alleles.empty());
            std::sort(row.alleles.begin(),row.alleles.end(),mutated_t_comparator());
            row.alleles.emplace_back(0,0xf);
            Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx, row.mut,
row.alleles, output);
        }
    }
};
*/
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
void get_sample_mut(const char *input_path,std::vector<Sample_Pos_Mut>& all_samples, MAT::Tree &tree){
    sample_pos_mut_local_t sample_pos_mut_local;
    All_Sample_Appender appender{tree,sample_pos_mut_local};
    load_mutations(input_path, 80, appender);
    for (auto &sample_block : sample_pos_mut_local) {
        all_samples.insert(all_samples.end(),
                           std::make_move_iterator(sample_block.begin()),
                           std::make_move_iterator(sample_block.end()));
    }
}
void print_mut(const Sample_Pos_Mut& to_print){
    for (const auto& mut : to_print.not_Ns) {
        fprintf(stderr, "%d:%d\t",mut.position,mut.mut);
    }
    for (const auto& mut : to_print.Ns) {
        fprintf(stderr, "%d-%d:N\t",mut.first,mut.second);
    }
    fprintf(stderr, "\n");
}
void assign_state(std::vector<Sample_Pos_Mut>& all_samples, MAT::Tree &tree,std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>& output,const std::vector<MAT::Node*> bfs_ordered_nodes){
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    Fitch_Sankoff_prep(bfs_ordered_nodes, child_idx_range, parent_idx);
    //size_t idx = 0;
    std::vector<int> pos_mut_idx;
    pos_mut_idx.reserve(MAT::Mutation::refs.size());
    /*fprintf(stderr, "241:%d\n",MAT::Mutation::refs[241]);
    fprintf(stderr, "3037:%d\n",MAT::Mutation::refs[3037]);
    auto p_idx=tree.get_node("LR882438.1|260113|20-08-26")->bfs_index;
    for (const auto& sample : all_samples) {
        if (sample.bfs_idx==p_idx) {
            print_mut(sample);
        }
    }*/
    for (size_t idx = 0; idx < MAT::Mutation::refs.size(); idx++) {
        if (MAT::Mutation::refs[idx]) {
            pos_mut_idx.push_back(idx);
        }
    }
    tbb::flow::graph g;
    assigner_t assigner(g,tbb::flow::unlimited,Assigner{child_idx_range,parent_idx,output,bfs_ordered_nodes});
    std::vector<std::thread> iterator_threads;
    size_t last_end_idx=0;
    int iter_thread_count=num_threads/6;
    fprintf(stderr, "Using %d transposer threads",iter_thread_count);
    size_t chunk_size=pos_mut_idx.size()/(iter_thread_count+1);
    for (int thread_idx=0; thread_idx<iter_thread_count; thread_idx++) {
        auto this_end_idx=last_end_idx+chunk_size;
        iterator_threads.emplace_back(output_vcf_rows,last_end_idx,this_end_idx,std::ref(pos_mut_idx),std::ref(assigner),std::ref(all_samples));
        last_end_idx=this_end_idx;
    }
    iterator_threads.emplace_back(output_vcf_rows,last_end_idx,pos_mut_idx.size(),std::ref(pos_mut_idx),std::ref(assigner),std::ref(all_samples));
    fprintf(stderr, "last chunk size %zu, other chunk size %zu",pos_mut_idx.size()-last_end_idx,chunk_size);
    for (auto& thread : iterator_threads) {
        thread.join();
    }
    g.wait_for_all();
    /*
    std::vector<mut_iterator> iters;
    iters.reserve(all_samples.size());
    for (const auto &samp : all_samples) {
        iters.emplace_back(samp, pos_mut_idx[0]);
    }
    
    tbb::parallel_pipeline(
        num_threads,
        tbb::make_filter<void, std::vector<row_t> *>(tbb::filter::serial,
                                        Transposer{iters, pos_mut_idx, idx}) &
            tbb::make_filter<std::vector<row_t> *, void>(
                tbb::filter::parallel,
                Assigner{child_idx_range, parent_idx, output,bfs_ordered_nodes}));
    */
}

void add_ambiguous_mutation(const char *input_path, MAT::Tree &tree) {
    tree.uncondense_leaves();
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    //p_idx=tree.get_node("LR882438.1|260113|20-08-26")->bfs_index;
    p_idx=-1;
    std::vector<Sample_Pos_Mut> all_samples;
    get_sample_mut(input_path, all_samples, tree);
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    assign_state(all_samples, tree, output, bfs_ordered_nodes);
    /*tbb::parallel_for(tbb::blocked_range<size_t>(0,pos_mut_idx.size(),100),Output_Genotypes{pos_mut_idx,all_samples,child_idx_range,parent_idx,output});*/
    tbb::affinity_partitioner ap;
    // sort and fill
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                auto &to_refill = output[i];
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size());
                to_refill.clear();
                to_refill.shrink_to_fit();
            }
        },
        ap);
    auto par_score = tree.get_parsimony_score();
    fprintf(stderr, "Before condensing %zu\n", par_score);
    recondense_tree(tree);
}
