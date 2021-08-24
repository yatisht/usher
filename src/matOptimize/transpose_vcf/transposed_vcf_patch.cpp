#include <cstddef>
#include <functional>
#include <tbb/task_group.h>
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
//#include <malloc.h>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <utility>
#include <vector>
#include <sys/types.h>
#include <signal.h>

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
        auto node=tree.get_node(name);
        if (node) {
            bfs_idx = node->bfs_index;
        }else {
            bfs_idx=-1;
        }
        
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
                fprintf(stderr, "At %d, N incremented to %d - %d, nuc %d\n", position, N_iter==N_end?INT_MAX:N_iter->first, N_iter==N_end?INT_MAX:N_iter->second, ret_val);
            }
        }
        assert(N_iter <= N_end);
        assert(N_iter == N_begin || (N_iter - 1)->second < position);
        if(!(N_iter == N_end || N_iter->first >= position ||
                N_iter->second >= position)) {
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
    mutated_t alleles;
    row_t() {}
    row_t(int pos) : mut(pos) {}
};
typedef tbb::flow::function_node<std::vector<row_t> *> assigner_t;
#define CHUNK_SIZ 64
struct output_vcf_rows {
    size_t start_idx;
    size_t end_idx;
    const std::vector<int>& positions;
    assigner_t& out;
    const std::vector<Sample_Pos_Mut>& all_samples;
    void operator()()const {
        size_t idx=start_idx;
        std::vector<mut_iterator> iters;
        for (const auto &samp : all_samples) {
            if (samp.bfs_idx!=-1) {                
                iters.emplace_back(samp, positions[idx]);
            }
        }
        while (idx<end_idx) {
            std::vector<row_t> *rows = new std::vector<row_t>;
            rows->reserve(CHUNK_SIZ);
            for (int count = 0; count < CHUNK_SIZ; count++) {
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
};
struct Assigner {
    const std::vector<backward_pass_range> &child_idx_range;
    const std::vector<forward_pass_range> &parent_idx;
    FS_result_per_thread_t& fs_result;
    const std::vector<MAT::Node*>& bfs_ordered_nodes;
    void operator()(std::vector<row_t> *rows) const {
        auto& this_content=fs_result.local();
        this_content.init(child_idx_range.size());
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
                                     row.alleles, this_content);
        }
        delete rows;
    }
};
void get_sample_mut(const char *input_path,std::vector<Sample_Pos_Mut>& all_samples, MAT::Tree &tree) {
    sample_pos_mut_local_t sample_pos_mut_local;
    All_Sample_Appender appender{tree,sample_pos_mut_local};
    load_mutations(input_path, 80, appender);
    for (auto &sample_block : sample_pos_mut_local) {
        all_samples.insert(all_samples.end(),
                           std::make_move_iterator(sample_block.begin()),
                           std::make_move_iterator(sample_block.end()));
    }
}
void print_mut(const Sample_Pos_Mut& to_print) {
    for (const auto& mut : to_print.not_Ns) {
        fprintf(stderr, "%d:%d\t",mut.position,mut.mut);
    }
    for (const auto& mut : to_print.Ns) {
        fprintf(stderr, "%d-%d:N\t",mut.first,mut.second);
    }
    fprintf(stderr, "\n");
}
void assign_state(std::vector<Sample_Pos_Mut>& all_samples, MAT::Tree &tree,FS_result_per_thread_t&output,const std::vector<MAT::Node*> bfs_ordered_nodes) {
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
    tbb::task_group transposer_group;
    size_t last_end_idx=0;
    int iter_thread_count=num_threads/6;
    fprintf(stderr, "Using %d transposer threads",iter_thread_count);
    size_t chunk_size=pos_mut_idx.size()/(iter_thread_count+1);
    for (int thread_idx=0; thread_idx<iter_thread_count; thread_idx++) {
        auto this_end_idx=last_end_idx+chunk_size;
        transposer_group.run(output_vcf_rows{last_end_idx,this_end_idx,pos_mut_idx,assigner,all_samples});
        last_end_idx=this_end_idx;
    }
    transposer_group.run(output_vcf_rows{last_end_idx,pos_mut_idx.size(),pos_mut_idx,assigner,all_samples});
    fprintf(stderr, "last chunk size %zu, other chunk size %zu",pos_mut_idx.size()-last_end_idx,chunk_size);
    transposer_group.wait();
    g.wait_for_all();
    deallocate_FS_cache(output);
}
void asign_and_fill(std::vector<Sample_Pos_Mut>& all_samples, MAT::Tree &tree,std::vector<MAT::Node*> bfs_ordered_nodes) {
    FS_result_per_thread_t output;
    assign_state(all_samples, tree, output, bfs_ordered_nodes);
    //raise(SIGUSR1);
    print_memory();
    //malloc_stats();
    fill_muts(output, bfs_ordered_nodes);
}
void add_ambiguous_mutation(const char *input_path, MAT::Tree &tree) {
    tree.uncondense_leaves();
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    //p_idx=tree.get_node("LR882438.1|260113|20-08-26")->bfs_index;
    p_idx=-1;
    std::vector<Sample_Pos_Mut> all_samples;
    get_sample_mut(input_path, all_samples, tree);
    /*tbb::parallel_for(tbb::blocked_range<size_t>(0,pos_mut_idx.size(),100),Output_Genotypes{pos_mut_idx,all_samples,child_idx_range,parent_idx,output});*/
    asign_and_fill(all_samples, tree, bfs_ordered_nodes);
    auto par_score = tree.get_parsimony_score();
    fprintf(stderr, "Before condensing %zu\n", par_score);
    //raise(SIGUSR1);
    print_memory();
    //malloc_stats();
    recondense_tree(tree);
}
