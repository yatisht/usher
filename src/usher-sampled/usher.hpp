#include <atomic>
#include <cstddef>
#define USHER
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/Fitch_Sankoff.hpp"
#include <cassert>
#include <cstdint>
#include "src/matOptimize/check_samples.hpp"
#include <vector>
#define LEVEL_T uint8_t
#define MAX_LEVLEL UINT8_MAX
#define IDX_TREE_IDX_T uint32_t
#define EMPTY_POS UINT32_MAX
#define SIZE_MULT 8
#define POSITION_TAG 30
#define FS_RESULT_TAG 31
#pragma once
namespace MAT = Mutation_Annotated_Tree;
//#define MPI_TRACE
#ifdef MPI_TRACE
#define mpi_trace_print(...) fprintf(stderr, __VA_ARGS__ )
#else
#define mpi_trace_print(...)
#endif
void check_order(MAT::Mutations_Collection& in);
struct To_Place_Sample_Mutation {
    int position;
    uint8_t chrom_idx;
    uint8_t mut_nuc;
    union {
        struct {
            uint8_t par_nuc;
            uint8_t descendent_possible_nuc;
        };
        uint16_t range;
    };
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc), range(0) {
        assert(mut_nuc == 0xf);
    }
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc,
                             uint8_t par_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc),par_nuc(par_nuc)
        ,descendent_possible_nuc(0xf) {
        assert(mut_nuc != 0xf);
    }
    To_Place_Sample_Mutation(int position, uint8_t chrom_idx, uint8_t mut_nuc,
                             uint8_t par_nuc,uint8_t descendant_nuc)
        : position(position), chrom_idx(chrom_idx), mut_nuc(mut_nuc),par_nuc(par_nuc)
        ,descendent_possible_nuc(descendant_nuc) {
        assert(mut_nuc != 0xf);
    }
    To_Place_Sample_Mutation()=default;
    int get_end_range() const {
        if (mut_nuc==0xf) {
            return position+range;
        } else {
            return position;
        }
    }
};
struct Sample_Muts {
    size_t sample_idx;
    std::vector<To_Place_Sample_Mutation> muts;
    int sorting_key1;
    int sorting_key2;
};
struct Clade_info {
    std::vector<std::string> best_clade_assignment;
    std::vector<std::vector<std::string>> clade_assignments;
    bool valid;
    Clade_info():valid(false) {}
};
void Sample_Input(const char *name, std::vector<Sample_Muts> &sample_mutations,
                  MAT::Tree &tree,std::vector<mutated_t>&,bool,
                  std::vector<std::string>& fields
                  ,const std::unordered_set<std::string>& samples_in_condensed_nodes);
#ifndef NDEBUG
Mutation_Set get_mutations(const MAT::Node *main_tree_node);
void check_descendant_nuc(const MAT::Node* node);
#endif
struct output_options {
    bool print_uncondensed_tree;
    std::string outdir;
    bool retain_original_branch_len;
    bool only_one_tree;
    size_t print_subtrees_single;
    size_t print_subtrees_size;
    std::string dout_filename;
    bool detailed_clades;
    bool redo_FS_Min_Back_Mutations;
};
bool final_output(MAT::Tree& T,const output_options& options,int t_idx,std::vector<Clade_info>& assigned_clades,
                  size_t sample_start_idx,size_t sample_end_idx,std::vector<std::string>& low_confidence_samples,std::vector<mutated_t>& position_wise_out);
void place_sample_leader(std::vector<Sample_Muts> &sample_to_place,
                         MAT::Tree &main_tree, int batch_size,
                         std::atomic_size_t &curr_idx,
                         int parsimony_increase_threshold,bool dry_run,
                         FILE *placement_stats_file,
                         int max_parsimony,size_t max_uncertainty,
                         std::vector<std::string>& low_confidence_samples,
                         std::vector<Clade_info>& samples_clade,
                         size_t sample_start_idx,std::vector<size_t>* idx_map,
                         bool do_print=false
                        ) ;
void fix_parent(Mutation_Annotated_Tree::Tree &tree);
void convert_mut_type(const std::vector<MAT::Mutation> &in,
                      std::vector<To_Place_Sample_Mutation> &out);
void assign_descendant_muts(MAT::Tree &in);
void assign_levels(MAT::Node* root);
void follower_place_sample(MAT::Tree &main_tree,int batch_size,bool dry_run);
void check_parent(MAT::Node* root,MAT::Tree& tree);
void find_moved_node_neighbors(int radius,size_t start_idx, MAT::Tree& tree, size_t cur_idx,std::vector<size_t>& node_to_search_idx);
int follower_recieve_positions( std::vector<mutated_t>& to_recieve);
void get_pos_samples_old_tree(MAT::Tree& tree,std::vector<mutated_t>& output);
void MPI_reassign_states(MAT::Tree& tree,const std::vector<mutated_t>& mutations,int start_position,bool initial=false);
struct Leader_Thread_Options {
    std::string protobuf_in;
    std::string tree_in;
    bool override_mutations;
    bool print_parsimony_scores;
    bool sort_before_placement_1;
    bool sort_before_placement_2;
    int max_uncertainty;
    int max_parsimony;
    bool sort_by_ambiguous_bases;
    bool reverse_sort;
    int parsimony_threshold;
    size_t first_n_samples;
    size_t keep_n_tree;
    bool collapse_tree;
    output_options out_options;
    float desired_optimization_msec;
    int initial_optimization_radius;
    int last_optimization_minutes;
    std::string vcf_filename;
    bool no_add;
    std::string diff_file_name;
    std::string reference_file_name;
};
int set_descendant_count(MAT::Node* root);
void discretize_mutations(const std::vector<To_Place_Sample_Mutation> &in,
                          const MAT::Mutations_Collection &shared_mutations,
                          MAT::Node *parent_node,
                          MAT::Mutations_Collection &out);
bool sort_samples(const Leader_Thread_Options& options,std::vector<Sample_Muts>& samples_to_place, MAT::Tree& tree,size_t sample_start_idx);
void place_sample_multiple_tree(
    std::vector<Sample_Muts> &sample_to_place,
    std::vector<MAT::Tree>& trees,
    FILE *placement_stats_file, int max_trees);
void distribute_positions(std::vector<mutated_t>& output);
void reassign_state_local(MAT::Tree& tree,const std::vector<mutated_t>& mutations,bool initial=false);
void remove_absent_leaves(MAT::Tree& tree,std::unordered_set<std::string>& present);
void print_annotation(const MAT::Tree &T, const output_options &options,
                      const std::vector<Clade_info> &assigned_clades,
                      size_t sample_start_idx, size_t sample_end_idx,
                      size_t num_annotations);
void min_back_reassign_state_local(MAT::Tree& tree,const std::vector<mutated_t>& mutations);
void MPI_min_back_reassign_states(MAT::Tree &tree,const std::vector<mutated_t> &mutations,int start_position);
int count_back_mutation(const MAT::Tree& tree);
void place_sample_sequential(
    std::vector<Sample_Muts> &sample_to_place, MAT::Tree &main_tree,
    bool dry_run, FILE *placement_stats_file, int max_parsimony,
    size_t max_uncertainty, std::vector<std::string> &low_confidence_samples,
    std::vector<Clade_info> &samples_clade, size_t sample_start_idx,
    bool do_print, FILE *printer_out);
void clean_up_leaf(std::vector<MAT::Node*>& dfs);
void prep_tree(MAT::Tree &tree);
void load_diff_for_usher(
    const char *input_path,std::vector<Sample_Muts>& all_samples,
    std::vector<mutated_t>& position_wise_out, MAT::Tree &tree,
    const std::string& fasta_fname,std::vector<std::string> & samples);