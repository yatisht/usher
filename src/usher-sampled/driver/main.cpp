#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "src/usher-sampled/usher.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <chrono>
#include <climits>
#include <complex>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <random>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include "version.hpp"
#ifdef __linux
#include <sys/prctl.h>
#endif
bool use_bound;
int process_count;
int this_rank;
unsigned int num_threads;
std::atomic_bool interrupted(false);
void fix_condensed_nodes(MAT::Tree *tree);
namespace po = boost::program_options;
static void clean_tree_for_placement(MAT::Tree& tree){
    auto dfs = tree.depth_first_expansion();
    for (auto node : dfs) {
        if(node->is_leaf()&&!node->is_root()) {
            for (auto& mut : node->mutations) {
                mut.set_mut_one_hot(mut.get_all_major_allele());
            }
        } else {
            node->mutations.remove_invalid();
        }
    }
}
void leader_thread_optimization(MAT::Tree& tree,std::vector<mutated_t>& position_wise_out,
                                std::atomic_size_t& curr_idx,int& optimization_radius, size_t start_idx,FILE* ignored_file,float desired_optimization_msec,bool is_last) {
    /*auto nodes=tree.depth_first_expansion();
    clean_up_leaf(nodes);
    fprintf(stderr, "init parsimony %zu\n",tree.get_parsimony_score());
    std::vector<mutated_t> other(MAT::Mutation::refs.size());
    get_pos_samples_old_tree(tree, other);
    for (size_t idx=0; idx<other.size(); idx++) {
        bool have_err=false;
        std::unordered_map<long, nuc_one_hot> old_map(position_wise_out[idx].begin(),position_wise_out[idx].end());
        for (const auto& new_mut : other[idx]) {
            auto iter=old_map.find(new_mut.first);
            if (iter==old_map.end()) {
                have_err=true;
                fprintf(stderr, "pos %zu, sample %s was ref %c, but changed to %c\n",
                    idx,tree.get_node_name_for_log_output(new_mut.first).c_str(),
                    MAT::get_nuc(MAT::Mutation::refs[idx]),MAT::get_nuc(new_mut.second));
            }else {
                if(iter->second!=new_mut.second&&iter->second!=0xf){
                    have_err=true;
                    fprintf(stderr, "pos %zu, sample %s was %c, but changed to %c\n",
                    idx,tree.get_node_name_for_log_output(new_mut.first).c_str(),
                    MAT::get_nuc(iter->second),MAT::get_nuc(new_mut.second));
                }
                old_map.erase(iter);
            }
        }
        for (const auto & new_mut : old_map) {
            if(new_mut.second==0xf){
                continue;
            }
            have_err=true;
            fprintf(stderr, "pos %zu, sample %s was %c, but changed to ref %c\n",
                idx,tree.get_node_name_for_log_output(new_mut.first).c_str(),
                MAT::get_nuc(new_mut.second),MAT::get_nuc(MAT::Mutation::refs[idx]));
        }
        if(have_err){
            raise(SIGTRAP);
        }
    }*/
    size_t last_parsimony_score=SIZE_MAX;
    std::default_random_engine g;
    tree.max_level=tree.get_max_level();
    auto optimiation_start=std::chrono::steady_clock::now();
    auto optimization_end=optimiation_start+std::chrono::milliseconds((long)desired_optimization_msec);
    bool timeout=false;
    if (is_last) {
        optimization_radius=4;
    }
    bool is_first=true;
    do  {
        bool distributed = process_count > 1;
        fprintf(stderr, "Main sent optimization prep\n");
        if (process_count == 1) {
            if(is_first){
            reassign_state_local(tree,position_wise_out);
            tree.populate_ignored_range();
            }
        } else {
            MPI_Bcast(&optimization_radius, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (is_first) {
                MPI_reassign_states(tree, position_wise_out, 0);
                tree.populate_ignored_range();
            } else {
                tree.MPI_send_tree();
            }
        }
        is_first=false;
        fprintf(stderr, "Main parsimony score %zu",tree.get_parsimony_score());
        fprintf(stderr, "Main sent optimization prep done\n");
        std::vector<size_t> node_to_search_idx;
        find_moved_node_neighbors(optimization_radius,
                                  start_idx, tree,
                                  curr_idx.load(), node_to_search_idx);
        fprintf(stderr, "Main found nodes to move\n");
        std::shuffle(node_to_search_idx.begin(),node_to_search_idx.end(),g);

        while (!node_to_search_idx.empty()) {
            std::vector<size_t> deferred_nodes_out;
            adjust_all(tree);
            fprintf(stderr, "Main sent tree_optimizing\n");
            optimize_tree_main_thread(
                node_to_search_idx, tree, optimization_radius, ignored_file,
                false, 1, deferred_nodes_out, distributed, optimization_end, true,
                true, true);
            node_to_search_idx.clear();
            auto dfs = tree.depth_first_expansion();
            node_to_search_idx.reserve(deferred_nodes_out.size());
            for (auto idx : deferred_nodes_out) {
                auto node=tree.get_node(idx);
                if(node) {
                    node_to_search_idx.push_back(node->dfs_index);
                }
            }
            distributed = false;
            auto new_parsimony_score=tree.get_parsimony_score();
            fprintf(stderr, "Last parsimony score %lu\n",new_parsimony_score);
            if(new_parsimony_score>last_parsimony_score
                    || std::chrono::steady_clock::now()>optimization_end) {
                timeout=true;
                break;
            }
            if (new_parsimony_score==last_parsimony_score) {
                if(optimization_radius>tree.max_level){
                    timeout=true;
                }
                break;
            }
            last_parsimony_score=new_parsimony_score;
        }
        if (is_last) {
            optimization_radius=2*optimization_radius;
        }
    } while(is_last&&!timeout);
    if (is_last) {
        int to_send=0;
        MPI_Bcast(&to_send, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    auto finish_time=std::chrono::steady_clock::now();
    int next_optimization_radius;
    if (finish_time>optimization_end) {
        next_optimization_radius=optimization_radius-1;
    } else {
        float actual_msec=std::chrono::duration_cast<std::chrono::milliseconds>(finish_time-optimiation_start).count();
        float time_ratio=sqrt(desired_optimization_msec/actual_msec);
        next_optimization_radius=(optimization_radius*time_ratio);
        fprintf(stderr,"Next radius %d, ratio %f",next_optimization_radius,time_ratio);
    }
    optimization_radius=std::max(next_optimization_radius,2);
    clean_tree_for_placement(tree);
}

static int leader_thread(
    int batch_size_per_process,
    Leader_Thread_Options& options
) {
    int optimization_radius=options.initial_optimization_radius;
    boost::filesystem::path path(options.out_options.outdir);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(options.out_options.outdir);
    options.out_options.outdir = path.generic_string();
    fprintf(stderr, "Output path %s\n",options.out_options.outdir.c_str());
    options.override_mutations|=options.tree_in!="";
    auto start_time=std::chrono::steady_clock::now();
    MAT::Tree tree;
    if (options.tree_in!="") {
        tree=Mutation_Annotated_Tree::create_tree_from_newick(options.tree_in);
    } else {
        if(!MAT::load_mutation_annotated_tree(options.protobuf_in,tree)) {
            exit(EXIT_FAILURE);
        }
    }
    tree.uncondense_leaves();
    std::vector<Sample_Muts> samples_to_place;
    std::vector<mutated_t> position_wise_out;
    std::vector<std::string> samples;
    const std::unordered_set<std::string> samples_in_condensed_nodes;
    if(options.diff_file_name!=""&&options.reference_file_name!=""){
        load_diff_for_usher(options.diff_file_name.c_str(), samples_to_place, position_wise_out,tree,options.reference_file_name,samples);
    }else {
        if(options.vcf_filename==""){
            fprintf(stderr, "Expect either VCF file or MAPLE file\n");
            exit(EXIT_FAILURE);
        }
        Sample_Input(options.vcf_filename.c_str(),samples_to_place,tree,position_wise_out,options.override_mutations,samples,samples_in_condensed_nodes,options.duplicate_prefix);
    }
    samples_to_place.resize(std::min(samples_to_place.size(),options.first_n_samples));
    if(samples_to_place.empty()){
        fprintf(stderr,"No samples to place\n");
        exit(EXIT_FAILURE);
    }
    size_t sample_start_idx=samples_to_place[0].sample_idx;
    size_t sample_end_idx=samples_to_place.back().sample_idx+1;
    fprintf(stderr, "Sample start idx %zu, end index %zu\n",sample_start_idx,sample_end_idx);
    if (options.tree_in!="") {
        std::unordered_set<std::string> sample_set(samples.begin(),samples.end());
        remove_absent_leaves(tree, sample_set);
        if(tree.root->children.size()){
        if (process_count==1) {
            reassign_state_local(tree,position_wise_out,true);
        } else {
            distribute_positions(position_wise_out);
            MPI_reassign_states(tree, position_wise_out, 0,true);
        }
        }
        clean_tree_for_placement(tree);
    }
    tree.condense_leaves();
    fix_parent(tree);
    tree.check_leaves();
    bool have_ambiguous_ref=false;
    for (int position=1; position<MAT::Mutation::refs.size(); position++) {
        auto nuc=MAT::Mutation::refs[position];
        if (nuc&(nuc-1)) {
            fprintf(stderr, "\nWARNING: Ref nuc @ %d : %c is ambiguous\n", position,MAT::get_nt(MAT::Mutation::refs[position]));
            have_ambiguous_ref=true;
        }
    }
    if (have_ambiguous_ref) {
        fprintf(stderr, "WARNING: Reference contain ambiguous nucleotide, optimization is disabled\n");
        options.initial_optimization_radius=0;
        optimization_radius=0;
    }
    if (options.initial_optimization_radius>0) {
        for (auto& pos : position_wise_out) {
            pos.erase(std::remove_if(pos.begin(), pos.end(), [sample_start_idx](const std::pair<long, nuc_one_hot>& in) {
                return in.first<(long)sample_start_idx;
            }),pos.end());
        }
        get_pos_samples_old_tree(tree, position_wise_out);

    }
    for(size_t idx=0;idx<position_wise_out.size();idx++){
        bool informative=false;
        for (const auto &samp : position_wise_out[idx]) {
            if (samp.second!=0xf) {
                informative=true;
                break;
            }
        }
        if (!informative) {
            position_wise_out[idx].clear();
        }
    }
    std::vector<std::string> low_confidence_samples;
    std::vector<Clade_info> samples_clade(samples_to_place.size());
    for (auto& temp : samples_clade) {
        temp.valid=false;
    }
    if (samples_to_place.empty()) {
        fprintf(stderr, "Found no new samples\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fprintf(stderr, "Found %zu missing samples.\n\n", samples_to_place.size());
    std::string placement_stats_filename = options.out_options.outdir + "/placement_stats.tsv";
    FILE *placement_stats_file = fopen(placement_stats_filename.c_str(), "w");

    if (options.no_add) {
        std::atomic_size_t curr_idx(0);
        assign_descendant_muts(tree);
        if (process_count>1) {
            tree.MPI_send_tree();
        }
        tree.breadth_first_expansion();
        place_sample_leader(samples_to_place, tree, 100, curr_idx, INT_MAX,
                            true, placement_stats_file, INT_MAX, INT_MAX,
                            low_confidence_samples, samples_clade,
                            sample_start_idx, nullptr, true);
        print_annotation(tree, options.out_options, samples_clade,
                         sample_start_idx, sample_end_idx,
                         tree.get_num_annotations());
        MPI_Finalize();
        return 0;
    }
    auto reordered=sort_samples(options, samples_to_place, tree,sample_start_idx);
    fprintf(stderr, "sorting done\n");
    std::vector<size_t> idx_map;
    std::vector<size_t>* idx_map_ptr=nullptr;
    if (reordered) {
        idx_map.resize(samples_to_place.size());
        for(size_t idx=0; idx<samples_to_place.size(); idx++) {
            idx_map[samples_to_place[idx].sample_idx-sample_start_idx]=idx;
        }
        idx_map_ptr=&idx_map;
    }
    std::atomic_size_t curr_idx(0);
    FILE* ignored_file=fopen("/dev/null", "w");
    use_bound=true;
    //samples_to_place.resize(1000);
    options.out_options.only_one_tree=options.keep_n_tree==1;
    if(options.keep_n_tree>1) {
        std::vector<MAT::Tree> trees{tree};
        place_sample_multiple_tree(samples_to_place, trees, placement_stats_file, options.keep_n_tree);
        for (size_t t_idx=0; t_idx<trees.size(); t_idx++) {
            std::vector<Clade_info> assigned_clades;
            std::vector<std::string> low_confidence_samples;
            final_output(trees[t_idx], options.out_options, t_idx, assigned_clades, sample_start_idx, sample_end_idx,low_confidence_samples,position_wise_out);
        }
        return 0;
    }
    while (true) {
        clean_tree_for_placement(tree);
        auto tree_size=prep_tree(tree);
        switch_to_serial_threshold=std::max((int)(tree_size*batch_size_per_process/(2*num_threads)),10);
        fprintf(stderr, "switch to serial search when there are less than %d descendants\n", switch_to_serial_threshold);
        if (process_count>1) {
            fprintf(stderr, "Main sending tree\n");
            tree.MPI_send_tree();
        }
        place_sample_leader(samples_to_place, tree, batch_size_per_process,
                            curr_idx, options.parsimony_threshold, false,
                            placement_stats_file,
                            options.max_parsimony, options.max_uncertainty,
                            low_confidence_samples, samples_clade,sample_start_idx,idx_map_ptr);
        bool is_last=false;
        tree.check_leaves();
        if (curr_idx >= samples_to_place.size()) {
            curr_idx.store(samples_to_place.size());
            is_last=true;
        }
        if (options.initial_optimization_radius > 0) {
            leader_thread_optimization(tree, position_wise_out, curr_idx, optimization_radius,
                                       sample_start_idx, ignored_file,is_last?60000*options.last_optimization_minutes:options.desired_optimization_msec,is_last);
            tree.check_leaves();
        }
        if (curr_idx<samples_to_place.size()) {
            int to_board_cast=0;
            if (process_count>1) {
                MPI_Bcast(&to_board_cast, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
        } else {
            int to_board_cast=1;
            if (process_count>1) {
                MPI_Bcast(&to_board_cast, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            break;
        }
    }
    fclose(placement_stats_file);
    fprintf(stderr, "Main finised place\n");
    auto dfs=tree.depth_first_expansion();
    clean_up_leaf(dfs);
    final_output(
        tree, options.out_options, 0, samples_clade, sample_start_idx, sample_end_idx, low_confidence_samples,position_wise_out);
    auto duration=std::chrono::steady_clock::now()-start_time;
    fprintf(stderr, "Took %ld msec\n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
    return 0;
}
void wait_debug();
int main(int argc, char **argv) {

    //signal(SIGSEGV, handler);
    int ignored;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    int batch_size_per_process;
    Leader_Thread_Options options;
    po::options_description desc{"Options"};
    int optimiation_minutes;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_cores) +
                                      " detected on this machine]";
    bool ignored_options;
    //std::vector<int> gdb_pids;
    desc.add_options()
    ("vcf,v", po::value<std::string>(&options.vcf_filename),"Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
    ("tree,t", po::value<std::string>(&options.tree_in)->default_value(""), "Input tree file")
    ("outdir,d", po::value<std::string>(&options.out_options.outdir)->default_value("."), "Output directory to dump output and log files [DEFAULT uses current directory]")
    ("load-mutation-annotated-tree,i",po::value<std::string>(&options.protobuf_in)->default_value(""),"Load mutation-annotated tree object")
    ("save-mutation-annotated-tree,o",po::value<std::string>(&options.out_options.dout_filename)->default_value(""),"Save output mutation-annotated tree object to the specified filename")
    ("sort-before-placement-1,s", po::bool_switch(&options.sort_before_placement_1)->default_value(false), \
     "Sort new samples based on computed parsimony score and then number of optimal placements before the actual placement [EXPERIMENTAL].")
    ("sort-before-placement-2,S", po::bool_switch(&options.sort_before_placement_2)->default_value(false), \
     "Sort new samples based on the number of optimal placements and then the parsimony score before the actual placement [EXPERIMENTAL].")
    ("sort-before-placement-3,A", po::bool_switch(&options.sort_by_ambiguous_bases)->default_value(false), \
     "Sort new samples based on the number of ambiguous bases [EXPERIMENTAL].")
    ("reverse-sort,r", po::bool_switch(&options.reverse_sort)->default_value(false), \
     "Reverse the sorting order of sorting options (sort-before-placement-1 or sort-before-placement-2) [EXPERIMENTAL]")
    ("collapse-tree,c", po::bool_switch(&options.collapse_tree)->default_value(false), \
     "Collapse internal nodes of the input tree with no mutations and condense identical sequences in polytomies into a single node and the save the tree to file condensed-tree.nh in outdir")
    ("collapse-output-tree,C", po::bool_switch(&ignored_options), \
     "Collapse internal nodes of the output tree with no mutations before the saving the tree to file final-tree.nh in outdir")
    ("max-uncertainty-per-sample,e", po::value(&options.max_uncertainty)->default_value(1e6), \
     "Maximum number of equally parsimonious placements allowed per sample beyond which the sample is ignored")
    ("max-parsimony-per-sample,E", po::value(&options.max_parsimony)->default_value(1e6), \
     "Maximum parsimony score of the most parsimonious placement(s) allowed per sample beyond which the sample is ignored")
    ("write-uncondensed-final-tree,u", po::bool_switch(&options.out_options.print_uncondensed_tree)->default_value(false), "Write the final tree in uncondensed format and save to file uncondensed-final-tree.nh in outdir")
    ("write-subtrees-size,k", po::value<size_t>(&options.out_options.print_subtrees_size)->default_value(0), \
     "Write minimum set of subtrees covering the newly added samples of size equal to this value")
    ("write-single-subtree,K", po::value<size_t>(&options.out_options.print_subtrees_single)->default_value(0), \
     "Similar to write-subtrees-size but produces a single subtree with all newly added samples along with random samples up to the value specified by this argument")
    ("multiple-placements,M", po::value(&options.keep_n_tree)->default_value(1), \
     "Create a new tree up to this limit for each possibility of parsimony-optimal placement")
    ("retain-input-branch-lengths,l", po::bool_switch(&options.out_options.retain_original_branch_len)->default_value(false), \
     "Retain the branch lengths from the input tree in out newick files instead of using number of mutations for the branch lengths.")
    ("no-add,n", po::bool_switch(&options.no_add), \
     "Do not add new samples to the tree")
    ("detailed-clades,D", po::bool_switch(&options.out_options.detailed_clades), \
     "In clades.txt, write a histogram of annotated clades and counts across all equally parsimonious placements")
    ("diff",po::value<std::string>(&options.diff_file_name),"diff file from MAPLE, to be used with reference sequence")
    ("ref",po::value<std::string>(&options.reference_file_name),"reference sequence, only needed for MAPLE")
    ("threads,T",po::value<uint32_t>(&num_threads)->default_value(num_cores),num_threads_message.c_str())
    ("reduce-back-mutation,B",po::bool_switch(&options.out_options.redo_FS_Min_Back_Mutations)->default_value(false),
     "Reassign states of internal nodes to reduce back mutation count.")
    ("version", "Print version number")
    ("help,h", "Print help messages")
    ("optimization_radius", po::value(&options.initial_optimization_radius)->default_value(4),
     "The search radius for optimization when parsimony score increase exceeds the threshold"
     "Set to 0 to disable optimiation"
     "Only newly placed samples and nodes within this radius will be searched")
    ("optimization_minutes", po::value(&optimiation_minutes)->default_value(5),
     "Optimization time of each iterations in minutes"
     "Set to 0 to disable optimiation"
     "Only newly placed samples and nodes within this radius will be searched")
    ("last_optimization_minutes", po::value(&options.last_optimization_minutes)->default_value(120),
     "Optimization radius for the last round")
    ("batch_size_per_process",po::value(&batch_size_per_process)->default_value(5),
     "The number of samples each process search simultaneously")
    ("parsimony_threshold",po::value(&options.parsimony_threshold)->default_value(100000),
     "Optimize after the parsimony score increase by this amount")
    ("first_n_samples",po::value(&options.first_n_samples)->default_value(SIZE_MAX),"[TESTING ONLY] Only place first n samples")
    ("no-ignore-prefix",po::value<std::string>(&options.duplicate_prefix),"prefix samples already in the tree to force placement")
    //("gdb_pid,g",po::value(&gdb_pids)->multitoken(),"gdb pids for attaching")
    ;
    po::variables_map vm;
    //wait_debug();
    bool have_error=false;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        // Return with error code 1 unless the user specifies help
        have_error=true;
    }
    if (have_error||(options.vcf_filename==""&&(options.diff_file_name==""||options.reference_file_name==""))) {
        if (this_rank==0) {
            if (vm.count("version")) {
                std::cout << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
            } else {
                std::cout << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
                std::cerr << desc << std::endl;
            }
        }
        if (this_rank==0) {
            if(vm.count("help")) {
                return 0;
            } else
                return 1;
        }
    }
    if (options.initial_optimization_radius<=0) {
        options.parsimony_threshold=INT_MAX;
    }
    num_threads=std::max(2u,num_threads);
    options.desired_optimization_msec=optimiation_minutes*60000;
    fprintf(stderr, "Num threads %d\n",num_threads);
#ifdef __linux
    prctl(0x59616d61,-1);
    fprintf(stderr, "rand %d of pid %d ",this_rank,getpid());
#endif
    tbb::task_scheduler_init init(num_threads);
    if (this_rank==0) {
        if (options.keep_n_tree>1&&process_count>1) {
            fprintf(stderr, "Multi-host parallelization of multiple placement is not supported\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        return leader_thread(batch_size_per_process,options);
    } else {
        MAT::Tree tree;
        std::vector<mutated_t> position_wise_mutations;
        int start_idx=0;
        if (!((options.tree_in==""&&options.no_add)||options.initial_optimization_radius==0)) {
            start_idx=follower_recieve_positions(position_wise_mutations);
        }
        if (options.tree_in!="") {
            MPI_reassign_states(tree, position_wise_mutations, start_idx,true);
            position_wise_mutations.clear();
            fprintf(stderr, "follwer recieving position again\n");
            start_idx=follower_recieve_positions(position_wise_mutations);
            fprintf(stderr, "follwer finished recieving position again\n");
        }
        if (options.sort_before_placement_1||options.sort_before_placement_2||options.no_add) {
            tree.MPI_receive_tree();
            if (!options.no_add) {
                assign_levels(tree.root);
                set_descendant_count(tree.root);
            }
            follower_place_sample(tree,100,true);
        }
        if (options.no_add) {
            MPI_Finalize();
            return 0;
        }
        while (true) {
            fprintf(stderr, "follwer recieving tree\n");
            tree.MPI_receive_tree();
            assign_levels(tree.root);
            set_descendant_count(tree.root);
            auto dfs=tree.depth_first_expansion();
            for (auto node : dfs) {
                //check_order(node->mutations);
#ifdef NDEBUG
                node->children.reserve(SIZE_MULT*node->children.size());
#endif
            }
            init.terminate();
            init.initialize(num_threads-1);
            follower_place_sample(tree,batch_size_per_process,false);
            init.terminate();
            init.initialize(num_threads);
            bool done=false;
            if (options.initial_optimization_radius>0) {
                bool is_last=false;
                do {
                    int optimization_radius;
                    MPI_Bcast(&optimization_radius, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    if (optimization_radius==0) {
                        done=true;
                        break;
                    }
                    fprintf(stderr, "Optimizing with radius %d\n",optimization_radius);
                    fprintf(stderr, "follower recieving trees\n");
                    if (is_last) {
                        tree.delete_nodes();
                        tree.MPI_receive_tree();
                    } else {
                        MPI_reassign_states(tree, position_wise_mutations, start_idx);
                    }
                    is_last=optimization_radius<0;
                    if (is_last) {
                        optimization_radius=-optimization_radius;
                    }
                    //tree.delete_nodes();
                    //tree.MPI_receive_tree();
                    adjust_all(tree);
                    use_bound=true;
                    fprintf(stderr, "follower start optimizing\n");
                    optimize_tree_worker_thread(tree, optimization_radius, false, true);
                } while (is_last);
            }
            int stop;
            MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (stop||done) {
                fprintf(stderr, "follower stoping %s\n", stop?"recieved stop":"last done");
                if (options.out_options.redo_FS_Min_Back_Mutations) {
                    MPI_min_back_reassign_states(tree, position_wise_mutations, start_idx);
                }
                break;
            }
            tree.delete_nodes();
        }
        MPI_Finalize();
    }
}
