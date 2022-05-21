#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/usher-sampled/usher.hpp"
#include <asm-generic/errno-base.h>
#include <atomic>
#include <bits/types/FILE.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <tbb/task_scheduler_init.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace po = boost::program_options;

struct tree_info {
    MAT::Tree tree;
    std::unordered_set<std::string> condensed_nodes;
    tree_info(const tree_info& )=delete;
    tree_info()=default;
    ~tree_info(){
        fprintf(stderr, "deleting nodes\n");
        tree.delete_nodes();
    }
};
bool use_bound = true;
int process_count = 1;
int this_rank = 0;
unsigned int num_threads;
std::atomic_bool interrupted(false);
bool prep_single_tree(std::string path, tree_info &out) {
    if (!MAT::load_mutation_annotated_tree(path, out.tree)) {
        return false;
    }
    fix_parent(out.tree);
    prep_tree(out.tree);
    for (const auto &temp : out.tree.condensed_nodes) {
        for (const auto &str : temp.second) {
            out.condensed_nodes.emplace(str);
        }
    }
    return true;
}

typedef std::shared_ptr<std::vector<tree_info>> TreeCollectionPtr;
void reload_trees(TreeCollectionPtr &to_replace, const std::vector<std::string>& paths){
    tbb::task_scheduler_init init(num_threads);
    auto next = new std::vector<tree_info>(paths.size());
    for (size_t idx=0; idx<paths.size(); idx++) {
        if(!prep_single_tree(paths[idx], (*next)[idx])){
            fprintf(stderr, "Not reloaded\n");
            delete next;
            return;
        }
    }
    to_replace.reset(next);
}
void refresh_tree(TreeCollectionPtr &to_replace, std::fstream &tree_paths) {
    std::string buf;
    std::vector<std::string> paths;
    while (std::getline(tree_paths, buf)) {
        if (buf == "") {
            break;
        }
        paths.push_back(buf);        
    }
    reload_trees(to_replace,paths);
}
static void mgr_thread(TreeCollectionPtr &to_replace, std::string cmd_fifo_name,
                       std::string socket_name) {
    unlink(cmd_fifo_name.c_str());
    auto err=mkfifo(cmd_fifo_name.c_str(), S_IRWXU);
    if (err!=0) {
        perror("cannot create mgr fifo");
        exit(EXIT_FAILURE);
    }
    std::fstream mgr_f(cmd_fifo_name);
    while (true) {
        std::string buf;
        std::getline(mgr_f, buf);
        if (!mgr_f) {
            mgr_f=std::fstream(cmd_fifo_name);
        }
        if (buf == "stop") {
            interrupted.store(true);
            int sock_fd = socket(AF_UNIX, SOCK_STREAM, 0);
            struct sockaddr_un addr;
            addr.sun_family = AF_UNIX;
            strncpy(addr.sun_path, socket_name.c_str(), socket_name.size());
            auto res = connect(sock_fd, (struct sockaddr *)&addr, sizeof(addr));
            if (res != 0) {
                perror("Unable to turn off cleanly");
                exit(EXIT_FAILURE);
            }
            close(sock_fd);
            unlink(socket_name.c_str());
            return;
        } else if (buf == "reload") {
            refresh_tree(to_replace, mgr_f);
        } else {
            int new_thread_count;
            auto scaned = sscanf(buf.c_str(), "thread %d", &new_thread_count);
            if (scaned) {
                fprintf(stderr, "setting thread count to %d\n",
                        new_thread_count);
                num_threads = new_thread_count;
            }
        }
    }
}

static int create_socket(std::string socket_name) {
    auto sock_fd = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock_fd == -1) {
        perror("unable to create socket");
        exit(EXIT_FAILURE);
    }
    struct sockaddr_un addr;
    addr.sun_family = AF_UNIX;
    strncpy(addr.sun_path, socket_name.c_str(), socket_name.size()+1);
    auto bind_ret = bind(sock_fd, (struct sockaddr *)&addr, sizeof(addr));
    if (bind_ret != 0) {
        perror("unable to bind");
        exit(EXIT_FAILURE);
    }
    listen(sock_fd, 128);
    return sock_fd;
}

static void collect_done(std::unordered_map<int, int> &pid_to_fd_map,
                         bool blocking) {
    while (true) {
        int ignored;
        auto done_pid = waitpid(-1, &ignored, blocking ? 0 : WNOHANG);
        if (done_pid == -1||done_pid==0) {
            if (errno == ECHILD || (!blocking && errno == 0)) {
                return;
            } else {
                perror("error waiting");
                continue;
            }
        }
        auto iter = pid_to_fd_map.find(done_pid);
        if (iter == pid_to_fd_map.end()) {
            fprintf(stderr, "Cannot find corresponding fd of proc %d\n",
                    done_pid);
        } else {
            close(iter->second);
            pid_to_fd_map.erase(iter);
        }
    }
}
char cmd[] = "usher";
size_t get_options(FILE *f, Leader_Thread_Options &options) {
    std::vector<char *> args({cmd});
    while (true) {
        char *buf = NULL;
        size_t len = 0;
        auto char_read = getline(&buf, &len, f);
        if (char_read == -1) {
            return -1;
        } else if (char_read == 1) {
            break;
        } else {
            buf[char_read-1]=0;
            args.push_back(buf);
            fputs(buf,stderr);
            fputc('\n', stderr);
        }
    }
    po::options_description desc{"Options"};
    bool ignored_options;
    options.out_options.redo_FS_Min_Back_Mutations = false;
    options.initial_optimization_radius = 0;
    // std::vector<int> gdb_pids;
    size_t mat_idx = 0;
    options.print_parsimony_scores=false;
    desc.add_options()(
        "vcf,v", po::value<std::string>(&options.vcf_filename)->required(),
        "Input VCF file (in uncompressed or gzip-compressed .gz format) "
        "[REQUIRED]")(
        "outdir,d",
        po::value<std::string>(&options.out_options.outdir)->default_value("."),
        "Output directory to dump output and log files [DEFAULT uses current "
        "directory]")("mat-index,i", po::value<size_t>(&mat_idx)->required(),
                      "Load mutation-annotated tree object")(
        "save-mutation-annotated-tree,o",
        po::value<std::string>(&options.out_options.dout_filename)
            ->default_value(""),
        "Save output mutation-annotated tree object to the specified filename")(
        "sort-before-placement-1,s",
        po::bool_switch(&options.sort_before_placement_1)->default_value(false),
        "Sort new samples based on computed parsimony score and then number of "
        "optimal placements before the actual placement [EXPERIMENTAL].")(
        "sort-before-placement-2,S",
        po::bool_switch(&options.sort_before_placement_2)->default_value(false),
        "Sort new samples based on the number of optimal placements and then "
        "the parsimony score before the actual placement [EXPERIMENTAL].")(
        "sort-before-placement-3,A",
        po::bool_switch(&options.sort_by_ambiguous_bases)->default_value(false),
        "Sort new samples based on the number of ambiguous bases "
        "[EXPERIMENTAL].")(
        "reverse-sort,r",
        po::bool_switch(&options.reverse_sort)->default_value(false),
        "Reverse the sorting order of sorting options (sort-before-placement-1 "
        "or sort-before-placement-2) [EXPERIMENTAL]")(
        "collapse-tree,c",
        po::bool_switch(&options.collapse_tree)->default_value(false),
        "Collapse internal nodes of the input tree with no mutations and "
        "condense identical sequences in polytomies into a single node and the "
        "save the tree to file condensed-tree.nh in outdir")(
        "collapse-output-tree,C", po::bool_switch(&ignored_options),
        "Collapse internal nodes of the output tree with no mutations before "
        "the saving the tree to file final-tree.nh in outdir")(
        "max-uncertainty-per-sample,e",
        po::value(&options.max_uncertainty)->default_value(1e6),
        "Maximum number of equally parsimonious placements allowed per sample "
        "beyond which the sample is ignored")(
        "max-parsimony-per-sample,E",
        po::value(&options.max_parsimony)->default_value(1e6),
        "Maximum parsimony score of the most parsimonious placement(s) allowed "
        "per sample beyond which the sample is ignored")(
        "write-uncondensed-final-tree,u",
        po::bool_switch(&options.out_options.print_uncondensed_tree)
            ->default_value(false),
        "Write the final tree in uncondensed format and save to file "
        "uncondensed-final-tree.nh in outdir")(
        "write-subtrees-size,k",
        po::value<size_t>(&options.out_options.print_subtrees_size)
            ->default_value(0),
        "Write minimum set of subtrees covering the newly added samples of "
        "size equal to this value")(
        "write-single-subtree,K",
        po::value<size_t>(&options.out_options.print_subtrees_single)
            ->default_value(0),
        "Similar to write-subtrees-size but produces a single subtree with all "
        "newly added samples along with random samples up to the value "
        "specified by this argument")(
        "retain-input-branch-lengths,l",
        po::bool_switch(&options.out_options.retain_original_branch_len)
            ->default_value(false),
        "Retain the branch lengths from the input tree in out newick files "
        "instead of using number of mutations for the branch lengths.")(
        "detailed-clades,D",
        po::bool_switch(&options.out_options.detailed_clades),
        "In clades.txt, write a histogram of annotated clades and counts "
        "across all equally parsimonious placements")(
        "version", "Print version number")("help,h", "Print help messages");
    po::variables_map vm;
    // wait_debug();
    try {
        po::store(po::command_line_parser(args.size(), args.data())
                      .options(desc)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        return -1;
    }

    for (size_t idx = 1; idx < args.size(); idx++) {
        free(args[idx]);
    }
    boost::filesystem::path path(options.out_options.outdir);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(options.out_options.outdir);
    options.out_options.outdir = path.generic_string();
    options.override_mutations = false;
    return mat_idx;
}
static void child_proc(int fd, TreeCollectionPtr &trees_ptr) {
    FILE *f = fdopen(fd, "a+");
    Leader_Thread_Options options;
    auto idx = get_options(f, options);
    if (idx >= trees_ptr->size()) {
        fprintf(f, "got idx %zu but only have %zu trees\n",idx,trees_ptr->size());
        fputc(4, f);
        fputc('\n', f);
        fclose(f);
        exit(EXIT_FAILURE);
    }
    MAT::Tree &tree = (*trees_ptr)[idx].tree;
    std::vector<Sample_Muts> samples_to_place;
    tbb::task_scheduler_init init(num_threads);
    std::vector<mutated_t> position_wise_out;
    std::vector<mutated_t> position_wise_out_dup;
    std::vector<std::string> samples;
    std::unordered_set<std::string> &samples_in_condensed_nodes =
        (*trees_ptr)[idx].condensed_nodes;
    Sample_Input(options.vcf_filename.c_str(), samples_to_place, tree,
                 position_wise_out, false, samples, samples_in_condensed_nodes);
    samples_to_place.resize(
        std::min(samples_to_place.size(), options.first_n_samples));
    size_t sample_start_idx = samples_to_place[0].sample_idx;
    size_t sample_end_idx = samples_to_place.back().sample_idx + 1;
    std::vector<std::string> low_confidence_samples;
    std::vector<Clade_info> samples_clade(samples_to_place.size());
    for (auto &temp : samples_clade) {
        temp.valid = false;
    }
    if (samples_to_place.empty()) {
        fprintf(f, "Found no new samples\n");
        fputc(4, f);
        fputc('\n', f);
        fclose(f);
        exit(EXIT_FAILURE);
    }
    std::string placement_stats_filename =
        options.out_options.outdir + "/placement_stats.tsv";
    FILE *placement_stats_file = fopen(placement_stats_filename.c_str(), "w");
    //auto reordered =
        sort_samples(options, samples_to_place, tree, sample_start_idx);
    place_sample_sequential(samples_to_place, tree, false, placement_stats_file,
                            options.max_parsimony, options.max_uncertainty,
                            low_confidence_samples, samples_clade,
                            sample_start_idx, true, f);
    auto dfs = tree.depth_first_expansion();
    clean_up_leaf(dfs);
    final_output(tree, options.out_options, 0, samples_clade, sample_start_idx,
                 sample_end_idx, low_confidence_samples, position_wise_out);
    fputc(4, f);
    fputc('\n', f);
    fclose(f);
    exit(EXIT_SUCCESS);
}
static void accept_fork_loop(int socket_fd, TreeCollectionPtr &trees_ptr) {
    std::unordered_map<int, int> pid_to_fd_map;
    while (true) {
        if (interrupted) {
            break;
        }
        auto conn_fd = accept(socket_fd, NULL, 0);
        if (interrupted) {
            break;
        }
        if (conn_fd == -1) {
            perror("cannot accept connection");
            continue;
        }
        {
            TreeCollectionPtr local_copy = trees_ptr;
            auto pid = fork();
            if (pid == 0) {
                child_proc(conn_fd, trees_ptr);
            } else {
                local_copy.reset();
                pid_to_fd_map.emplace(pid, conn_fd);
            }
        }
        collect_done(pid_to_fd_map, false);
    }
    collect_done(pid_to_fd_map, true);
}

int main(int argc, char** argv){
    po::options_description desc{"Options"};
    std::string mgr_fifo;
    std::string socket_path;
    num_threads=tbb::task_scheduler_init::default_num_threads();
    std::vector<std::string> init_pb_to_load;
    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_threads) +
                                      " detected on this machine]";
    desc.add_options()("manager-fifo-path,m",po::value<std::string>(&mgr_fifo),
    "path to a fifo taking commands,existing file will be deleted currently support:\n"
    "stop: stop taking new connections and exit\n"
    "thread [int] : reset the number of threads for each job\n"
    "reload : replace the cached trees, expecting one per line\n")
    ("socket-patch,s",po::value<std::string>(&socket_path),
    "path to socket,existing file will be deleted\n"
    "Expect usher cmd args separated by new line, terminated with an empty line"
    "then the server will reply with the output of usher\n"
    "")
    ("threads-per-process,t",po::value<unsigned int>(&num_threads),num_threads_message.c_str())
    ("pb-to-load,l",po::value<std::vector<std::string>>(&init_pb_to_load)->multitoken()->composing(),"initial list of protobufs to load")
    ;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        return -1;
    }
    if (socket_path.length()>=107) {
        fprintf(stderr, "socket path length %lu is too long, cannot exceed 107 bytes\n",socket_path.length());
        exit(EXIT_FAILURE);
    }
    if (socket_path=="") {
        fprintf(stderr, "socket path is empty\n");
        exit(EXIT_FAILURE);
    }
    
    TreeCollectionPtr trees;
    if (!init_pb_to_load.empty()) {
        reload_trees(trees, init_pb_to_load);    
    }
    auto socket_fd=create_socket(socket_path);
    std::thread mgr(mgr_thread,std::ref(trees), mgr_fifo, socket_path);
    accept_fork_loop(socket_fd, trees);
    mgr.join();
    return EXIT_SUCCESS;
}