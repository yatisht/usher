#include "ripples.hpp"
#include "src/usher_graph.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <array>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <random>
#include <time.h>
#include <vector>
#define CHECK_MAPPER
Timer timer;
struct idx_hash {
    size_t operator()(const MAT::Node* in) const {
        return in->dfs_idx;
    }
};
struct idx_eq {
    bool operator()(const MAT::Node* a,const MAT::Node* b) const {
        return a->dfs_idx==b->dfs_idx;
    }
};
MAT::Node* get_node_cstr (MAT::Tree& tree,char* name) {
    return tree.get_node(std::string(name));
}
struct interval_sorter {
    bool operator()(Recomb_Interval& a,Recomb_Interval& b) {
        if (a.start_range_high<b.start_range_high) {
            return true;
        } else if(a.start_range_high==b.start_range_high&&a.end_range_low<b.end_range_low) {
            return true;
        }
        return false;
    }
};
struct Ripple_Result_Pack {
    std::vector<Recomb_Interval> intervals;
    MAT::Node* node_to_consider;
    int orig_parsimony;
};
struct next_node {
    std::vector<MAT::Node*>::const_iterator& iter;
    std::vector<MAT::Node*>::const_iterator end;
    MAT::Node* operator()(tbb::flow_control& fc) const {
        if (iter==end) {
            fc.stop();
        }
        auto to_resturn=*iter;
        iter++;
        return to_resturn;
    }
};
struct Ripple_Pipeline {
    const std::vector<MAT::Node*>& nodes_to_seach;
    const std::vector<bool>& do_parallel;
    const std::vector<int>& index_map;
    uint32_t branch_len;
    int min_range;
    int max_range;
    uint32_t num_threads;
    int parsimony_improvement;
    MAT::Tree& T;
    const std::vector<Mapper_Info> &traversal_track;
    const unsigned short tree_height;
    Ripple_Result_Pack* operator()(MAT::Node* node_to_consider) const;
};
struct Ripple_Finalizer {
    FILE *desc_file;
    FILE *recomb_file;
    size_t& num_done;
    size_t total_size;
    MAT::Tree& T;

    void operator()(Ripple_Result_Pack*) const;
};
int main(int argc, char **argv) {
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string outdir = vm["outdir"].as<std::string>();
    std::string samples_filename = vm["samples-filename"].as<std::string>();
    uint32_t branch_len = vm["branch-length"].as<uint32_t>();
    int parsimony_improvement = vm["parsimony-improvement"].as<int>();
    int max_range = vm["max-coordinate-range"].as<int>();
    int min_range = vm["min-coordinate-range"].as<int>();
    uint32_t num_descendants = vm["num-descendants"].as<uint32_t>();
    int start_idx = vm["start-index"].as<int>();
    int end_idx = vm["end-index"].as<int>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);
    srand(time(NULL));

    static tbb::affinity_partitioner ap;

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    //get_node_cstr(T,"");
    timer.Start();

    fprintf(stderr,
            "Finding the branches with number of mutations equal to or "
            "exceeding %u.\n",
            branch_len);

    tbb::concurrent_unordered_set<MAT::Node *,idx_hash,idx_eq> nodes_to_consider;
    auto dfs = T.depth_first_expansion();

    if (samples_filename != "") {
        std::ifstream infile(samples_filename);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the samples file: %s!\n",
                    samples_filename.c_str());
            exit(1);
        }
        std::string line;

        fprintf(stderr, "Reading samples from the file %s.\n",
                samples_filename.c_str());
        timer.Start();
        while (std::getline(infile, line)) {
            std::vector<std::string> words;
            MAT::string_split(line, words);
            if (words.size() != 1) {
                fprintf(stderr,
                        "ERROR: Incorrect format for samples file: %s!\n",
                        samples_filename.c_str());
                exit(1);
            }
            auto n = T.get_node(words[0]);
            if (n == NULL) {
                fprintf(stderr, "ERROR: Node id %s not found!\n",
                        words[0].c_str());
                exit(1);
            } else {
                for (auto anc : T.rsearch(n->identifier, true)) {
                    nodes_to_consider.insert(anc);
                }
                //nodes_to_consider.insert(n);
            }
        }
        infile.close();
    } else {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, dfs.size()),
        [&](const tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                auto n = dfs[i];
                if (n == T.root) {
                    continue;
                }
                if (n->mutations.size() >= branch_len) {
                    if (T.get_num_leaves(n) >= num_descendants) {
                        nodes_to_consider.insert(n);
                    }
                }
            }
        },
        ap);
    }

    std::vector<MAT::Node *> nodes_to_consider_vec;
    for (auto elem : nodes_to_consider) {
        nodes_to_consider_vec.emplace_back(elem);
    }
    std::sort(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end());
    std::shuffle(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end(),
                 std::default_random_engine(0));

    fprintf(stderr, "Found %zu long branches\n", nodes_to_consider.size());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    fprintf(stderr, "Creating output files.\n");
    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        boost::filesystem::create_directory(path);
    }

    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();

    auto desc_filename = outdir + "/descendants.tsv";
    fprintf(stderr,
            "Creating file %s to write descendants of recombinant nodes\n",
            desc_filename.c_str());
    FILE *desc_file = fopen(desc_filename.c_str(), "w");
    fprintf(desc_file, "#node_id\tdescendants\n");

    auto recomb_filename = outdir + "/recombination.tsv";
    fprintf(stderr, "Creating file %s to write recombination events\n",
            recomb_filename.c_str());
    FILE *recomb_file = fopen(recomb_filename.c_str(), "w");
    fprintf(
        recomb_file,
        "#recomb_node_id\tbreakpoint-1_interval\tbreakpoint-2_interval\tdonor_"
        "node_id\tdonor_is_sibling\tdonor_parsimony\tacceptor_node_"
        "id\tacceptor_is_sibling\tacceptor_parsimony\toriginal_parsimony\tmin_"
        "starting_parsimony\trecomb_parsimony\n");
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    std::vector<MAT::Node*> nodes_to_search;
    nodes_to_search.reserve(dfs.size());
    for (auto &node : dfs) {
        if ((node->dfs_end_idx - node->dfs_idx) >= num_descendants) {
            nodes_to_search.push_back(node);
        }
    }
    std::vector<bool> do_parallel(dfs.size(),false);
    std::vector<Mapper_Info> traversal_track;
    unsigned short tree_height=0;
    check_parallelizable(T.root,do_parallel,nodes_to_search.size()/num_threads,num_descendants,tree_height,traversal_track,0);
    std::vector<int> index_map;
    int node_to_search_idx=0;
    index_map.reserve(dfs.size());
    for (int dfs_idx = 0; dfs_idx <(int) dfs.size(); dfs_idx++) {
        if (node_to_search_idx!=(int)nodes_to_search.size()&&(int)nodes_to_search[node_to_search_idx]->dfs_idx==dfs_idx) {
            index_map.push_back(node_to_search_idx);
            node_to_search_idx++;
        } else {
            index_map.push_back(-node_to_search_idx);
        }
    }

    fprintf(stderr, "%zu out of %zu nodes have enough descendant to be donor/acceptor",nodes_to_search.size(),dfs.size());
    size_t s = 0, e = nodes_to_consider.size();

    if ((start_idx >= 0) && (end_idx >= 0)) {
        s = start_idx;
        if (end_idx <= (int)e) {
            e = end_idx;
        }
    }

    fprintf(stderr,
            "Running placement individually for %zu branches to identify "
            "potential recombination events.\n",
            e - s);

    size_t num_done = 0;
    //FILE* before_joining_fh=fopen("before_join_test","w");
    std::vector<MAT::Node*>::const_iterator cur_iter=nodes_to_consider_vec.begin()+s;
    std::vector<MAT::Node*>::const_iterator end=nodes_to_consider_vec.begin()+e;
    tbb::parallel_pipeline(
        4, tbb::make_filter<void, MAT::Node *>(tbb::filter::serial_in_order,
                next_node{cur_iter, end}) &
        tbb::make_filter<MAT::Node *, Ripple_Result_Pack *>(
            tbb::filter::parallel,
            Ripple_Pipeline{nodes_to_search, do_parallel, index_map,
                            branch_len, min_range, max_range,
                            num_threads, parsimony_improvement, T,traversal_track,tree_height}) &
        tbb::make_filter<Ripple_Result_Pack *, void>(
            tbb::filter::serial_in_order,
            Ripple_Finalizer{desc_file, recomb_file, num_done,
                             nodes_to_consider_vec.size(), T}));
    fclose(desc_file);
    fclose(recomb_file);

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
Ripple_Result_Pack* Ripple_Pipeline::operator()(MAT::Node* node_to_consider) const {
    fprintf(stderr, "At node id: %s\n",
            node_to_consider->identifier.c_str());

    int orig_parsimony = (int)node_to_consider->mutations.size();

    Pruned_Sample pruned_sample(node_to_consider);
    // Find mutations on the node to prune
    auto node_to_root = T.rsearch(node_to_consider->identifier, true);
    for (auto curr : node_to_root) {
        for (auto m : curr->mutations) {
            pruned_sample.add_mutation(m);
        }
    }
    //==== new mapper
    Ripples_Mapper_Output_Interface mapper_out;
    ripples_mapper(pruned_sample, mapper_out, nodes_to_seach.size(),index_map,do_parallel, traversal_track,tree_height,T.root,node_to_consider);
    //==== END new mapper
    tbb::concurrent_vector<Recomb_Interval> valid_pairs_con;
    ripplrs_merger(pruned_sample, index_map,nodes_to_seach, nodes_to_seach.size(),
                   orig_parsimony - parsimony_improvement, T,
                   valid_pairs_con, mapper_out, num_threads, branch_len,
                   min_range, max_range);
    std::vector<Recomb_Interval> temp(std::vector<Recomb_Interval>(valid_pairs_con.begin(),valid_pairs_con.end()));
    std::sort(temp.begin(),temp.end(),interval_sorter());
    /*       for(auto p: temp) {
               std::string end_range_high_str = (p.end_range_high == 1e9) ? "GENOME_SIZE" : std::to_string(p.end_range_high);
                           fprintf(
                   before_joining_fh,
                   "%s\t(%i,%i)\t(%i,%s)\t%s\t%s\n",
                   node_to_consider->identifier.c_str(), p.start_range_low,
                   p.start_range_high, p.end_range_low, end_range_high_str.c_str(),
                   p.d.node->identifier.c_str(), p.a.node->identifier.c_str());
               fflush(before_joining_fh);
           }*/
    return (new Ripple_Result_Pack{combine_intervals(temp),node_to_consider,orig_parsimony});
}
void Ripple_Finalizer::operator()(Ripple_Result_Pack* result) const {
    // print combined pairs
    auto & valid_pairs=result->intervals;
    auto node_to_consider=result->node_to_consider;
    auto orig_parsimony=result->orig_parsimony;
    for (auto p : valid_pairs) {
        std::string end_range_high_str =
            (p.end_range_high == 1e9) ? "GENOME_SIZE"
            : std::to_string(p.end_range_high);
        auto donor_adj_parsimony=p.d.node_parsimony+!p.d.is_sibling;
        auto acceptor_adj_parsimony=p.a.node_parsimony+!p.a.is_sibling;
        fprintf(
            recomb_file,
            "%s\t(%i,%i)\t(%i,%s)\t%s\t%c\t%i\t%s\t%c\t%i\t%i\t%i\t%i\n",
            node_to_consider->identifier.c_str(), p.start_range_low,
            p.start_range_high, p.end_range_low, end_range_high_str.c_str(),
            p.d.node->identifier.c_str(), p.d.is_sibling?'y':'n',
            donor_adj_parsimony, p.a.node->identifier.c_str(),
            p.a.is_sibling?'y':'n', acceptor_adj_parsimony, orig_parsimony,
            std::min(
        {orig_parsimony, donor_adj_parsimony, acceptor_adj_parsimony}),
        p.d.parsimony + p.a.parsimony);
        fflush(recomb_file);
    }

    if (!valid_pairs.empty()) {
        fprintf(desc_file, "%s\t", node_to_consider->identifier.c_str());
        for (auto l : T.get_leaves(node_to_consider->identifier)) {
            fprintf(desc_file, "%s,", l->identifier.c_str());
        }
        fprintf(desc_file, "\n");
        fflush(desc_file);
        fprintf(stderr, "Done %zu/%zu branches [RECOMBINATION FOUND!]\n\n",
                ++num_done, total_size);
    } else {
        fprintf(stderr, "Done %zu/%zu branches\n\n", ++num_done,
                total_size);
    }
    delete result;
}