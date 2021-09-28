#include <time.h>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "tbb/concurrent_unordered_set.h"
#include "../usher_graph.hpp"


namespace po = boost::program_options;

Timer timer;

po::variables_map check_options(int argc, char** argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::options_description desc("optimize options");
    desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file to optimize [REQUIRED].")
    ("branch-length,l", po::value<uint32_t>()->default_value(3), \
     "Minimum length of the branch to consider to recombination events")
    ("min-coordinate-range,r", po::value<int>()->default_value(1e3), \
     "Minimum range of the genomic coordinates of the mutations on the recombinant branch")
    ("max-coordinate-range,R", po::value<int>()->default_value(1e7), \
     "Maximum range of the genomic coordinates of the mutations on the recombinant branch")
    ("outdir,d", po::value<std::string>()->default_value("."),
     "Output directory to dump output files [DEFAULT uses current directory]")
    ("samples-filename,s", po::value<std::string>()->default_value(""),
     "Restrict the search to the ancestors of the samples specified in the input file")

    ("parsimony-improvement,p", po::value<int>()->default_value(3), \
     "Minimum improvement in parsimony score of the recombinant sequence during the partial placement")
    ("num-descendants,n", po::value<uint32_t>()->default_value(10), \
     "Minimum number of leaves that node should have to be considered for recombination.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("start-index,S", po::value<int>()->default_value(-1), "start index [EXPERIMENTAL]")
    ("end-index,E", po::value<int>()->default_value(-1), "end index [EXPERIMENTAL]")
    ("help,h", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(p)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless
        // the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

struct Pruned_Sample {
    std::string sample_name;
    std::vector<MAT::Mutation> sample_mutations;
    std::unordered_set<uint32_t> positions;

    // Assumes mutations are added in reverse chrono order
    void add_mutation (MAT::Mutation mut) {
        // If not reversal to reference allele
        if ((mut.ref_nuc != mut.mut_nuc) && (positions.find(mut.position) == positions.end())) {
            auto iter = std::lower_bound(sample_mutations.begin(), sample_mutations.end(), mut);
            auto m = mut.copy();
            m.par_nuc = m.ref_nuc;
            sample_mutations.insert(iter, m);
        }
        positions.insert(mut.position);
    }

    Pruned_Sample (std::string name) {
        sample_name = name;
        sample_mutations.clear();
        positions.clear();
    }
};

struct Recomb_Node {
    std::string name;
    int node_parsimony;
    int parsimony;
    char is_sibling;
    Recomb_Node() {
        name = "";
        node_parsimony = -1;
        parsimony = -1;
        is_sibling = false;
    }
    Recomb_Node(std::string n, int np, int p, char s) : name(n), node_parsimony(np), parsimony(p), is_sibling(s)
    {}
    inline bool operator< (const Recomb_Node& n) const {
        return (((*this).parsimony < n.parsimony) || (((*this).name <  n.name) && ((*this).parsimony == n.parsimony)));
    }

};

struct Recomb_Interval {
    Recomb_Node d;//donor
    Recomb_Node a;//acceptor
    int start_range_low;
    int start_range_high;
    int end_range_low;
    int end_range_high;
    Recomb_Interval(Recomb_Node donor, Recomb_Node acceptor, int sl, int sh, int el, int eh):
        d(donor), a(acceptor), start_range_low(sl), start_range_high(sh), end_range_low(el), end_range_high(eh)
    {}
    bool operator<(const Recomb_Interval& other) const { //compare second interval
        return end_range_low < other.end_range_low;
    }
};

struct Comp_First_Interval {
    inline bool operator()(const Recomb_Interval& one, const Recomb_Interval& other) {
        return one.start_range_low < other.start_range_low;
    }
};


std::vector<Recomb_Interval> combine_intervals(std::vector<Recomb_Interval> pair_list) {
    //combine second interval
    std::vector<Recomb_Interval> pairs(pair_list);
    std::sort(pairs.begin(), pairs.end());//sorts by beginning of second interval
    for(size_t i = 0; i < pairs.size(); i++) {
        for(size_t j = i+1; j < pairs.size(); j++) {
            //check everything except first interval is same and first interval of pairs[i] ends where it starts for pairs[j]
            if((pairs[i].d.name == pairs[j].d.name) && (pairs[i].a.name == pairs[j].a.name) &&
                    (pairs[i].start_range_low == pairs[j].start_range_low) && (pairs[i].start_range_high == pairs[j].start_range_high) &&
                    (pairs[i].end_range_high == pairs[j].end_range_low) && ((pairs[i].d.parsimony+pairs[i].a.parsimony) == (pairs[j].d.parsimony+pairs[j].a.parsimony))) {
                pairs[i].end_range_high = pairs[j].end_range_high;
                pairs.erase(pairs.begin()+j);//remove the combined element
                j--;
            }
        }
    }
    //combine first interval
    std::sort(pairs.begin(), pairs.end(), Comp_First_Interval());
    for(size_t i = 0; i < pairs.size(); i++) {
        for(size_t j = i + 1; j < pairs.size(); j++) {
            //check everything except second interval is same and second interval of pairs[i] ends where it starts for pairs[j]
            if((pairs[i].d.name == pairs[j].d.name) && (pairs[i].a.name == pairs[j].a.name) &&
                    (pairs[i].end_range_low == pairs[j].end_range_low) && (pairs[i].end_range_high == pairs[j].end_range_high) &&
                    (pairs[i].start_range_high == pairs[j].start_range_low) && ((pairs[i].d.parsimony+pairs[i].a.parsimony) == (pairs[j].d.parsimony+pairs[j].a.parsimony))) {
                pairs[i].start_range_high = pairs[j].start_range_high;
                pairs.erase(pairs.begin()+j);
                j--;
            }
        }
    }
    return pairs;
}


int main(int argc, char** argv) {
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
    srand (time(NULL));

    static tbb::affinity_partitioner ap;

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();

    fprintf(stderr, "Finding the branches with number of mutations equal to or exceeding %u.\n", branch_len);
    auto bfs = T.breadth_first_expansion();

    tbb::concurrent_unordered_set<std::string> nodes_to_consider;

    if (samples_filename != "") {
        std::ifstream infile(samples_filename);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the samples file: %s!\n", samples_filename.c_str());
            exit(1);
        }
        std::string line;

        fprintf(stderr, "Reading samples from the file %s.\n", samples_filename.c_str());
        timer.Start();
        while (std::getline(infile, line)) {
            std::vector<std::string> words;
            MAT::string_split(line, words);
            if (words.size() != 1) {
                fprintf(stderr, "ERROR: Incorrect format for samples file: %s!\n", samples_filename.c_str());
                exit(1);
            }
            auto n = T.get_node(words[0]);
            if (n == NULL) {
                fprintf(stderr, "ERROR: Node id %s not found!\n", words[0].c_str());
                exit(1);
            } else {
                for (auto anc: T.rsearch(n->identifier, true)) {
                    nodes_to_consider.insert(anc->identifier);
                }
            }
        }
        infile.close();
    } else {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, bfs.size()),
        [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i) {
                auto n = bfs[i];
                if (n == T.root) {
                    continue;
                }
                if (n->mutations.size() >= branch_len) {
                    if (T.get_num_leaves(n) >= num_descendants) {
                        nodes_to_consider.insert(n->identifier);
                    }
                }
            }
        }, ap);
    }

    std::vector<std::string> nodes_to_consider_vec;
    for (auto elem: nodes_to_consider) {
        nodes_to_consider_vec.emplace_back(elem);
    }
    std::sort(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end());
    std::shuffle(nodes_to_consider_vec.begin(), nodes_to_consider_vec.end(), std::default_random_engine(0));

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
    fprintf(stderr, "Creating file %s to write descendants of recombinant nodes\n",
            desc_filename.c_str());
    FILE* desc_file = fopen(desc_filename.c_str(), "w");
    fprintf(desc_file, "#node_id\tdescendants\n");

    auto recomb_filename = outdir + "/recombination.tsv";
    fprintf(stderr, "Creating file %s to write recombination events\n",
            recomb_filename.c_str());
    FILE* recomb_file = fopen(recomb_filename.c_str(), "w");
    fprintf(recomb_file, "#recomb_node_id\tbreakpoint-1_interval\tbreakpoint-2_interval\tdonor_node_id\tdonor_is_sibling\tdonor_parsimony\tacceptor_node_id\tacceptor_is_sibling\tacceptor_parsimony\toriginal_parsimony\tmin_starting_parsimony\trecomb_parsimony\n");
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();

    std::unordered_map<MAT::Node*, size_t> tree_num_leaves;

    for (int i=int(bfs.size())-1; i>=0; i--) {
        auto n = bfs[i];
        size_t desc = 1;
        for (auto child: n->children) {
            desc += tree_num_leaves[child];
        }
        tree_num_leaves[n] = desc;
    }

    size_t s = 0, e = nodes_to_consider.size();

    if ((start_idx >= 0) && (end_idx >= 0)) {
        s = start_idx;
        if (end_idx <= (int) e) {
            e = end_idx;
        }
    }

    fprintf(stderr, "Running placement individually for %zu branches to identify potential recombination events.\n", e-s);

    size_t num_done = 0;
    for (size_t idx = s; idx < e; idx++) {
        auto nid_to_consider = nodes_to_consider_vec[idx];
        fprintf(stderr, "At node id: %s\n", nid_to_consider.c_str());

        int orig_parsimony = (int) T.get_node(nid_to_consider)->mutations.size();

        Pruned_Sample pruned_sample(nid_to_consider);
        // Find mutations on the node to prune
        auto node_to_root = T.rsearch(nid_to_consider, true);
        for (auto curr: node_to_root) {
            for (auto m: curr->mutations) {
                pruned_sample.add_mutation(m);
            }
        }
        size_t num_mutations = pruned_sample.sample_mutations.size();

        size_t total_nodes = bfs.size();

        std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
        std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);

        std::vector<int> node_set_difference(total_nodes);

        size_t best_node_num_leaves = 0;
        int best_set_difference = 1e9;

        std::vector<bool> node_has_unique(total_nodes);
        size_t best_j = 0;
        bool best_node_has_unique = false;

        size_t best_distance = 1e9;

        std::vector<size_t> best_j_vec;

        size_t num_best = 1;
        MAT::Node* best_node = T.root;
        best_j_vec.emplace_back(0);

        std::vector<size_t> node_distance(total_nodes, 0);

        tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k) {
                if (tree_num_leaves[bfs[k]] < num_descendants) {
                    continue;
                }

                node_has_unique[k] = false;

                mapper2_input inp;
                inp.T = &T;
                inp.node = bfs[k];
                inp.missing_sample_mutations = &pruned_sample.sample_mutations;
                inp.excess_mutations = &node_excess_mutations[k];
                inp.imputed_mutations = &imputed_mutations[k];
                inp.best_node_num_leaves = &best_node_num_leaves;
                inp.best_set_difference = &best_set_difference;
                inp.best_node = &best_node;
                inp.best_j =  &best_j;
                inp.num_best = &num_best;
                inp.j = k;
                inp.has_unique = &best_node_has_unique;

                inp.set_difference = &node_set_difference[k];

                inp.distance = node_distance[k];
                inp.best_distance = &best_distance;

                inp.best_j_vec = &best_j_vec;
                inp.node_has_unique = &(node_has_unique);

                mapper2_body(inp, true);

            }
        }, ap);

        std::vector<Recomb_Interval> valid_pairs;
        bool has_recomb = false;

        size_t at = 0;
        for (size_t i = 0; i<num_mutations; i++) {
            for (size_t j=i; j<num_mutations; j++) {
                fprintf(stderr, "\rTrying %zu of %zu breakpoint pairs.", ++at, num_mutations*(1+num_mutations)/2);
                Pruned_Sample donor("donor");
                Pruned_Sample acceptor("acceptor");

                donor.sample_mutations.clear();
                acceptor.sample_mutations.clear();

                for (size_t k=0; k<num_mutations; k++) {
                    if ((k>=i) && (k<j)) {
                        donor.add_mutation(pruned_sample.sample_mutations[k]);
                    } else {
                        acceptor.add_mutation(pruned_sample.sample_mutations[k]);
                    }
                }

                int start_range_high = pruned_sample.sample_mutations[i].position;
                int start_range_low = (i>=1) ? pruned_sample.sample_mutations[i-1].position : 0;

                //int end_range_high = pruned_sample.sample_mutations[j].position;
                int end_range_high = 1e9;
                int end_range_low = (j>=1) ? pruned_sample.sample_mutations[j-1].position : 0;

                //int end_range_low = pruned_sample.sample_mutations[j].position;
                //int end_range_high = (j+1<num_mutations) ? pruned_sample.sample_mutations[j+1].position : 1e9;

                if ((donor.sample_mutations.size() < branch_len) || (acceptor.sample_mutations.size() < branch_len) ||
                        (end_range_low-start_range_high < min_range) || (end_range_low-start_range_high > max_range)) {
                    continue;
                }

                tbb::concurrent_vector<Recomb_Node> donor_nodes;
                tbb::concurrent_vector<Recomb_Node> acceptor_nodes;

                donor_nodes.clear();
                acceptor_nodes.clear();
                // find acceptor(s)
                {
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {


                        for (size_t k=r.begin(); k<r.end(); ++k) {
                            size_t num_mut = 0;
                            if ((tree_num_leaves[bfs[k]] < num_descendants) || (T.is_ancestor(nid_to_consider,bfs[k]->identifier))) {
                                continue;
                            }
                            // Is placement as sibling
                            if (bfs[k]->is_leaf() || node_has_unique[k]) {
                                std::vector<MAT::Mutation> l2_mut;

                                // Compute l2_mut
                                for (auto m1: node_excess_mutations[k]) {
                                    bool found = false;
                                    for (auto m2: bfs[k]->mutations) {
                                        if (m1.is_masked()) {
                                            break;
                                        }
                                        if (m1.position == m2.position) {
                                            if (m1.mut_nuc == m2.mut_nuc) {
                                                found = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!found) {
                                        if ((m1.position < start_range_high) || (m1.position > end_range_low)) {
                                            num_mut++;
                                        }
                                    }
                                    if (num_mut + parsimony_improvement > (size_t) orig_parsimony) {
                                        break;
                                    }
                                }

                                if (num_mut + parsimony_improvement <= (size_t) orig_parsimony) {
                                    acceptor_nodes.emplace_back(Recomb_Node(bfs[k]->identifier, node_set_difference[k], num_mut, 'y'));
                                }
                            }
                            // Else placement as child
                            else {
                                for (auto m1: node_excess_mutations[k]) {
                                    bool found = false;
                                    for (auto m2: bfs[k]->mutations) {
                                        if (m1.is_masked()) {
                                            break;
                                        }
                                        if (m1.position == m2.position) {
                                            if (m1.mut_nuc == m2.mut_nuc) {
                                                found = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!found) {
                                        if ((m1.position < start_range_high) || (m1.position > end_range_low)) {
                                            num_mut++;
                                        }
                                    }
                                    if (num_mut + parsimony_improvement > (size_t) orig_parsimony) {
                                        break;
                                    }
                                }

                                if (num_mut + parsimony_improvement <= (size_t) orig_parsimony) {
                                    acceptor_nodes.emplace_back(Recomb_Node(bfs[k]->identifier, node_set_difference[k], num_mut, 'n'));
                                }
                            }
                        }
                    }, ap);
                }

                if (acceptor_nodes.size() == 0) {
                    continue;
                }

                // find donor(s)
                {
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {
                        for (size_t k=r.begin(); k<r.end(); ++k) {
                            if ((tree_num_leaves[bfs[k]] < num_descendants) || (T.is_ancestor(nid_to_consider,bfs[k]->identifier))) {
                                continue;
                            }

                            size_t num_mut = 0;

                            // Is placement as sibling
                            if (bfs[k]->is_leaf() || node_has_unique[k]) {
                                std::vector<MAT::Mutation> l2_mut;

                                // Compute l2_mut
                                for (auto m1: node_excess_mutations[k]) {
                                    bool found = false;
                                    for (auto m2: bfs[k]->mutations) {
                                        if (m1.is_masked()) {
                                            break;
                                        }
                                        if (m1.position == m2.position) {
                                            if (m1.mut_nuc == m2.mut_nuc) {
                                                found = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!found) {
                                        if ((m1.position >= start_range_high) && (m1.position <= end_range_low)) {
                                            num_mut++;
                                        }
                                    }

                                    if (num_mut + parsimony_improvement > (size_t) orig_parsimony) {
                                        break;
                                    }
                                }

                                if (num_mut + parsimony_improvement <= (size_t) orig_parsimony) {
                                    donor_nodes.emplace_back(Recomb_Node(bfs[k]->identifier, node_set_difference[k], num_mut, 'y'));
                                }
                            }
                            // Else placement as child
                            else {
                                for (auto m1: node_excess_mutations[k]) {
                                    bool found = false;
                                    for (auto m2: bfs[k]->mutations) {
                                        if (m1.is_masked()) {
                                            break;
                                        }
                                        if (m1.position == m2.position) {
                                            if (m1.mut_nuc == m2.mut_nuc) {
                                                found = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!found) {
                                        if ((m1.position >= start_range_high) && (m1.position <= end_range_low)) {
                                            num_mut++;
                                        }
                                    }
                                    if (num_mut + parsimony_improvement > (size_t) orig_parsimony) {
                                        break;
                                    }
                                }

                                if (num_mut + parsimony_improvement <= (size_t) orig_parsimony) {
                                    donor_nodes.emplace_back(Recomb_Node(bfs[k]->identifier, node_set_difference[k], num_mut, 'n'));
                                }
                            }
                        }
                    }, ap);
                }

                tbb::parallel_sort (donor_nodes.begin(), donor_nodes.end());
                tbb::parallel_sort (acceptor_nodes.begin(), acceptor_nodes.end());


                if (donor_nodes.size() > 1000) {
                    donor_nodes.resize(1000);
                }

                if (acceptor_nodes.size() > 1000) {
                    acceptor_nodes.resize(1000);
                }


                // to print any pair of breakpoint interval exactly once for
                // multiple donor-acceptor pairs
                bool has_printed = false;

                for (auto d: donor_nodes) {
                    if (T.is_ancestor(nid_to_consider,d.name)) {
                        continue;
                    }
                    for (auto a:acceptor_nodes) {
                        if (T.is_ancestor(nid_to_consider,a.name)) {
                            continue;
                        }
                        // Ensure donor and acceptor are not the same and
                        // neither of them is a descendant of the recombinant
                        // node total parsimony is less than the maximum allowed
                        if ((d.name!=a.name) && (d.name!=nid_to_consider) && (a.name!=nid_to_consider) &&
                                //(!T.is_ancestor(nid_to_consider,d.name)) && (!T.is_ancestor(nid_to_consider,a.name)) &&
                                (orig_parsimony >= d.parsimony + a.parsimony + parsimony_improvement)) {
                            Pruned_Sample donor("curr-donor");
                            donor.sample_mutations.clear();
                            acceptor.sample_mutations.clear();

                            for (auto anc: T.rsearch(d.name, true)) {
                                for (auto mut: anc->mutations) {
                                    donor.add_mutation(mut);
                                }
                            }

                            for (auto mut: donor.sample_mutations) {
                                if ((mut.position > start_range_low) && (mut.position <= start_range_high)) {
                                    bool in_pruned_sample = false;
                                    for (auto mut2: pruned_sample.sample_mutations) {
                                        if (mut.position == mut2.position) {
                                            in_pruned_sample = true;
                                        }
                                    }
                                    if (!in_pruned_sample) {
                                        start_range_low = mut.position;
                                    }
                                }
                                if ((mut.position > end_range_low) && (mut.position <= end_range_high)) {
                                    bool in_pruned_sample = false;
                                    for (auto mut2: pruned_sample.sample_mutations) {
                                        if (mut.position == mut2.position) {
                                            in_pruned_sample = true;
                                        }
                                    }
                                    if (!in_pruned_sample) {
                                        end_range_high = mut.position;
                                    }
                                }
                            }

                            for (auto mut: pruned_sample.sample_mutations) {
                                if ((mut.position > start_range_low) && (mut.position <= start_range_high)) {
                                    bool in_pruned_sample = false;
                                    for (auto mut2: donor.sample_mutations) {
                                        if (mut.position == mut2.position) {
                                            in_pruned_sample = true;
                                        }
                                    }
                                    if (!in_pruned_sample) {
                                        start_range_low = mut.position;
                                    }
                                }
                                if ((mut.position > end_range_low) && (mut.position <= end_range_high)) {
                                    bool in_pruned_sample = false;
                                    for (auto mut2: donor.sample_mutations) {
                                        if (mut.position == mut2.position) {
                                            in_pruned_sample = true;
                                        }
                                    }
                                    if (!in_pruned_sample) {
                                        end_range_high = mut.position;
                                    }
                                }
                            }

                            //tbb_lock.lock();
                            valid_pairs.push_back(Recomb_Interval(d, a, start_range_low, start_range_high, end_range_low, end_range_high));
                            //tbb_lock.unlock();

                            has_recomb = true;
                            has_printed = true;
                            break;
                        }
                    }
                    if (has_printed) {
                        break;
                    }
                }
            }
        }
        fprintf(stderr, "\n");

        valid_pairs = combine_intervals(valid_pairs);
        //print combined pairs
        for(auto p: valid_pairs) {
            std::string end_range_high_str = (p.end_range_high == 1e9) ? "GENOME_SIZE" : std::to_string(p.end_range_high);
            fprintf(recomb_file, "%s\t(%i,%i)\t(%i,%s)\t%s\t%c\t%i\t%s\t%c\t%i\t%i\t%i\t%i\n", nid_to_consider.c_str(), p.start_range_low,
                    p.start_range_high, p.end_range_low, end_range_high_str.c_str(), p.d.name.c_str(), p.d.is_sibling, p.d.node_parsimony,
                    p.a.name.c_str(), p.a.is_sibling, p.a.node_parsimony, orig_parsimony,
                    std::min({orig_parsimony, p.d.node_parsimony, p.a.node_parsimony}), p.d.parsimony+p.a.parsimony);
            fflush(recomb_file);
        }

        if (has_recomb) {
            fprintf(desc_file, "%s\t", nid_to_consider.c_str());
            for (auto l: T.get_leaves(nid_to_consider)) {
                fprintf(desc_file, "%s,", l->identifier.c_str());
            }
            fprintf(desc_file, "\n");
            fflush(desc_file);
            fprintf(stderr, "Done %zu/%zu branches [RECOMBINATION FOUND!]\n\n", ++num_done, nodes_to_consider.size());
        } else {
            fprintf(stderr, "Done %zu/%zu branches\n\n", ++num_done, nodes_to_consider.size());
        }
    }

    fclose(desc_file);
    fclose(recomb_file);

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

