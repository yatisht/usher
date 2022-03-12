#include "ripples.hpp"
#include <algorithm>
#include <cstdint>

void Pruned_Sample::add_mutation(MAT::Mutation mut) {
    // If not reversal to reference allele
    if ((mut.ref_nuc != mut.mut_nuc) &&
            (positions.find(mut.position) == positions.end())) {
        auto iter = std::lower_bound(sample_mutations.begin(),
                                     sample_mutations.end(), mut);
        auto m = mut.copy();
        m.par_nuc = m.ref_nuc;
        sample_mutations.insert(iter, m);
    }
    positions.insert(mut.position);
}
Pruned_Sample::Pruned_Sample(MAT::Node* name) {
    sample_name = name;
    sample_mutations.clear();
    positions.clear();
}

std::vector<Recomb_Interval>
combine_intervals(std::vector<Recomb_Interval> pair_list) {
    // combine second interval
    std::vector<Recomb_Interval> pairs(pair_list);
    std::sort(pairs.begin(),
              pairs.end()); // sorts by beginning of second interval
    for (size_t i = 0; i < pairs.size(); i++) {
        for (size_t j = i + 1; j < pairs.size(); j++) {
            // check everything except first interval is same and first interval
            // of pairs[i] ends where it starts for pairs[j]
            if ((pairs[i].d.node->identifier == pairs[j].d.node->identifier) &&
                    (pairs[i].a.node->identifier == pairs[j].a.node->identifier) &&
                    (pairs[i].start_range_low == pairs[j].start_range_low) &&
                    (pairs[i].start_range_high == pairs[j].start_range_high) &&
                    (pairs[i].end_range_high == pairs[j].end_range_low) &&
                    ((pairs[i].d.parsimony + pairs[i].a.parsimony) ==
                     (pairs[j].d.parsimony + pairs[j].a.parsimony))) {
                pairs[i].end_range_high = pairs[j].end_range_high;
                pairs.erase(pairs.begin() + j); // remove the combined element
                j--;
            }
        }
    }
    // combine first interval
    std::sort(pairs.begin(), pairs.end(), Comp_First_Interval());
    for (size_t i = 0; i < pairs.size(); i++) {
        for (size_t j = i + 1; j < pairs.size(); j++) {
            // check everything except second interval is same and second
            // interval of pairs[i] ends where it starts for pairs[j]
            if ((pairs[i].d.node->identifier == pairs[j].d.node->identifier) &&
                    (pairs[i].a.node->identifier == pairs[j].a.node->identifier) &&
                    (pairs[i].end_range_low == pairs[j].end_range_low) &&
                    (pairs[i].end_range_high == pairs[j].end_range_high) &&
                    (pairs[i].start_range_high == pairs[j].start_range_low) &&
                    ((pairs[i].d.parsimony + pairs[i].a.parsimony) ==
                     (pairs[j].d.parsimony + pairs[j].a.parsimony))) {
                pairs[i].start_range_high = pairs[j].start_range_high;
                pairs.erase(pairs.begin() + j);
                j--;
            }
        }
    }
    return pairs;
}
po::variables_map check_options(int argc, char **argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_cores) +
                                      " detected on this machine]";

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
                  .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
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
size_t check_parallelizable(const MAT::Node *root,
                            std::vector<bool> &do_parallel,
                            size_t parallel_threshold,
                            size_t check_threshold,
                            unsigned short& tree_height,
                            std::vector<Mapper_Info>& traversal_track,unsigned short level) {
    size_t child_counted_size = 0;
    if ((root->dfs_end_idx - root->dfs_idx) >= check_threshold) {
        child_counted_size++;
        auto cur_idx=traversal_track.size();
        traversal_track.push_back(
            Mapper_Info{root->mutations.data(),
                        root->mutations.data() + root->mutations.size(),
                        UINT32_MAX, level, root->children.empty()});
        for (const auto child : root->children) {
            child_counted_size += check_parallelizable(
                                      child, do_parallel, parallel_threshold, check_threshold,tree_height,traversal_track,level+1);
        }
        if (child_counted_size > parallel_threshold) {
            do_parallel[root->dfs_idx] = true;
        }
        traversal_track[cur_idx].sibling_start_idx=traversal_track.size();
        tree_height=std::max(tree_height,level);
    }
    return child_counted_size;
}