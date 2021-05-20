#include "optimize.hpp"
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "boost/filesystem.hpp"
#include <time.h>  
#include <set>  

Timer timer;

po::variables_map check_options(int argc, char** argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::options_description desc("optimize options");
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file to optimize [REQUIRED].")
        ("output-mat,o", po::value<std::string>()->required(),
         "Output optimized mutation-annotated tree file [REQUIRED].") 
        ("sample-names,f", po::value<std::string>()->default_value(""),
         "File containing sample names whose ancestral nodes would be considered for optimization.") 
        ("input-vcf,v", po::value<std::string>()->default_value(""),
         "Input vcf filename to reassign resolved bases using Fitch-Sankoff initially.")
        ("radius,r", po::value<uint32_t>()->default_value(10), \
         "Radius in which to restrict the SPR moves.")
        ("optimization-seconds,s", po::value<uint32_t>()->default_value(360000), \
         "Approximate number of seconds to run a single iteration of the tree optimization stage. The iteration terminates as soon as the elapsed time exceeds this value.")
        ("save-every-seconds,S", po::value<uint32_t>()->default_value(300), \
         "Periodically save the optimized tree after every specified number of seconds have elapsed since the last save.") 
        ("max-iterations,N", po::value<int>()->default_value(1000), \
         "Maximum number of optimization iterations to perform.")
        ("min-improvement,m", po::value<float>()->default_value(0.001), \
         "Minimum improvement in the parsimony score as a fraction of the previous score in ordder to perform another iteration.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(p)
                .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
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
            mut.par_nuc = mut.ref_nuc;
            sample_mutations.insert(iter, mut);
        }
        positions.insert(mut.position);
    }

    Pruned_Sample (std::string name) {
        sample_name = name;
        sample_mutations.clear();
        positions.clear();
    }
};


inline std::string get_reversal_mutation(const std::string& mut_string) {
    std::vector<std::string> words;
    MAT::string_split(mut_string, words);
    assert(words.size() == 4);
    return words[0] + "\t" + words[3] + "\t" + words[2] + "\t" + words[1]; 
}

size_t get_node_distance (const MAT::Tree& T, MAT::Node* source, MAT::Node* dest) {
    auto s_anc = T.rsearch(source->identifier, true);
    auto d_anc = T.rsearch(dest->identifier, true);

    size_t s_off=0, d_off=0, num_iter=0;
    if (s_anc.size() <= d_anc.size()) {
        d_off = d_anc.size() - s_anc.size();
        num_iter = s_anc.size();
    }
    else {
        s_off = s_anc.size() - d_anc.size();
        num_iter = d_anc.size();
    }
    
    auto distance = s_off+d_off;
    for (size_t i=0; i<num_iter; i++) {
        if (s_anc[s_off+i] == d_anc[d_off+i]) {
            break;
        }
        else {
            distance+=2;
        }
    }

    return distance;
}

/*
void optimize_main_old(po::parsed_options parsed) {
    po::variables_map vm = parse_optimize_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    uint32_t max_seconds = vm["optimization-seconds"].as<uint32_t>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    tbb::task_scheduler_init init(num_threads);
    srand (time(NULL));
    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "The parsimony score for this MAT is %zu\n", T.get_parsimony_score());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    timer.Start();
    fprintf(stderr, "Starting tree optimization.\n\n"); 
    
    fprintf(stderr, "Finding the nodes with recurrent or reversal mutations.\n");
    auto dfs = T.depth_first_expansion();
    
    std::map<std::string, int> mutation_counts;
    tbb::mutex tbb_lock;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
            [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
            auto n = dfs[i];
            n->clear_annotations();
            
            for (auto m: n->mutations) {
               std::string mut_string = m.chrom + "\t" + std::to_string(m.par_nuc) + "\t" + 
                   std::to_string(m.position) + "\t" + std::to_string(m.mut_nuc);
               tbb_lock.lock();
               if (mutation_counts.find(mut_string) == mutation_counts.end()) {
                   mutation_counts[mut_string] = 1;
               }
               else {
                   mutation_counts[mut_string] += 1;
               }
               tbb_lock.unlock();
           }
        }
    }, ap);
    std::vector<std::string> nodes_to_prune;
    for (auto n: dfs) {
        if ((n == T.root) || (T.get_num_leaves(n) > 1000)) {
            continue;
        }
        for (auto m: n->mutations) {
            std::string mut_string = m.chrom + "\t" + std::to_string(m.par_nuc) + "\t" + 
                std::to_string(m.position) + "\t" + std::to_string(m.mut_nuc);
            // Contains recurrent mutation
            if (mutation_counts[mut_string] > 1) {
                nodes_to_prune.emplace_back(n->identifier);
                break;
            }
            auto reversal_string = get_reversal_mutation(mut_string);
            // reversal mutation found
            if (mutation_counts.find(reversal_string) != mutation_counts.end()) {
                nodes_to_prune.emplace_back(n->identifier);
                break;
            }
        }
    }
    // Shuffle the nodes to be pruned
    std::random_shuffle(nodes_to_prune.begin(), nodes_to_prune.end());
    fprintf(stderr, "%zu nodes found with recurrent or reversal mutations.\n\n", nodes_to_prune.size());
    auto best_parsimony_score = T.get_parsimony_score();
    for (auto node_to_prune: nodes_to_prune) {
        if (T.get_node(node_to_prune) == NULL) {
            continue;
        }
        auto copy = MAT::get_tree_copy(T);
        fprintf(stderr, "Pruning subtree at node %s.\n", node_to_prune.c_str());
        auto leaves_to_prune = T.get_leaves(node_to_prune);
        std::vector<Pruned_Sample> pruned_samples;
        for (auto l: leaves_to_prune) {
            Pruned_Sample to_prune(l->identifier);
            auto root_to_node = T.rsearch(l->identifier, true); 
            std::reverse(root_to_node.begin(), root_to_node.end());
            for (auto curr: root_to_node) {
                for (auto m: curr->mutations) {
                    to_prune.add_mutation(m);
                }
            }
            // move_level set to false as moving it would be an unnecessary overhead
            T.remove_node(l->identifier, false);
            pruned_samples.emplace_back(to_prune);
        }
        fprintf(stderr, "Pruned %zu samples at node %s.\n", pruned_samples.size(), node_to_prune.c_str());
        fprintf(stderr, "Placing %zu samples at node %s in random order.\n", pruned_samples.size(), node_to_prune.c_str());
        std::random_shuffle(pruned_samples.begin(), pruned_samples.end());
        for (auto s: pruned_samples) {
            dfs = T.depth_first_expansion();
            size_t total_nodes = dfs.size();
            // Stores the excess mutations to place the sample at each
            // node of the tree in DFS order. When placement is as a
            // child, it only contains parsimony-increasing mutations in
            // the sample. When placement is as a sibling, it contains 
            // parsimony-increasing mutations as well as the mutations
            // on the placed node in common with the new sample. Note
            // guaranteed to be corrrect only for optimal nodes since
            // the mapper can terminate the search early for non-optimal
            // nodes
            std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
            // Stores the imputed mutations for ambiguous bases in the
            // sampled in order to place the sample at each node of the 
            // tree in DFS order. Again, guaranteed to be corrrect only 
            // for pasrimony-optimal nodes 
            std::vector<std::vector<MAT::Mutation>> node_imputed_mutations(total_nodes);
            std::vector<int> node_set_difference;
            size_t best_node_num_leaves = 0;
            // The maximum number of mutations is bound by the number
            // of mutations in the missing sample (place at root)
            //int best_set_difference = 1e9;
            // TODO: currently number of root mutations is also added to
            // this value since it forces placement as child but this
            // could be changed later 
            int best_set_difference = s.sample_mutations.size() + T.root->mutations.size() + 1;
            size_t best_j = 0;
            bool best_node_has_unique = false;
            std::vector<bool> node_has_unique(total_nodes, false);
            std::vector<size_t> best_j_vec;
            size_t num_best = 1;
            MAT::Node* best_node = T.root;
            best_j_vec.emplace_back(0);
            // Parallel for loop to search for most parsimonious
            // placements. Real action happens within mapper2_body
            static tbb::affinity_partitioner ap;
            tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k){
                    mapper2_input inp;
                    inp.T = &T;
                    inp.node = dfs[k];
                    inp.missing_sample_mutations = &s.sample_mutations;
                    inp.excess_mutations = &node_excess_mutations[k];
                    inp.imputed_mutations = &node_imputed_mutations[k];
                    inp.best_node_num_leaves = &best_node_num_leaves;
                    inp.best_set_difference = &best_set_difference;
                    inp.best_node = &best_node;
                    inp.best_j =  &best_j;
                    inp.num_best = &num_best;
                    inp.j = k;
                    inp.has_unique = &best_node_has_unique;
                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);
                    mapper2_body(inp, false);
                    }       
                    }, ap); 
            // Ensure sample not already in the tree
            if (T.get_node(s.sample_name) == NULL) {
                // Is placement as sibling
                if (best_node->is_leaf() || best_node_has_unique) {
                    std::string nid = std::to_string(++T.curr_internal_node);
                    T.create_node(nid, best_node->parent->identifier);
                    T.create_node(s.sample_name, nid);
                    T.move_node(best_node->identifier, nid);
                    // common_mut stores mutations common to the
                    // best node branch and the sample, l1_mut
                    // stores mutations unique to best node branch
                    // and l2_mut stores mutations unique to the
                    // sample not in best node branch
                    std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                    std::vector<MAT::Mutation> curr_l1_mut;
                    // Compute current best node branch mutations
                    for (auto m1: best_node->mutations) {
                        MAT::Mutation m = m1.copy();
                        curr_l1_mut.emplace_back(m);
                    }
                    // Clear mutations on the best node branch which
                    // will be later replaced by l1_mut
                    best_node->clear_mutations();
                    // Compute l1_mut
                    for (auto m1: curr_l1_mut) {
                        bool found = false;
                        for (auto m2: node_excess_mutations[best_j]) {
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
                            MAT::Mutation m = m1.copy();
                            l1_mut.emplace_back(m);
                        }
                    }
                    // Compute l2_mut
                    for (auto m1: node_excess_mutations[best_j]) {
                        bool found = false;
                        for (auto m2: curr_l1_mut) {
                            if (m1.is_masked()) {
                                break;
                            }
                            if (m1.position == m2.position) {
                                if (m1.mut_nuc == m2.mut_nuc) {
                                    found = true;
                                    MAT::Mutation m = m1.copy();
                                    common_mut.emplace_back(m);
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            MAT::Mutation m = m1.copy();
                            l2_mut.emplace_back(m);
                        }
                    }
                    // Add mutations to new node using common_mut
                    for (auto m: common_mut) {
                        T.get_node(nid)->add_mutation(m);
                    }
                    // Add mutations to best node using l1_mut
                    for (auto m: l1_mut) {
                        T.get_node(best_node->identifier)->add_mutation(m);
                    }
                    // Add new sample mutations using l2_mut
                    for (auto m: l2_mut) {
                        T.get_node(s.sample_name)->add_mutation(m);
                    }
                }
                // Else placement as child
                else {
                    MAT::Node* node = T.create_node(s.sample_name, best_node->identifier);
                    std::vector<MAT::Mutation> node_mut;
                    std::vector<MAT::Mutation> curr_l1_mut;
                    for (auto m1: best_node->mutations) {
                        MAT::Mutation m = m1.copy();
                        curr_l1_mut.emplace_back(m);
                    }
                    for (auto m1: node_excess_mutations[best_j]) {
                        bool found = false;
                        for (auto m2: curr_l1_mut) {
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
                            MAT::Mutation m = m1.copy();
                            node_mut.emplace_back(m);
                        }
                    }
                    for (auto m: node_mut) {
                        node->add_mutation(m);
                    }
                }
            }
        }
        fprintf(stderr, "Done pruning and placing samples at node %s\n", node_to_prune.c_str());
        auto new_parsimony_score = T.get_parsimony_score();
        if (new_parsimony_score >= best_parsimony_score) {
            fprintf(stderr, "Placement resulted in a parsimony score of %zu. Retaining old tree.\n", new_parsimony_score);
            MAT::clear_tree(T);
            T = copy;
        }
        else {
            MAT::clear_tree(copy);
            best_parsimony_score = new_parsimony_score; 
            fprintf(stderr, "Placement lowered parsimony score to %zu!\n", best_parsimony_score);
        }
        auto elapsed_time = timer.Stop()/1000;
        fprintf(stderr, "Elapsed time (tree optimization stage): %ld seconds.\n\n", elapsed_time);
        if (elapsed_time >= max_seconds) {
            break;
        }
    }
    
    fprintf(stderr, "Optimization complete. Saving the final MAT with a parsimony score of %zu\n", T.get_parsimony_score());
    T.collapse_tree();
    T.condense_leaves();
    MAT::save_mutation_annotated_tree(T, output_mat_filename);
    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
*/

int main(int argc, char** argv) {
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string input_vcf_filename = vm["input-vcf"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string sample_name_filename = vm["sample-names"].as<std::string>();
    uint32_t max_seconds = vm["optimization-seconds"].as<uint32_t>();
    uint32_t save_every = vm["save-every-seconds"].as<uint32_t>();
    uint32_t radius = vm["radius"].as<uint32_t>();
    int max_iterations = vm["max-iterations"].as<int>();
    float min_improvement = vm["min-improvement"].as<float>();

    uint32_t num_threads = vm["threads"].as<uint32_t>();


    if (min_improvement > 1.0) {
        fprintf(stderr, "ERROR! The value of --min-improvement should be less than or equal to 1.0.\n");
        exit(1);
    }

    tbb::task_scheduler_init init(num_threads);
    srand (time(NULL));

    static tbb::affinity_partitioner ap;

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    
    std::unordered_map<std::string, bool> nid_has_changed;
    float last_improvement = 1.0;

    for (int iter=0; iter < max_iterations; iter++) {
        T.uncondense_leaves();
    
        // If VCF specified, re-compute the assignments using Fitch-Sankoff
        if (input_vcf_filename != "") {
            fprintf(stderr, "Computing parsimonious assignments for input variants.\n"); 
            timer.Start();

            // Variables below used to store the different fields of the input VCF file 
            bool header_found = false;
            std::vector<std::string> variant_ids;
            std::vector<Missing_Sample> missing_samples;

            // Vector used to store all tree nodes in breadth-first search (BFS) order
            std::vector<MAT::Node*> bfs;
            // Map the node identifier string to index in the BFS traversal
            std::unordered_map<std::string, size_t> bfs_idx;

            // Breadth-first expansion to populate bfs and bfs_idx
            bfs = T.breadth_first_expansion();
            for (size_t idx = 0; idx < bfs.size(); idx++) {
                bfs_idx[bfs[idx]->identifier] = idx;
                // clear the node mutations
                bfs[idx]->clear_mutations();
            }

            // Boost library used to stream the contents of the input VCF file in
            // uncompressed or compressed .gz format
            std::ifstream infile(input_vcf_filename, std::ios_base::in | std::ios_base::binary);
            if (!infile) {
                fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", input_vcf_filename.c_str());
                exit(1);
            }
            boost::iostreams::filtering_istream instream;
            try {
                if (input_vcf_filename.find(".gz\0") != std::string::npos) {
                    instream.push(boost::iostreams::gzip_decompressor());
                }
                instream.push(infile);
            }
            catch(const boost::iostreams::gzip_error& e) {
                std::cout << e.what() << '\n';
            }


            // A TBB flow graph containing a single source_node (reader) connected
            // to several mappers. The source_node sequentially reads in the different 
            // lines of the input VCF file and constructs a mapper task for each
            // VCF line. Each mapper task takes a mapper_input as input, which stores
            // the alternate alleles, ambiguous bases and missing data (Ns) for
            // different tree samples at the corresponding VCF line/position. The 
            // mappers use Fitch-Sankoff algorithm to assign mutations at different
            // branches of the tree and update the mutation-annotated tree (T)
            // accordingly. 
            tbb::flow::graph mapper_graph;

            tbb::flow::function_node<mapper_input, int> mapper(mapper_graph, tbb::flow::unlimited, mapper_body());
            tbb::flow::source_node <mapper_input> reader (mapper_graph,
                    [&] (mapper_input &inp) -> bool {

                    //check if reached end-of-file
                    int curr_char = instream.peek();
                    if(curr_char == EOF)
                    return false;

                    std::string s;
                    std::getline(instream, s);
                    std::vector<std::string> words;
                    MAT::string_split(s, words);
                    inp.variant_pos = -1;

                    // Header is found when "POS" is the second word in the line
                    if ((not header_found) && (words.size() > 1)) {
                    if (words[1] == "POS") {
                    // Sample names start from the 10th word in the header
                    for (size_t j=9; j < words.size(); j++) {
                    variant_ids.emplace_back(words[j]);
                    // If sample name not in tree, add it to missing_samples
                    if (bfs_idx.find(words[j]) == bfs_idx.end()) {
                        missing_samples.emplace_back(Missing_Sample(words[j]));
                    }
                    }
                    header_found = true;
                    }
                    }
                    else if (header_found) {
                        if (words.size() != 9+variant_ids.size()) {
                            fprintf(stderr, "ERROR! Incorrect VCF format.\n");
                            exit(1);
                        }
                        std::vector<std::string> alleles;
                        alleles.clear();
                        inp.variant_pos = std::stoi(words[1]); 
                        MAT::string_split(words[4], ',', alleles);
                        // T will be modified by the mapper with mutation
                        // annotations
                        inp.T = &T;
                        inp.chrom = words[0];
                        inp.bfs = &bfs;
                        inp.bfs_idx = &bfs_idx;
                        inp.variant_ids = &variant_ids;
                        inp.missing_samples = &missing_samples;
                        // Ref nuc id uses one-hot encoding (A:0b1, C:0b10, G:0b100,
                        // T:0b1000)
                        inp.ref_nuc = MAT::get_nuc_id(words[3][0]);
                        assert((inp.ref_nuc & (inp.ref_nuc-1)) == 0); //check if it is power of 2
                        inp.variants.clear();
                        for (size_t j=9; j < words.size(); j++) {
                            if (isdigit(words[j][0])) {
                                int allele_id = std::stoi(words[j]);
                                if (allele_id > 0) { 
                                    std::string allele = alleles[allele_id-1];
                                    inp.variants.emplace_back(std::make_tuple(j-9, MAT::get_nuc_id(allele[0])));
                                }
                            }
                            else {
                                inp.variants.emplace_back(std::make_tuple(j-9, MAT::get_nuc_id('N')));
                            }
                        }
                    }
                    return true;
                    }, true );
            tbb::flow::make_edge(reader, mapper);
            mapper_graph.wait_for_all();

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }

        // Collapse tree for optimal performance and results
        T.collapse_tree();

        fprintf(stderr, "The parsimony score for this MAT is %zu\n", T.get_parsimony_score());
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        float last_update = 0;
        float curr_parsimony = (float) T.get_parsimony_score();

        timer.Start();
        fprintf(stderr, "Finding the subtrees to prune.\n");
        auto bfs = T.breadth_first_expansion();

        if (nid_has_changed.size() == 0) {
            for (auto n: bfs) {
                nid_has_changed.insert({n->identifier, true});
            }
        }

        struct Profitable_Node {
            std::string nid;
            int profit;

            Profitable_Node (std::string id, int p) {
                nid = id;
                profit = p;
            }

            inline bool operator<(const Profitable_Node& n) {
                  return (*this).nid < n.nid;
            }
            
            inline bool operator==(const Profitable_Node& n) {
                  return (*this).nid == n.nid;
            }
        };
        
        std::vector<Profitable_Node> profitable_nodes;

        if (sample_name_filename != "") {
            std::ifstream infile(sample_name_filename);
            if (!infile) {
                fprintf(stderr, "ERROR: Could not open the sample name file: %s!\n", sample_name_filename.c_str());
                exit(1);
            }    
            std::string line;

            fprintf(stderr, "Reading sample name file.\n"); 
            timer.Start();
            while (std::getline(infile, line)) {
                std::vector<std::string> words;
                MAT::string_split(line, words);
                if (words.size() != 1) {
                    fprintf(stderr, "ERROR: Incorrect format for sample name file: %s!\n", sample_name_filename.c_str());
                    exit(1);
                }
                auto n = T.get_node(words[0]);
                if (n == NULL) {
                    fprintf(stderr, "WARNING: Node id %s not found in the tree!\n", words[0].c_str());
                    exit(1);
                }
                else {
                    for (auto anc: T.rsearch(words[0], true)) {
                        if (anc->level >= 2) { 
                            profitable_nodes.emplace_back(Profitable_Node(anc->identifier, 0));
                        }
                    }
                }
            }
            infile.close();
        }

        else {
            std::unordered_map<int, std::vector<std::string>> pos_to_nid;
            for (auto n: bfs) {
                for (auto m: n->mutations) {
                    if (pos_to_nid.find(m.position) != pos_to_nid.end()) {
                        pos_to_nid[m.position].emplace_back(n->identifier);
                    }
                    else {
                        pos_to_nid[m.position] = std::vector<std::string>();
                        pos_to_nid[m.position].emplace_back(n->identifier);
                    }
                }
            }

            size_t at = 0;
            size_t total = 0;
            for (auto d: pos_to_nid) {
                total+=d.second.size();
            }
            tbb::mutex tbb_lock;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, pos_to_nid.size()),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t it = r.begin(); it < r.end(); it++) {
                    auto cn = pos_to_nid.begin();
                    std::advance(cn, it);
                    size_t num_elem = cn->second.size();
                    for (size_t i=0; i<num_elem; i++) {
                    __sync_fetch_and_add(&at, 1);
                    tbb_lock.lock();
                    fprintf(stderr, "\r%zu/%zu", at, total); 
//                    if (std::find(profitable_nodes.begin(), profitable_nodes.end(), 
//                                Profitable_Node(cn->second[i],0)) != profitable_nodes.end()) {
//                            tbb_lock.unlock();
//                            continue;
//                    }
                    tbb_lock.unlock();
                    auto n = T.get_node(cn->second[i]);
                    if (n->level < 2) {
                        continue;
                    }

                    for (auto m: n->mutations) {
                        if (nid_has_changed[n->identifier] && (m.ref_nuc == m.mut_nuc)) {
                            tbb_lock.lock();
                            profitable_nodes.emplace_back(Profitable_Node(n->identifier, 0));
                            tbb_lock.unlock();
                        }
                    }

                    int max_profit = 0;

                    //bool has_inserted = false;
                    for (size_t j=0; j<num_elem; j++) {
                        auto n2 = T.get_node(cn->second[j]);
                        if (n == n2) {
                            continue;
                        }
                        if ((nid_has_changed[n2->identifier] || nid_has_changed[n->identifier]) && 
                               //(!T.is_ancestor(n->identifier, n2->identifier)) && 
                               (get_node_distance(T, n, n2) <= radius)) {
                            
                           // Find mutations on the node to prune
                           Pruned_Sample pruned_sample(n->identifier);

                           auto node_to_root = T.rsearch(n->identifier, true); 
                           for (auto curr: node_to_root) {
                               for (auto m: curr->mutations) {
                                   pruned_sample.add_mutation(m);
                               }
                           }

                           std::vector<MAT::Mutation> node_excess_mutations;
                           std::vector<MAT::Mutation> imputed_mutations;

                           size_t best_node_num_leaves = 0;
                           int best_set_difference = 1e9;

                           std::vector<bool> node_has_unique(1);
                           size_t best_j = 0;
                           bool best_node_has_unique = false;

                           size_t best_distance = 1e9;

                           std::vector<size_t> best_j_vec;

                           size_t num_best = 1;
                           MAT::Node* best_node = NULL;
                           best_j_vec.emplace_back(0);

                           size_t node_distance=0;

                           mapper2_input inp;
                           inp.T = &T;
                           inp.node = n2;
                           inp.missing_sample_mutations = &pruned_sample.sample_mutations;
                           inp.excess_mutations = &node_excess_mutations;
                           inp.imputed_mutations = &imputed_mutations;
                           inp.best_node_num_leaves = &best_node_num_leaves;
                           inp.best_set_difference = &best_set_difference;
                           inp.best_node = &best_node;
                           inp.best_j =  &best_j;
                           inp.num_best = &num_best;
                           inp.j = 0;
                           inp.has_unique = &best_node_has_unique;

                           inp.distance = node_distance;
                           inp.best_distance = &best_distance;

                           inp.best_j_vec = &best_j_vec;
                           inp.node_has_unique = &(node_has_unique);

                           mapper2_body(inp, false);
                           
                           int profit = n->mutations.size() - best_set_difference;
                           if (profit > max_profit) {
                               max_profit = profit;
                           }

                        }
                    }
                    
                    if (max_profit > 0) {
                        tbb_lock.lock();
                        profitable_nodes.emplace_back(Profitable_Node(n->identifier, max_profit));
                        //profitable_nodes.emplace_back(Profitable_Node(n2->identifier, max_profit));
                        tbb_lock.unlock();
                    }
                    }
                    }
                    }, ap);
            fprintf(stderr, "\n"); 

        }

        tbb::parallel_sort(profitable_nodes.begin(), profitable_nodes.end(), 
                [](const Profitable_Node& n1, const Profitable_Node& n2) {return n1.profit > n2.profit;}); 


        std::vector<Profitable_Node> unique_profitable_nodes;
        std::set<std::string> nodes_to_prune;

        for (auto pn: profitable_nodes) {
            if (nodes_to_prune.find(pn.nid) == nodes_to_prune.end()) {
                unique_profitable_nodes.emplace_back(pn);
                nodes_to_prune.insert(pn.nid);
            }
        }
        
        fprintf(stderr, "%zu subtrees found to prune.\n", nodes_to_prune.size());
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        timer.Start();
        fprintf(stderr, "Starting tree optimization.\n\n"); 

        auto best_parsimony_score = T.get_parsimony_score();

        size_t curr=0;
        for (auto elem: unique_profitable_nodes) {
            auto nid_to_prune = elem.nid;

            fprintf(stderr, "At %zu of %zu\n", ++curr, nodes_to_prune.size());

            auto node_to_prune = T.get_node(nid_to_prune);
            if (node_to_prune == NULL) {
                continue;
            }

            fprintf(stderr, "Pruning subtree at node %s.\n", nid_to_prune.c_str());

            // Find mutations on the node to prune
            Pruned_Sample pruned_sample(nid_to_prune);

            auto node_to_root = T.rsearch(nid_to_prune, true); 
            //std::reverse(root_to_node.begin(), root_to_node.end());
            for (auto curr: node_to_root) {
                for (auto m: curr->mutations) {
                    pruned_sample.add_mutation(m);
                }
            }

            // Now prune the node_to_prune from the tree
            auto curr_parent = node_to_prune->parent; 
            auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), node_to_prune);
            assert (iter != curr_parent->children.end());
            curr_parent->children.erase(iter);

            // Set source to current parent
            auto source = curr_parent;

            // Remove curr_parent if it has no children
            if (curr_parent->children.size() == 0) {
                auto ancestors = T.rsearch(curr_parent->identifier, true);
                T.remove_node(curr_parent->identifier, true);

                // Since remove_node can remove multiple levels of ancestors,
                // reassign source to nearest ancestor currently in the tree
                for (auto anc: ancestors) {
                    if (T.get_node(anc->identifier) != NULL) {
                        source = anc;
                        break;
                    }
                }
            }

            // Move the remaining child one level up if it is the only child of its parent 
            if (curr_parent->children.size() == 1) {
                auto child = curr_parent->children[0];
                if (curr_parent->parent != NULL) {
                    child->parent = curr_parent->parent;
                    // Adjust levels
                    for (auto n: T.breadth_first_expansion(child->identifier)) {
                        if (n == T.root) {
                            n->level = 1;
                        }
                        else {
                            n->level = n->parent->level+1;
                        }
                    }
                    for (auto m: curr_parent->mutations) {
                        child->add_mutation(m.copy());
                    }

                    curr_parent->parent->children.push_back(child);

                    curr_parent->children.clear();
                    T.remove_node(curr_parent->identifier, false);

                    source = child;
                }
            }

            bfs = T.breadth_first_expansion();
            size_t total_nodes = bfs.size();

            std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
            std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);

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

            std::vector<size_t> node_distance(total_nodes);

            static tbb::affinity_partitioner ap;
            tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k){
                    node_distance[k] = get_node_distance(T, source, bfs[k]);
                    if (node_distance[k] > radius) {
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

                    inp.distance = node_distance[k];
                    inp.best_distance = &best_distance;

                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);

                    mapper2_body(inp, false);
                    }       
                    }, ap); 

            auto distance = node_distance[best_j]; 

            // Clear current mutations of the node to prune
            node_to_prune->clear_mutations();

            // Is placement as sibling
            if (best_node->is_leaf() || best_node_has_unique) {
                std::string nid = std::to_string(++T.curr_internal_node);
                auto new_node = T.create_node(nid, best_node->parent->identifier);
                node_to_prune->parent = new_node;
                new_node->children.emplace_back(node_to_prune);
                T.move_node(best_node->identifier, nid);

                std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                std::vector<MAT::Mutation> curr_l1_mut;

                // Compute current best node branch mutations
                for (auto m1: best_node->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }
                // Clear mutations on the best node branch which
                // will be later replaced by l1_mut
                best_node->clear_mutations();

                // Compute l1_mut
                for (auto m1: curr_l1_mut) {
                    bool found = false;
                    for (auto m2: node_excess_mutations[best_j]) {
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
                        MAT::Mutation m = m1.copy();
                        l1_mut.emplace_back(m);
                    }
                }
                // Compute l2_mut
                for (auto m1: node_excess_mutations[best_j]) {
                    bool found = false;
                    for (auto m2: curr_l1_mut) {
                        if (m1.is_masked()) {
                            break;
                        }
                        if (m1.position == m2.position) {
                            if (m1.mut_nuc == m2.mut_nuc) {
                                found = true;
                                MAT::Mutation m = m1.copy();
                                common_mut.emplace_back(m);
                                break;
                            }
                        }
                    }
                    if (!found) {
                        MAT::Mutation m = m1.copy();
                        l2_mut.emplace_back(m);
                    }
                }

                // Add mutations to new node using common_mut
                for (auto m: common_mut) {
                    new_node->add_mutation(m);
                }
                // Add mutations to best node using l1_mut
                for (auto m: l1_mut) {
                    best_node->add_mutation(m);
                }
                // Add new sample mutations using l2_mut
                for (auto m: l2_mut) {
                    node_to_prune->add_mutation(m);
                }
            }
            // Else placement as child
            else {
                best_node->children.emplace_back(node_to_prune);
                node_to_prune->parent = best_node;

                std::vector<MAT::Mutation> node_mut;

                std::vector<MAT::Mutation> curr_l1_mut;

                for (auto m1: best_node->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }

                for (auto m1: node_excess_mutations[best_j]) {
                    bool found = false;
                    for (auto m2: curr_l1_mut) {
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
                        MAT::Mutation m = m1.copy();
                        node_mut.emplace_back(m);
                    }
                }
                for (auto m: node_mut) {
                    node_to_prune->add_mutation(m);
                }
            }

            fprintf(stderr, "Done pruning and placing samples at node %s\n", best_node->identifier.c_str());
            fprintf(stderr, "Distance from original placement: %zu\n", distance);

            auto new_parsimony_score = T.get_parsimony_score();
            assert(new_parsimony_score <= best_parsimony_score);
            if (new_parsimony_score == best_parsimony_score) {
                fprintf(stderr, "Parsimony score of %zu unchanged after placement.\n", new_parsimony_score);
                nid_has_changed[node_to_prune->identifier] = false;
            }
            else {
                best_parsimony_score = new_parsimony_score; 
                fprintf(stderr, "Placement lowered parsimony score to %zu!\n", best_parsimony_score);
                nid_has_changed[node_to_prune->identifier] = true;
                nid_has_changed[source->identifier] = true;
            }

            float elapsed_sec = static_cast<float>(timer.Stop())/1000;
            fprintf(stderr, "Elapsed time (current optimization iteration): %.3f sec\n\n", elapsed_sec);

            if (elapsed_sec >= max_seconds) {
                break;
            }

            if (elapsed_sec-last_update >= save_every) {
                // Levels need to be adjusted after placement 
                bfs = T.depth_first_expansion();
                for (auto n: bfs) {
                    if (n == T.root) {
                        n->level = 1;
                    }
                    else {
                        n->level = n->parent->level+1;
                    }
                }

                fprintf(stderr, "Saving the current MAT with a parsimony score of %zu\n\n", T.get_parsimony_score());
                T.condense_leaves();
                MAT::save_mutation_annotated_tree(T, output_mat_filename);
                last_update = elapsed_sec;

                // Uncondense for the next iteration
                T.uncondense_leaves();
            }
        }

        // Levels need to be adjusted after placement 
        bfs = T.breadth_first_expansion();
        for (auto n: bfs) {
            if (n == T.root) {
                n->level = 1;
            }
            else {
                n->level = n->parent->level+1;
            }
        }

        fprintf(stderr, "Iteration %i complete. Saving the final MAT with a parsimony score of %zu\n", iter+1, T.get_parsimony_score());
        T.condense_leaves();
        MAT::save_mutation_annotated_tree(T, output_mat_filename);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        // Quit if parsimony score not improving by minimum threshold
        float new_parsimony = (float) T.get_parsimony_score();
        last_improvement = (curr_parsimony - new_parsimony)/curr_parsimony;
        if (last_improvement <= min_improvement) {
            fprintf(stderr, "\nLast iteration improvement (%f) smaller than minimum required improvement (%f). Quitting.\n\n", last_improvement, min_improvement);
            break;
        }
    }
}

