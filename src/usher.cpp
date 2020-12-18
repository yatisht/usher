#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <memory>
#include <limits>
#include "boost/filesystem.hpp"
#include "usher_graph.hpp"
#include "parsimony.pb.h"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

int main(int argc, char** argv) {

    //Variables to load command-line options using Boost program_options
    std::string tree_filename;
    std::string din_filename;
    std::string dout_filename;
    std::string outdir;
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    uint32_t max_trees;
    uint32_t max_uncertainty;
    bool sort_before_placement_1 = false;
    bool sort_before_placement_2 = false;
    bool reverse_sort = false;
    bool collapse_tree=false;
    bool print_uncondensed_tree = false;
    bool print_parsimony_scores = false;
    bool retain_original_branch_len = false;
    size_t print_subtrees_size=0;
    po::options_description desc{"Options"};

    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
        ("vcf,v", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
        ("tree,t", po::value<std::string>(&tree_filename)->default_value(""), "Input tree file")
        ("outdir,d", po::value<std::string>(&outdir)->default_value("."), "Output directory to dump output and log files [DEFAULT uses current directory]")
        ("load-mutation-annotated-tree,i", po::value<std::string>(&din_filename)->default_value(""), "Load mutation-annotated tree object")
        ("save-mutation-annotated-tree,o", po::value<std::string>(&dout_filename)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
        ("sort-before-placement-1,s", po::bool_switch(&sort_before_placement_1), \
         "Sort new samples based on computed parsimony score and then number of optimal placements before the actual placement [EXPERIMENTAL].")
        ("sort-before-placement-2,S", po::bool_switch(&sort_before_placement_2), \
         "Sort new samples based on the number of optimal placements and then the parsimony score before the actual placement [EXPERIMENTAL].")
        ("reverse-sort,r", po::bool_switch(&reverse_sort), \
         "Reverse the sorting order of sorting options (sort-before-placement-1 or sort-before-placement-2) [EXPERIMENTAL]")
        ("collapse-tree,c", po::bool_switch(&collapse_tree), \
         "Collapse internal nodes of the input tree with no mutations and condense identical sequences in polytomies into a single node and the save the tree to file condensed-tree.nh in outdir")
        ("max-uncertainty-per-sample,e", po::value<uint32_t>(&max_uncertainty)->default_value(1e6), \
         "Maximum number of equally parsimonious placements allowed per sample beyond which the sample is ignored")
        ("write-uncondensed-final-tree,u", po::bool_switch(&print_uncondensed_tree), "Write the final tree in uncondensed format and save to file uncondensed-final-tree.nh in outdir")
        ("write-subtrees-size,k", po::value<size_t>(&print_subtrees_size)->default_value(0), \
         "Write minimum set of subtrees covering the newly added samples of size equal to or larger than this value")
        ("write-parsimony-scores-per-node,p", po::bool_switch(&print_parsimony_scores), \
         "Write the parsimony scores for adding new samples at each existing node in the tree without modifying the tree in a file names parsimony-scores.tsv in outdir")
        ("multiple-placements,M", po::value<uint32_t>(&max_trees)->default_value(1), \
         "Create a new tree up to this limit for each possibility of parsimony-optimal placement")
        ("retain-input-branch-lengths,l", po::bool_switch(&retain_original_branch_len), \
         "Retain the branch lengths from the input tree in out newick files instead of using number of mutations for the branch lengths.")
        ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    if (print_subtrees_size == 1) {
        std::cerr << "ERROR: print-subtrees-size should be larger than 1\n";
        return 1;
    }

    if (sort_before_placement_1 && sort_before_placement_2) {
        std::cerr << "ERROR: Can't use sort-before-placement-1 and sort-before-placement-2 simultaneously. Please specify only one.\n";
        return 1;
    }

    if (sort_before_placement_1 || sort_before_placement_2) { 
        std::cerr << "WARNING: Using experimental option ";
        if (sort_before_placement_1) {
            std::cerr << "--sort-before-placement-1 (-s)\n";
        }
        if (sort_before_placement_2) {
            std::cerr << "--sort-before-placement-2 (-S)\n";
        }
    }
    else if (reverse_sort) {
        std::cerr << "ERROR: Can't use reverse-sort without sorting options (sort-before-placement-1 or sort-before-placement-2)\n";
        return 1;
    }

    if (print_parsimony_scores) {
        if (max_trees > 1) {
            std::cerr << "ERROR: cannot use --multiple-placements (-M) and --print_parsimony_scores (-p) options simulaneously.\n";
            return 1;
        }
    
        if (sort_before_placement_1 || sort_before_placement_2 || collapse_tree || print_uncondensed_tree || (print_subtrees_size > 0) || (dout_filename != "")) {
            fprintf (stderr, "WARNING: --print-parsimony-scores-per-node is set. Will terminate without modifying the original tree.\n");
        }
    }

    if (max_trees == 0) {
        std::cerr << "ERROR: Number of trees specified by --multiple-placements (-M) should be >= 1\n";
        return 1;
    }
    /*
    if (max_trees > 1) {
        if (max_trees < 256) {
            std::cerr << "WARNING: Using experimental option --multiple-placements (-M)\n";
        }
        else { 
            std::cerr << "ERROR: Number of trees specified by --multiple-placements (-M) should be <= 255\n";
            return 1; 
        }
    }
    */

    if (retain_original_branch_len) {
        fprintf(stderr, "Output newick files will retain branch lengths from the input tree (unspecified at branches modified during the placement).\n\n");
    }
    else {
        fprintf(stderr, "Output newick files will have branch lengths equal to the number of mutations of that branch.\n\n");
    }

    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();

    // timer object to be used to measure runtimes of individual stages
    Timer timer; 

    fprintf(stderr, "Initializing %u worker threads.\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);

#if SAVE_PROFILE == 1
    Instrumentor::Get().BeginSession("test-main", "p1.json");
#endif

    // Vector to store multiple trees, each corresponding to a different
    // possibility of a  parsimony-optimal placement, when --multiple-placements
    // is used. Otherwise, this vector maintains a single tree througout the
    // execution in which a tie-breaking strategy defined in usher_mapper is
    // used for multiple parsimony-optimal placements.
    std::vector<MAT::Tree> optimal_trees;
    // Tree pointer to point to some element in optimal_trees that would be
    // updated several times during the execution 
    MAT::Tree* T = NULL;

    // Since --multiple-placements can result in trees with different parsimony
    // scores, the vector below will be used to maintain the final parsimony
    // score of each tree 
    std::vector<size_t> tree_parsimony_scores;
    
    auto num_trees = optimal_trees.size();

    // Variables below used to store the different fields of the input VCF file 
    bool header_found = false;
    std::vector<std::string> variant_ids;
    std::vector<std::string> missing_samples;
    std::vector<std::vector<MAT::Mutation>> missing_sample_mutations;
    size_t num_missing = 0;

    // Vector used to store all tree nodes in breadth-first search (BFS) order
    std::vector<MAT::Node*> bfs;
    // Map the node identifier string to index in the BFS traversal
    std::unordered_map<std::string, size_t> bfs_idx;
    
    // Vectore to store the names of samples which have a high number of 
    // parsimony-optimal placements
    std::vector<std::string> low_confidence_samples;

    // If tree filename is specified, UShER needs to first load the tree and
    // a create a new mutation-annotated tree object (stored in optimal_trees) 
    // using the sample variants in the input VCF file. If the VCF contains
    // samples missing in the input tree, they get added to missing_samples
    if (tree_filename != "") {
        // Create a new tree from the input newick file and store it in
        // optimal_trees
        auto tmp_T = MAT::create_tree_from_newick(tree_filename);
        optimal_trees.emplace_back(std::move(tmp_T));

        // Nodes of T don't have mutations assigned yet
        T = &optimal_trees[0];

        // Breadth-first expansion to populate bfs and bfs_idx
        bfs = T->breadth_first_expansion();
        for (size_t idx = 0; idx < bfs.size(); idx++) {
            bfs_idx[bfs[idx]->identifier] = idx;
        }

        // Boost library used to stream the contents of the input VCF file in
        // uncompressed or compressed .gz format
        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", vcf_filename.c_str());
            exit(1);
        }
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        }
        catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }

        fprintf(stderr, "Computing parsimonious assignments for input variants.\n"); 
        timer.Start();

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
                        missing_samples.emplace_back(words[j]);
                        num_missing++;
                      }
                    }
                    missing_sample_mutations.resize(num_missing);
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
                    inp.T = T;
                    inp.chrom = words[0];
                    inp.bfs = &bfs;
                    inp.bfs_idx = &bfs_idx;
                    inp.variant_ids = &variant_ids;
                    inp.missing_samples = &missing_samples;
                    inp.missing_sample_mutations = &missing_sample_mutations;
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

                //check if reached end-of-file
                int curr_char = instream.peek();
                if(curr_char == EOF)
                    return false;
                else
                    return true;
                }, true );
        tbb::flow::make_edge(reader, mapper);
        mapper_graph.wait_for_all();

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    // Check if input mutation-annotated tree filename is specified
    else if (din_filename != "") {
        timer.Start();

        fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", din_filename.c_str());
       
        // Load mutation-annotated tree and store it in optimal_trees
        auto tmp_T = MAT::load_mutation_annotated_tree(din_filename);
        optimal_trees.emplace_back(std::move(tmp_T));

        T = &optimal_trees[0];
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        fprintf(stderr, "Loading VCF file\n");
        timer.Start();

        // Boost library used to stream the contents of the input VCF file in
        // uncompressed or compressed .gz format
        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", vcf_filename.c_str());
            exit(1);
        }
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        }
        catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }

        std::vector<size_t> missing_idx;
        std::string s;
        // This while loop reads the VCF file line by line and populates
        // missing_samples and missing_sample_mutations based on the names and 
        // variants of missing samples. If a sample name in the VCF is already
        // found in the tree, it gets ignored with a warning message
        while (instream.peek() != EOF) {
            std::getline(instream, s);
            std::vector<std::string> words;
            MAT::string_split(s, words);
            if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                    for (size_t j=9; j < words.size(); j++) {
                        variant_ids.emplace_back(words[j]);
                        if ((T->get_node(words[j]) == NULL) && (T->condensed_leaves.find(words[j]) == T->condensed_leaves.end())) {
                            missing_samples.emplace_back(words[j]);
                            num_missing++;
                            missing_idx.emplace_back(j);
                        }
                        else {
                            fprintf(stderr, "WARNING: Ignoring sample %s as it is already in the tree.\n", words[j].c_str());
                        }
                    }
                    missing_sample_mutations.resize(num_missing);
                    header_found = true;
                }
            }
            else if (header_found) {
                if (words.size() != 9+variant_ids.size()) {
                    fprintf(stderr, "ERROR! Incorrect VCF format. Expected %zu columns but got %zu.\n", 9+variant_ids.size(), words.size());
                    exit(1);
                }
                std::vector<std::string> alleles;
                alleles.clear();
                MAT::string_split(words[4], ',', alleles);
                for (size_t k = 0; k < missing_idx.size(); k++) {
                    size_t j = missing_idx[k];
                    auto iter = missing_samples.begin();
                    std::advance(iter, k);
                    if (iter != missing_samples.end()) {
                        auto mutations_iter = missing_sample_mutations.begin() + (iter - missing_samples.begin());
                        MAT::Mutation m;
                        m.chrom = words[0];
                        m.position = std::stoi(words[1]);
                        m.ref_nuc = MAT::get_nuc_id(words[3][0]);
                        assert((m.ref_nuc & (m.ref_nuc-1)) == 0); //check if it is power of 2
                        m.par_nuc = m.ref_nuc;
                        // Alleles such as '.' should be treated as missing
                        // data. if the word is numeric, it is an index to one
                        // of the alleles 
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) { 
                                std::string allele = alleles[allele_id-1];
                                if (allele[0] == 'N') {
                                    m.is_missing = true;
                                    m.mut_nuc = MAT::get_nuc_id('N');
                                }
                                else {
                                    auto nuc = MAT::get_nuc_id(allele[0]);
                                    if (nuc == MAT::get_nuc_id('N')) {
                                        m.is_missing = true;
                                    }
                                    else {
                                        m.is_missing = false;
                                    }
                                    m.mut_nuc = nuc;
                                }
                                (*mutations_iter).emplace_back(m);
                            }
                        }
                        else {
                            m.is_missing = true;
                            m.mut_nuc = MAT::get_nuc_id('N');
                            (*mutations_iter).emplace_back(m);
                        }
                    }
                }
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    else {
        fprintf(stderr, "Error! No input tree or assignment file provided!\n");
        exit(1);
    }

    // Collapses the tree nodes not carrying a mutation and also condenses
    // identical sequences into a single node. 
    if (collapse_tree) {
        timer.Start();
        
        fprintf(stderr, "Collapsing input tree.\n");
        
        assert (optimal_trees.size() == 1);

        T->collapse_tree();

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        fprintf(stderr, "Condensing identical sequences. \n");

        T->condense_leaves();

        auto condensed_tree_filename = outdir + "/condensed-tree.nh";
        fprintf(stderr, "Writing condensed input tree to file %s\n", condensed_tree_filename.c_str());
        
        FILE* condensed_tree_file = fopen(condensed_tree_filename.c_str(), "w");
        fprintf(condensed_tree_file, "%s\n", MAT::get_newick_string(*T, true, true, retain_original_branch_len).c_str());
        fclose(condensed_tree_file);

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    fprintf(stderr, "Found %zu missing samples.\n\n", missing_samples.size()); 

    FILE* parsimony_scores_file = NULL; 
        
    // If samples found in VCF that are missing from the input tree, they are
    // now placed using maximum parsimony or if print_parsimony_score is set,
    // the parsimony scores for placning the sample is printed (without the
    // actual placement). 
    if (missing_samples.size() > 0) {
        
        // indexes stores the order in which the missing samples should be
        // placed sequentially on the tree. It is initialized with sequentially
        // increasing values 0,1,2,.. which will be modified if one of the
        // sorting options is used
        std::vector<size_t> indexes(missing_samples.size());
        std::iota(indexes.begin(), indexes.end(), 0);
        
        // Write current tree with internal nodes labelled. Parsimony score of
        // placement at each node (including internal nodes) will be printed later.
        if (print_parsimony_scores) {
            timer.Start();
            auto current_tree_filename = outdir + "/current-tree.nh";
            
            fprintf(stderr, "Writing current tree with internal nodes labelled to file %s \n", current_tree_filename.c_str());
            FILE* current_tree_file = fopen(current_tree_filename.c_str(), "w");
            fprintf(current_tree_file, "%s\n", MAT::get_newick_string(*T, true, true, retain_original_branch_len).c_str());
            fclose(current_tree_file);
            
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        } 
        else {
            if ((sort_before_placement_1 || sort_before_placement_2) && (missing_samples.size() > 1)) {
                timer.Start();
                fprintf(stderr, "Computing parsimony scores and number of parsimony-optimal placements for new samples and using them to sort the samples.\n"); 
                if (max_trees > 1) {
                    fprintf(stderr, "WARNING: --multiple-placements option is used but note that the samples will be sorted only once using the parsimony scores on input tree (without actual placements).\n"); 
                }

                // vectors to store the best parsimony scores and the number of
                // parsimony-optimal placements for each of the new samples to
                // be placed on the tree
                std::vector<int> best_parsimony_scores;
                std::vector<size_t> num_best_placements;

                for (size_t s=0; s<missing_samples.size(); s++) {
                    auto dfs = T->depth_first_expansion();
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

                    // Stores the parsimony score to place the sample at each
                    // node of the tree in DFS order.
                    std::vector<int> node_set_difference;
                    if (print_parsimony_scores) {
                        node_set_difference.resize(total_nodes);
                    }

                    size_t best_node_num_leaves = 0;
                    int best_set_difference = 1e9;
                    size_t best_j = 0;
                    size_t num_best = 1;
                    bool best_node_has_unique = false;
                    MAT::Node* best_node = T->root;

                    std::vector<bool> node_has_unique(total_nodes, false);
                    std::vector<size_t> best_j_vec;
                    best_j_vec.emplace_back(0);

                    // Parallel for loop to search for most parsimonious
                    // placements. Real action happens within mapper2_body
                    auto grain_size = 400; 
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes, grain_size),
                            [&](tbb::blocked_range<size_t> r) {
                            for (size_t k=r.begin(); k<r.end(); ++k){
                                mapper2_input inp;
                                inp.T = T;
                                inp.node = dfs[k];
                                inp.missing_sample_mutations = &missing_sample_mutations[s];
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
                    }); 

                    best_parsimony_scores.emplace_back(best_set_difference);
                    num_best_placements.emplace_back(num_best);
                }

                // Sort samples order in indexes based on parsimony scores
                // and number of parsimony-optimal placements
                if (sort_before_placement_1) {
                    std::stable_sort(indexes.begin(), indexes.end(),
                            [&num_best_placements, &best_parsimony_scores](size_t i1, size_t i2) 
                            {return ((best_parsimony_scores[i1] < best_parsimony_scores[i2]) || \
                                    ((best_parsimony_scores[i1] == best_parsimony_scores[i2]) && (num_best_placements[i1] < num_best_placements[i2])));});
                }
                else if (sort_before_placement_2) {
                    std::stable_sort(indexes.begin(), indexes.end(),
                            [&num_best_placements, &best_parsimony_scores](size_t i1, size_t i2) 
                            {return ((num_best_placements[i1] < num_best_placements[i2]) || \
                                    ((num_best_placements[i1] == num_best_placements[i2]) && (best_parsimony_scores[i1] < best_parsimony_scores[i2])));});
                }
                
                // Reverse sorted order if specified
                if (reverse_sort) {
                    std::reverse(indexes.begin(), indexes.end());
                }

                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }

            fprintf(stderr, "Adding missing samples to the tree.\n");  
        }

        // Traverse in sorted sample order
        for (size_t idx=0; idx<indexes.size(); idx++) {

            num_trees = optimal_trees.size();

            for (size_t t_idx=0; t_idx < num_trees; t_idx++) {
                timer.Start();
                
                T = &optimal_trees[t_idx];

                if (num_trees > 1) {
                    fprintf(stderr, "==Tree %zu=== \n", t_idx+1);
                }

                size_t s = indexes[idx];
                auto sample = missing_samples[s];

                if (print_parsimony_scores) {
                    auto parsimony_scores_filename = outdir + "/parsimony-scores.tsv";
                    if (s==0) {
                        fprintf(stderr, "\nNow computing branch parsimony scores for adding the missing samples at each of the %zu nodes in the existing tree without modifying the tree.\n", T->depth_first_expansion().size()); 
                        fprintf(stderr, "The branch parsimony scores will be written to file %s\n\n", parsimony_scores_filename.c_str());

                        parsimony_scores_file = fopen(parsimony_scores_filename.c_str(), "w");
                        fprintf (parsimony_scores_file, "#Sample\tTree node\tParsimony score\tOptimal (y/n)\tParsimony-increasing mutations (for optimal nodes)\n"); 
                    }
                }

                auto dfs = T->depth_first_expansion();
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

                if (print_parsimony_scores) {
                    node_set_difference.resize(total_nodes);
                }

                size_t best_node_num_leaves = 0;
                int best_set_difference = 1e9;
                
                size_t best_j = 0;
                bool best_node_has_unique = false;
                
                std::vector<bool> node_has_unique(total_nodes, false);
                std::vector<size_t> best_j_vec;
                
                size_t num_best = 1;
                MAT::Node* best_node = T->root;
                best_j_vec.emplace_back(0);

                // Parallel for loop to search for most parsimonious
                // placements. Real action happens within mapper2_body
                auto grain_size = 400; 
                tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes, grain_size),
                        [&](tbb::blocked_range<size_t> r) {
                        for (size_t k=r.begin(); k<r.end(); ++k){
                        mapper2_input inp;
                        inp.T = T;
                        inp.node = dfs[k];
                        inp.missing_sample_mutations = &missing_sample_mutations[s];
                        inp.excess_mutations = &node_excess_mutations[k];
                        inp.imputed_mutations = &node_imputed_mutations[k];
                        inp.best_node_num_leaves = &best_node_num_leaves;
                        inp.best_set_difference = &best_set_difference;
                        inp.best_node = &best_node;
                        inp.best_j =  &best_j;
                        inp.num_best = &num_best;
                        inp.j = k;
                        inp.has_unique = &best_node_has_unique;

                        if (print_parsimony_scores) {
                        inp.set_difference = &node_set_difference[k];
                        }
                        inp.best_j_vec = &best_j_vec;
                        inp.node_has_unique = &(node_has_unique);

                        mapper2_body(inp, print_parsimony_scores);
                        }       
                        }); 

                if (!print_parsimony_scores) {
                    fprintf(stderr, "Current tree size (#nodes): %zu\tSample name: %s\tParsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", total_nodes, sample.c_str(), \
                            best_set_difference, num_best);
                    // Prints a warning message if 2 or more
                    // parsimony-optimal placements found
                    if (num_best > 1) {
                        if (max_trees == 1) {
                            low_confidence_samples.emplace_back(sample);
                        }
                        if (num_best > max_uncertainty) {
                            fprintf(stderr, "WARNING: Number of parsimony-optimal placements exceeds maximum allowed value (%u). Ignoring sample %s.\n", max_uncertainty, sample.c_str());
                        }
                        else {
                            fprintf(stderr, "WARNING: Too many parsimony-optimal placements found. Placement done without high confidence.\n");
                        }
                    }
                }
                else {
                    fprintf(stderr, "Missing sample: %s\t Best parsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", sample.c_str(), \
                            best_set_difference, num_best);

                }

                // Debugging information to be printed if -DDEBUG compile-time
                // flag is set. This includes sample mutations, details of the
                // best node and the list of mutations at the best node
#if DEBUG == 1
                fprintf (stderr, "Sample mutations:\t");
                if (missing_sample_mutations[s].size() > 0) {
                    for (auto m: missing_sample_mutations[s]) {
                        if (m.is_missing) {
                            continue;
                        }
                        fprintf(stderr, "|%s", (MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc)).c_str());
                        fprintf(stderr, "| ");
                    }
                }
                fprintf (stderr, "\n");

                assert(num_best > 0);

                //best_node_vec.emplace_back(best_node);
                if ((num_best > 0) && (num_best <= max_uncertainty)) {
                    for (auto j: best_j_vec) {
                        auto node = dfs[j];

                        std::vector<std::string> muts;

                        fprintf(stderr, "Best node ");
                        if (node->is_leaf() || node_has_unique[j]) {
                            fprintf(stderr, "(sibling)");
                        }
                        else {
                            fprintf(stderr, "(child)");
                        }

                        if (node == best_node) {
                            fprintf(stderr, "*: %s\t", node->identifier.c_str());
                        }
                        else {
                            fprintf(stderr, ": %s\t", node->identifier.c_str());
                        }

                        std::string s = "|";
                        for (auto m: node->mutations) {
                            s += MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + '|';
                        }
                        if (node->mutations.size() > 0) {
                            muts.emplace_back(std::move(s));
                        }

                        for (auto anc: T->rsearch(node->identifier)) {
                            s = "|";
                            for (auto m: anc->mutations) {
                                s += MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + '|';
                            }
                            if (anc->mutations.size() > 0) {
                                muts.emplace_back(std::move(s));
                            }
                        }


                        std::reverse(muts.begin(), muts.end());

                        fprintf(stderr, "Mutations: "); 
                        for (size_t m = 0; m < muts.size(); m++) {
                            fprintf(stderr, "%s", muts[m].c_str());
                            if (m+1 < muts.size()) {
                                fprintf(stderr, " > "); 
                            }
                        }
                        fprintf(stderr, "\n"); 
                    }
                    fprintf(stderr, "\n"); 
                }

#endif
                
                // If number of parsimony-optimal trees is more than 1 and if
                // the number of trees has not already exceeded the maximum
                // limit, create a copy of the current tree in curr_tree
                MAT::Tree curr_tree;
                if ((max_trees > 1) && (num_best > 1) && (num_trees < max_trees)) {
                    curr_tree = MAT::get_tree_copy(*T);
                }

                if (print_parsimony_scores) {
                    for (size_t k = 0; k < total_nodes; k++) {
                        char is_optimal = (node_set_difference[k] == best_set_difference) ? 'y' : 'n';
                        fprintf (parsimony_scores_file, "%s\t%s\t%d\t\t%c\t", sample.c_str(), dfs[k]->identifier.c_str(), node_set_difference[k], is_optimal); 
                        if (node_set_difference[k] == best_set_difference) {
                            if (node_set_difference[k] == 0) {
                                fprintf(parsimony_scores_file, "*");
                            }
                            for (size_t idx = 0; idx < static_cast<size_t>(node_set_difference[k]); idx++) {
                                auto m = node_excess_mutations[k][idx];
                                assert ((m.mut_nuc & (m.mut_nuc-1)) == 0);
                                fprintf(parsimony_scores_file, "%s", (MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc)).c_str());
                                if (idx+1 < static_cast<size_t>(node_set_difference[k])) {
                                    fprintf(parsimony_scores_file, ",");
                                }
                            }
                        }
                        else {
                            fprintf(parsimony_scores_file, "N/A");
                        }
                        fprintf(parsimony_scores_file, "\n");
                    }
                }
                // Do placement only if number of parsimony-optimal placements
                // does not exceed the maximum allowed value
                else if (num_best <= max_uncertainty) {
                    if (num_best > 1) {
                        if (max_trees > 1) {
                            // Sorting by dfs order ensures reproducible results
                            // during multiple placements
                            std::sort(best_j_vec.begin(), best_j_vec.end());
                        }

                        // Update num_best so that the number of trees does
                        // not exceed maximum limit
                        if ((optimal_trees.size() <= max_trees) && (num_best + optimal_trees.size() > max_trees)) {
                            if ((num_best + optimal_trees.size() > max_trees+1) && (max_trees > 1))
                                fprintf (stderr, "%zu parsimony-optimal placements found but total trees has already exceed the max possible value (%i)!\n", num_best, max_trees);
                            num_best = 1 + max_trees - optimal_trees.size();
                        }
                    }

                    // Iterate over the number of parsimony-optimal placements
                    // for which a new tree will be created
                    for (size_t k = 0; k < num_best; k++) {

                        // best_j is updated using best_j_vec if multiple
                        // placements are allowed and the number of new trees
                        // for the given sample is greater than 1. If not, the
                        // default tie-breaking strategy used in mapper2_body has
                        // already chosen a single best_j
                        if ((max_trees > 1) && (num_best > 1)) {
                            if ((k==0) && (num_best > 1)) {
                                fprintf (stderr, "Creating %zu additional tree(s) for %zu parsimony-optimal placements.\n", num_best-1, num_best);
                            }
                            // If at second placement or higher, a new tree needs to
                            // be added to optimal_trees and T needs to point to its
                            // last element. If not, T is already pointing to the
                            // last element of optimal_trees on which placement will
                            // be carried out
                            if (k > 0) {
                                auto tmp_T = MAT::get_tree_copy(curr_tree);
                                optimal_trees.emplace_back(std::move(tmp_T));
                                T = &optimal_trees[optimal_trees.size()-1];
                                dfs = T->depth_first_expansion();
                            }

                            best_j = best_j_vec[k];                                                                                                                                                             
                            best_node_has_unique = node_has_unique[k];
                            best_node = dfs[best_j];
                        }
                        
                        // Ensure sample not already in the tree
                        if (T->get_node(sample) == NULL) {
                            // Is placement as sibling
                            if (best_node->is_leaf() || best_node_has_unique) {
                                std::string nid = std::to_string(++T->curr_internal_node);
                                T->create_node(nid, best_node->parent->identifier);
                                T->create_node(sample, nid);
                                T->move_node(best_node->identifier, nid);
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
                                    T->get_node(nid)->add_mutation(m);
                                }
                                // Add mutations to best node using l1_mut
                                for (auto m: l1_mut) {
                                    T->get_node(best_node->identifier)->add_mutation(m);
                                }
                                // Add new sample mutations using l2_mut
                                for (auto m: l2_mut) {
                                    T->get_node(sample)->add_mutation(m);
                                }
                            }
                            // Else placement as child
                            else {
                                T->create_node(sample, best_node->identifier);
                                MAT::Node* node = T->get_node(sample);
                                std::vector<MAT::Mutation> node_mut;

                                std::vector<MAT::Mutation> curr_l1_mut;

                                for (auto m1: best_node->mutations) {
                                    MAT::Mutation m = m1.copy();
                                    curr_l1_mut.emplace_back(m);
                                }

                                for (auto m1: node_excess_mutations[best_j]) {
                                    bool found = false;
                                    for (auto m2: curr_l1_mut) {
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

                            if (node_imputed_mutations[best_j].size() > 0) {
                                fprintf (stderr, "Imputed mutations:\t");
                                size_t tot = node_imputed_mutations[best_j].size();
                                for (size_t curr = 0; curr < tot; curr++) {
                                    if (curr < tot-1) {
                                        fprintf (stderr, "%i:%c;", node_imputed_mutations[best_j][curr].position, MAT::get_nuc(node_imputed_mutations[best_j][curr].mut_nuc));
                                    }
                                    else {
                                        fprintf (stderr, "%i:%c", node_imputed_mutations[best_j][curr].position, MAT::get_nuc(node_imputed_mutations[best_j][curr].mut_nuc));
                                    }
                                }
                                fprintf(stderr, "\n");
                            }
                        }
                        
                        if (max_trees == 1) {
                            break;
                        }
                    }
                }

                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }
        }
        
    }

    num_trees = optimal_trees.size();
            
    // If user specified print_parsimony_scores, close corresponding file and
    // terminate normally
    if (print_parsimony_scores) {
        if (parsimony_scores_file) {
            fclose(parsimony_scores_file);
        }
        return 0;
    }


    // Write final tree(s) to file(s)
    for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
        timer.Start();
            
        T = &optimal_trees[t_idx];
        
        auto final_tree_filename = outdir + "/final-tree.nh";
        if (num_trees > 1) {
            final_tree_filename = outdir + "/final-tree-" + std::to_string(t_idx+1) + ".nh";
            fprintf(stderr, "Writing final tree %zu to file %s \n", t_idx+1, final_tree_filename.c_str());
        }
        else {
            fprintf(stderr, "Writing final tree to file %s \n", final_tree_filename.c_str());
        }
        auto parsimony_score = T->get_parsimony_score();
        fprintf(stderr, "The parsimony score for this tree is: %zu \n", parsimony_score);
        FILE* final_tree_file = fopen(final_tree_filename.c_str(), "w");
        fprintf(final_tree_file, "%s\n", MAT::get_newick_string(*T, true, true, retain_original_branch_len).c_str());
        fclose(final_tree_file);

        tree_parsimony_scores.emplace_back(parsimony_score);
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    // If user need uncondensed tree output, write uncondensed tree(s) to
    // file(s)
    if (print_uncondensed_tree) {
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();
            
            T = &optimal_trees[t_idx];

            auto uncondensed_final_tree_filename = outdir + "/uncondensed-final-tree.nh";
            if (num_trees > 1) {
                uncondensed_final_tree_filename = outdir + "/uncondensed-final-tree-" + std::to_string(t_idx+1) + ".nh"; 
                fprintf(stderr, "Writing uncondensed final tree %zu to file %s \n", (t_idx+1), uncondensed_final_tree_filename.c_str());
            }
            else {
                fprintf(stderr, "Writing uncondensed final tree to file %s \n", uncondensed_final_tree_filename.c_str());
            }

            FILE* uncondensed_final_tree_file = fopen(uncondensed_final_tree_filename.c_str(), "w");

            if (T->condensed_nodes.size() > 0) {
                MAT::Tree T_to_print = MAT::get_tree_copy(*T);
                T_to_print.uncondense_leaves();
                fprintf(uncondensed_final_tree_file, "%s\n", MAT::get_newick_string(T_to_print, true, true, retain_original_branch_len).c_str());
            }
            else {
                fprintf(uncondensed_final_tree_file, "%s\n", MAT::get_newick_string(*T, true, true, retain_original_branch_len).c_str());
            }
            fclose(uncondensed_final_tree_file);

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    }

    if (missing_samples.size() > 0) {
        // For each final tree write the path of mutations from tree root to the
        // sample for each newly placed sample
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();

            T = &optimal_trees[t_idx];

            auto mutation_paths_filename = outdir + "/mutation-paths.txt";
            if (num_trees > 1) {
                mutation_paths_filename = outdir + "/mutation-paths-" + std::to_string(t_idx+1) + ".txt";
                fprintf(stderr, "Writing mutation paths for tree %zu to file %s \n", t_idx+1, mutation_paths_filename.c_str());
            }
            else {
                fprintf(stderr, "Writing mutation paths to file %s \n", mutation_paths_filename.c_str());
            }

            FILE* mutation_paths_file = fopen(mutation_paths_filename.c_str(), "w");

            for (size_t s=0; s<missing_samples.size(); s++) {
                auto sample = missing_samples[s];
                auto sample_node = T->get_node(sample);
                
                // If the missing sample is not found in the tree, it was not placed
                // because of max_uncertainty.  
                if (T->get_node(sample) == NULL) {
                    continue; 
                }

                // Stack for last-in first-out ordering
                std::stack<std::string> mutation_stack;
                std::string curr_node_mutation_string;

                // Mutations on the added sample
                auto curr_node_mutations = sample_node->mutations;
                if (curr_node_mutations.size() > 0) {
                    curr_node_mutation_string = sample + ":";
                    size_t num_mutations = curr_node_mutations.size();
                    for (size_t k = 0; k < num_mutations; k++) {
                        curr_node_mutation_string += MAT::get_nuc(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + MAT::get_nuc(curr_node_mutations[k].mut_nuc); 
                        if (k < num_mutations-1) {
                            curr_node_mutation_string += ',';
                        }
                        else {
                            curr_node_mutation_string += ' ';    
                        }
                    }
                    mutation_stack.push(curr_node_mutation_string);
                }

                // Mutations on the ancestors of added sample
                for (auto anc_node: T->rsearch(sample)) {
                    curr_node_mutations = anc_node->mutations;
                    if (curr_node_mutations.size() > 0) {
                        curr_node_mutation_string = anc_node->identifier + ":";
                        size_t num_mutations = curr_node_mutations.size();
                        for (size_t k = 0; k < num_mutations; k++) {
                            curr_node_mutation_string += MAT::get_nuc(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + MAT::get_nuc(curr_node_mutations[k].mut_nuc); 
                            if (k < num_mutations-1) {
                                curr_node_mutation_string += ',';
                            }
                            else {
                                curr_node_mutation_string += ' ';    
                            }
                        }
                        mutation_stack.push(curr_node_mutation_string);
                    }
                }

                fprintf(mutation_paths_file, "%s\t", sample.c_str()); 
                while (mutation_stack.size()) {
                    fprintf(mutation_paths_file, "%s", mutation_stack.top().c_str()); 
                    mutation_stack.pop();
                }
                fprintf(mutation_paths_file, "\n"); 
            }

            fclose(mutation_paths_file);

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    }

    if ((print_subtrees_size > 1) && (missing_samples.size() > 0)) {
        fprintf(stderr, "Computing subtrees for added samples. \n\n");

        // For each final tree, write a subtree of user-specified size around
        // each newly placed sample in newick format
        for (size_t t_idx = 0; t_idx < num_trees; t_idx++) {
            timer.Start();
            
            if (num_trees > 1) {
                fprintf(stderr, "==Tree %zu=== \n", t_idx+1);
            }

            T = &optimal_trees[t_idx];

            // Bool vector to mark which newly placed samples have already been
            // displayed in a subtree (initialized to false)
            std::vector<bool> displayed_mising_sample (missing_samples.size(), false);

            // If the missing sample is not found in the tree, it was not placed
            // because of max_uncertainty. Mark those samples as already
            // displayed. 
            for (size_t ms_idx = 0; ms_idx < missing_samples.size(); ms_idx++) {
                if (T->get_node(missing_samples[ms_idx]) == NULL) {
                    displayed_mising_sample[ms_idx] = true;
                }
            }

            int num_subtrees = 0;
            for (size_t i = 0; i < missing_samples.size(); i++) {
                if (displayed_mising_sample[i]) {
                    continue;
                }

                // Keep moving up the tree till a subtree of required size is
                // found
                for (auto anc: T->rsearch(missing_samples[i])) {
                    size_t num_leaves = T->get_num_leaves(anc);
                    if (num_leaves < print_subtrees_size) {
                        continue;
                    }

                    // if subtree of required size or more is found, create a
                    // copy of the subtree in new_T. Some tips may have to be
                    // pruned to meet the specified size constraints 
                    std::string newick = MAT::get_newick_string(*T, anc, false, true, retain_original_branch_len);
                    MAT::Tree new_T = MAT::create_tree_from_newick_string(newick);

                    // Map to store mutations for individual tree nodes in new_T
                    std::unordered_map<MAT::Node*, std::vector<MAT::Mutation>> subtree_node_mutations;
                    std::vector<MAT::Mutation> subtree_root_mutations;
                    auto dfs1 = T->depth_first_expansion(anc);
                    auto dfs2 = new_T.depth_first_expansion();

                    assert(dfs1.size() == dfs2.size());

                    for (size_t k = 0; k < dfs1.size(); k++) {
                        auto n1 = dfs1[k];
                        auto n2 = dfs2[k];
                        subtree_node_mutations[n2] = n1->mutations;

                        if (k == 0) {
                            for (auto p: T->rsearch(n1->identifier)) {
                                for (auto m: p->mutations) {
                                    subtree_root_mutations.emplace_back(m);
                                }
                            }
                            std::reverse(subtree_root_mutations.begin(), subtree_root_mutations.end());
                        }
                    }

                    // Prune tips if number of tips (leaves) currently is above
                    // the specified size in the subtree
                    if (num_leaves > print_subtrees_size) {
                        auto last_anc = new_T.get_node(missing_samples[i]);
                        auto ancestors = new_T.rsearch(missing_samples[i]);
                        if (ancestors.size() > 1) {
                            last_anc = ancestors[ancestors.size()-2];
                        }
                        std::vector<MAT::Node*> siblings;
                        for (auto child: new_T.root->children) {
                            if (child->identifier != last_anc->identifier) {
                                siblings.emplace_back(child);
                            }
                        }

                        for (size_t k=0; k<siblings.size(); k++) {
                            auto curr_sibling = siblings[k];
                            auto sibling_leaves = new_T.get_leaves(curr_sibling->identifier);
                            if (num_leaves-sibling_leaves.size() < print_subtrees_size) {
                                for (auto child: curr_sibling->children) {
                                    siblings.emplace_back(child);
                                }
                            }
                            else {
                                new_T.remove_node(curr_sibling->identifier, true);
                                num_leaves -= sibling_leaves.size();
                                if (num_leaves == print_subtrees_size) {
                                    break;
                                }
                            }
                        }

                        newick = MAT::get_newick_string(new_T, true, true, retain_original_branch_len);
                    }

                    tbb::parallel_for (tbb::blocked_range<size_t>(i+1, missing_samples.size(), 100),
                            [&](tbb::blocked_range<size_t> r) {
                            for (size_t j=r.begin(); j<r.end(); ++j){
                              for (size_t j = i+1; j < missing_samples.size(); j++) {
                                if (!displayed_mising_sample[j]) {
                                  if (new_T.get_node(missing_samples[j]) != NULL) {
                                    displayed_mising_sample[j] = true;
                                   }
                                 }
                               }
                             }
                            });

                    // Write subtree to file
                    ++num_subtrees;
                    auto subtree_filename = outdir + "/subtree-" + std::to_string(num_subtrees) + ".nh";
                    if (num_trees > 1) {
                        subtree_filename = outdir + "/" + "tree-" + std::to_string(t_idx+1) + "-subtree-" + std::to_string(num_subtrees) + ".nh";
                    }
                    fprintf(stderr, "Writing subtree %d to file %s.\n", num_subtrees, subtree_filename.c_str());
                    FILE* subtree_file = fopen(subtree_filename.c_str(), "w");
                    fprintf(subtree_file, "%s\n", newick.c_str());
                    fclose(subtree_file);

                    // Write list of mutations on the subtree to file
                    auto subtree_mutations_filename = outdir + "/subtree-" + std::to_string(num_subtrees) + "-mutations.txt";
                    if (num_trees > 1) {
                        subtree_mutations_filename = outdir + "/" + "tree-" + std::to_string(t_idx+1) + "-subtree-" + std::to_string(num_subtrees) + "-mutations.txt";;
                    }
                    fprintf(stderr, "Writing list of mutations at the nodes of subtree %d to file %s\n", num_subtrees, subtree_mutations_filename.c_str());
                    FILE* subtree_mutations_file = fopen(subtree_mutations_filename.c_str(), "w");

                    // Mutations on root
                    fprintf(subtree_mutations_file, "ROOT->%s: ", new_T.root->identifier.c_str());
                    size_t tot_mutations = subtree_root_mutations.size();
                    for (size_t idx = 0; idx < tot_mutations; idx++) {
                        auto m = subtree_root_mutations[idx];
                        fprintf(subtree_mutations_file, "%s", (MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc)).c_str());
                        if (idx+1 <tot_mutations) {
                            fprintf(subtree_mutations_file, ",");
                        }
                    }
                    fprintf(subtree_mutations_file, "\n");

                    // Mutations on remaining nodes
                    for (auto n: new_T.depth_first_expansion()) {
                        size_t tot_mutations = subtree_node_mutations[n].size();
                        fprintf(subtree_mutations_file, "%s: ", n->identifier.c_str());
                        for (size_t idx = 0; idx < tot_mutations; idx++) {
                            auto m = subtree_node_mutations[n][idx];
                            fprintf(subtree_mutations_file, "%s", (MAT::get_nuc(m.par_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc)).c_str());
                            if (idx+1 <tot_mutations) {
                                fprintf(subtree_mutations_file, ",");
                            }
                        }
                        fprintf(subtree_mutations_file, "\n");
                    }

                    fclose(subtree_mutations_file);

                    // Expand internal nodes that are condensed
                    bool has_condensed = false;
                    FILE* subtree_expanded_file = NULL;
                    for (auto l: new_T.get_leaves()) {
                        if (T->condensed_nodes.find(l->identifier) != T->condensed_nodes.end()) {
                            if (!has_condensed) {
                                auto subtree_expanded_filename = outdir + "/subtree-" + std::to_string(num_subtrees) + "-expanded.txt";
                                if (num_trees > 1) {
                                    subtree_expanded_filename = outdir + "/tree-" + std::to_string(1+t_idx) + "-subtree-" + std::to_string(num_subtrees) + "-expanded.txt";
                                }
                                fprintf(stderr, "Subtree %d has condensed nodes.\nExpanding the condensed nodes for subtree %d in file %s\n", num_subtrees, num_subtrees, subtree_expanded_filename.c_str());
                                subtree_expanded_file = fopen(subtree_expanded_filename.c_str(), "w");
                                has_condensed = true;
                            }
                            fprintf(subtree_expanded_file, "%s: ", l->identifier.c_str());
                            for (auto n: T->condensed_nodes[l->identifier]) {
                                fprintf(subtree_expanded_file, "%s ", n.c_str());
                            }
                            fprintf(subtree_expanded_file, "\n");
                        }
                    }
                    if (has_condensed) {
                        fclose(subtree_expanded_file);
                    }
                    break;
                }
            }
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
    }

    // Print warning message with a list of all samples placed with low
    // confidence (>=4 parsimony-optimal placements
    if (low_confidence_samples.size() > 0) {
        fprintf(stderr, "WARNING: Following samples had multiple possibilities of parsimony-optimal placements:\n");
        for (auto lcs: low_confidence_samples) { 
            fprintf(stderr, "%s\n", lcs.c_str());
        }
    }

    // Store mutation-annotated tree to a protobuf file if user has asked for it
    if (dout_filename != "") {

        timer.Start();

        fprintf(stderr, "Saving mutation-annotated tree object to file %s\n", dout_filename.c_str());
        if (num_trees > 1) {
            fprintf(stderr, "WARNING: --multiple-placements option was used but only the first mutation-annotated tree object will be saved to file.\n");
        }

        Parsimony::data data;

        T = &optimal_trees[0];
        MAT::save_mutation_annotated_tree(*T, dout_filename);
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    
    /*
    // When multiple placements was used, print the tree names with lowest total
    // parsimony scores
    if (max_trees > 1) {
        size_t min_parsimony_score = *std::min_element(tree_parsimony_scores.begin(), tree_parsimony_scores.end());
        if (optimal_trees.size() > 1) {
            fprintf(stderr, "Trees below have the lowest parsimony score of %zu among all trees: \n", min_parsimony_score);
            for (size_t t_idx = 0; t_idx < optimal_trees.size(); t_idx++) {
                if (tree_parsimony_scores[t_idx] == min_parsimony_score) {
                    fprintf(stderr, "tree-%zu ", t_idx+1);
                }
            }
            fprintf(stderr, "\n\n");
        }
        else {
            fprintf(stderr, "Single best tree with a parsimony score of %zu was found during multiple placements.\n\n", min_parsimony_score);;
        }
    }
    */

    google::protobuf::ShutdownProtobufLibrary();

#if SAVE_PROFILE == 1
    Instrumentor::Get().EndSession();
#endif

    return 0;
}

