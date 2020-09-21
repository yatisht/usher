#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/tbb.h>
#include "usher_graph.hpp"
#include "parsimony.pb.h"

namespace po = boost::program_options;

std::vector<int8_t> get_nuc_id (char c) {
    switch (c) {
        case 'a':
        case 'A': return std::vector<int8_t>{0};
        case 'c':
        case 'C': return std::vector<int8_t>{1};
        case 'g':
        case 'G': return std::vector<int8_t>{2};
        case 't':
        case 'T': return std::vector<int8_t>{3};
        case 'R': return std::vector<int8_t>{0,2};
        case 'Y': return std::vector<int8_t>{1,3};
        case 'S': return std::vector<int8_t>{1,2};
        case 'W': return std::vector<int8_t>{0,3};
        case 'K': return std::vector<int8_t>{2,3};
        case 'M': return std::vector<int8_t>{0,1};
        case 'B': return std::vector<int8_t>{1,2,3};
        case 'D': return std::vector<int8_t>{0,2,3};
        case 'H': return std::vector<int8_t>{0,1,3};
        case 'V': return std::vector<int8_t>{0,1,2};
        case 'n':
        case 'N': return std::vector<int8_t>{0,1,2,3};
        default: return std::vector<int8_t>{0,1,2,3};
    }
}

char get_nuc_char (int8_t nuc_id) {
    switch (nuc_id) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default : return 'N'; 
    }
}

int main(int argc, char** argv) {

    std::string tree_filename;
    std::string din_filename;
    std::string dout_filename;
    std::string vcf_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    bool collapse_tree=false;
    bool print_uncondensed_tree = false;
    bool print_parsimony_scores = false;
    size_t print_subtrees_size=0;
    po::options_description desc{"Options"};

    std::string num_threads_message = "Number of threads (default uses all available cores, " + std::to_string(num_cores) + " detected on this machine)";
    desc.add_options()
        ("vcf", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed .gz format)")
        ("tree", po::value<std::string>(&tree_filename)->default_value(""), "Input tree file")
        ("load-assignments", po::value<std::string>(&din_filename)->default_value(""), "Load mutation-annotated tree object")
        ("save-assignments", po::value<std::string>(&dout_filename)->default_value(""), "Save output mutation-annotated tree object")
        ("threads", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("collapse-final-tree", po::bool_switch(&collapse_tree), "Collapse internal nodes of the output tree with no mutations and condense identical sequences in polytomies into a single node")
        ("print-uncondensed-final-tree", po::bool_switch(&print_uncondensed_tree), "Print the final tree in uncondensed format")
        ("print-subtrees-size", po::value<size_t>(&print_subtrees_size)->default_value(0), \
         "Print minimum set of subtrees covering the newly added samples of size equal to or larger than this value")
        ("print-parsimony-scores-per-node", po::bool_switch(&print_parsimony_scores), \
         "Print the parsimony scores for adding new samples at each existing node in the tree without modifying the tree")
        ("help", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);

    //    po::positional_options_description p;
    //    p.add("tree", 1);
    //    p.add("vcf", 1);
    //    p.add("threads", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    if (print_subtrees_size == 1) {
        std::cerr << "ERROR: print-subtrees-size should be larger than 1\n";
        return 1;
    }

    if ((din_filename != "") && collapse_tree) {
        std::cerr << "ERROR: cannot load assignments and collapse tree simulaneously.\n";
        return 1;
    }

    if (print_parsimony_scores) {
        if (collapse_tree || print_uncondensed_tree || (print_subtrees_size > 0) || (dout_filename != "")) {
            fprintf (stderr, "WARNING: --print-parsimony-scores-per-node is set. Will terminate without modifying the original tree.\n");
        }
    }

    fprintf(stderr, "Initializing %u worker threads.\n\n", num_threads);

    Timer timer; 

    omp_lock_t omplock;
    omp_set_num_threads(num_threads);
    omp_init_lock(&omplock);
        
    tbb::task_scheduler_init init(num_threads);

#if SAVE_PROFILE == 1
    Instrumentor::Get().BeginSession("test-main", "p1.json");
#endif

    Tree T;

    bool header_found = false;
    std::vector<std::string> variant_ids;
    std::vector<std::string> missing_samples;
    std::unordered_map<Node*, std::vector<mutation>> node_mutations;
    std::vector<std::vector<mutation>> missing_sample_mutations;
    size_t num_missing = 0;

    std::unordered_map<std::string, std::vector<std::string>> condensed_nodes;
    Tree condensed_T; 
    std::unordered_map<Node*, std::vector<mutation>> condensed_node_mutations;
    
    std::vector<Node*> bfs;
    std::unordered_map<std::string, size_t> bfs_idx;
    
    std::vector<std::string> low_confidence_samples;

    if (tree_filename != "") {
        T = TreeLib::create_tree_from_newick(tree_filename);
        bfs = T.breadth_first_expansion();

        for (size_t idx = 0; idx < bfs.size(); idx++) {
            bfs_idx[bfs[idx]->identifier] = idx;
        }

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

        tbb::flow::graph mapper_graph;

        tbb::flow::function_node<mapper_input, int> mapper(mapper_graph, tbb::flow::unlimited, mapper_body());
        tbb::flow::source_node <mapper_input> reader (mapper_graph,
                [&] (mapper_input &inp) -> bool {
                std::string s;
                std::getline(instream, s);
                std::vector<std::string> words;
                TreeLib::split(s, words);
                inp.variant_pos = -1;
                if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                for (size_t j=9; j < words.size(); j++) {
                variant_ids.emplace_back(words[j]);
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
                    TreeLib::split(words[4], ',', alleles);
                    inp.T = &T;
                    inp.chrom = words[0];
                    inp.bfs = &bfs;
                    inp.bfs_idx = &bfs_idx;
                    inp.variant_ids = &variant_ids;
                    inp.missing_samples = &missing_samples;
                    inp.node_mutations = &node_mutations;
                    inp.missing_sample_mutations = &missing_sample_mutations;
                    auto ref_nucs = get_nuc_id(words[3][0]);
                    assert(ref_nucs.size() == 1);
                    inp.ref_nuc = ref_nucs[0]; 
                    inp.variants.clear();
                    for (size_t j=9; j < words.size(); j++) {
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) { 
                                std::string allele = alleles[allele_id-1];
                                inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id(allele[0])));
                            }
                        }
                        else {
                            inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id('N')));
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
    else if (din_filename != "") {
        fprintf(stderr, "Loading existing assignments\n");
        timer.Start();
        
        Parsimony::data data;

        std::ifstream inpfile(din_filename, std::ios::in | std::ios::binary);
        if (!inpfile) {
            fprintf(stderr, "ERROR: Could not load the assignments file: %s!\n", din_filename.c_str());
            exit(1);
        }
        data.ParseFromIstream(&inpfile);
        inpfile.close();

        T = TreeLib::create_tree_from_newick_string(data.newick());

        auto dfs = T.depth_first_expansion();
        std::unordered_set<std::string> condensed_leaves;

        for (size_t idx = 0; idx < dfs.size(); idx++) {
            auto mutation_list = data.node_mutations(idx);
            auto node = dfs[idx];
            node_mutations.insert(std::pair<Node*, std::vector<mutation>>(node, std::vector<mutation>()));  
            for (int k = 0; k < mutation_list.mutation_size(); k++) {
                auto mut = mutation_list.mutation(k);
                mutation m;
                m.chrom = mut.chromosome();
                m.position = mut.position();
                m.ref_nuc = mut.ref_nuc();
                m.par_nuc = mut.par_nuc();
                for (int n = 0; n < mut.mut_nuc_size(); n++) {
                    m.mut_nuc.emplace_back(mut.mut_nuc(n));
                }
                node_mutations[node].emplace_back(m);
            }
        }
        
        size_t num_condensed_nodes = static_cast<size_t>(data.condensed_nodes_size());
        for (size_t idx = 0; idx < num_condensed_nodes; idx++) {
            auto cn = data.condensed_nodes(idx);
            condensed_nodes.insert(std::pair<std::string, std::vector<std::string>>(cn.node_name(), std::vector<std::string>(cn.condensed_leaves_size())));
            for (int k = 0; k < cn.condensed_leaves_size(); k++) {
                condensed_nodes[cn.node_name()][k] =cn.condensed_leaves(k);
                condensed_leaves.insert(cn.condensed_leaves(k));
            }
        }
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        
        fprintf(stderr, "Loading VCF file\n");
        timer.Start();

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
        while (instream.peek() != EOF) {
            std::getline(instream, s);
            std::vector<std::string> words;
            TreeLib::split(s, words);
            if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                    for (size_t j=9; j < words.size(); j++) {
                        variant_ids.emplace_back(words[j]);
                        if ((T.get_node(words[j]) == NULL) && (condensed_leaves.find(words[j]) == condensed_leaves.end())) {
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
                TreeLib::split(words[4], ',', alleles);
                for (size_t k = 0; k < missing_idx.size(); k++) {
                    size_t j = missing_idx[k];
                    auto iter = missing_samples.begin();
                    std::advance(iter, k);
                    if (iter != missing_samples.end()) {
                        auto mutations_iter = missing_sample_mutations.begin() + (iter - missing_samples.begin());
                        mutation m;
                        m.position = std::stoi(words[1]);
                        auto ref_nucs = get_nuc_id(words[3][0]);
                        assert(ref_nucs.size() == 1);
                        m.ref_nuc = ref_nucs[0];
                        m.par_nuc = ref_nucs[0];
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) { 
                                std::string allele = alleles[allele_id-1];
                                if (allele[0] == 'N') {
                                    m.is_missing = true;
                                }
                                else {
                                    auto nucs = get_nuc_id(allele[0]);
                                    if (nucs.size() == 4) {
                                        m.is_missing = true;
                                    }
                                    else {
                                        m.is_missing = false;
                                        for (auto n: nucs) {
                                            m.mut_nuc.emplace_back(n);
                                        }
                                    }
                                }
                                (*mutations_iter).emplace_back(m);
                            }
                        }
                        else {
                            m.is_missing = true;
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

    fprintf(stderr, "Found %zu missing samples.\n", missing_samples.size()); 

    // Timer timer;

    if (missing_samples.size() > 0) {
        fprintf(stderr, "Sorting node mutations by positions.\n");  
        timer.Start();
#pragma omp parallel for
        for (size_t k=0; k<node_mutations.size(); k++) {
            auto iter = node_mutations.begin();
            std::advance(iter, k);
            if (!std::is_sorted(iter->second.begin(), iter->second.end(), compare_by_position)) {
                std::sort(iter->second.begin(), iter->second.end(), compare_by_position); 
            }
        }
#pragma omp parallel for
        for (size_t k=0; k<missing_samples.size(); k++) {
            if (!std::is_sorted(missing_sample_mutations[k].begin(), missing_sample_mutations[k].end(), compare_by_position)) {
                std::sort(missing_sample_mutations[k].begin(), missing_sample_mutations[k].end(), compare_by_position); 
            }
        }

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        if (print_parsimony_scores) {
            fprintf(stderr, "Printing current tree (with internal nodes labelled). \n");
            fprintf(stdout, "%s\n", TreeLib::get_newick_string(T, true, true).c_str());

            fprintf(stderr, "\nComputing parsimony score for adding the missing samples at each of the %zu nodes in the existing tree without modifying the tree.\n\n", T.depth_first_expansion().size()); 
        } 
        else {
            fprintf(stderr, "Adding missing samples to the tree.\n");  
        }

        for (size_t s=0; s<missing_samples.size(); s++) {
            timer.Start();
            
            auto sample = missing_samples[s];

            auto dfs = T.depth_first_expansion();
            size_t total_nodes = dfs.size();
            
            std::vector<std::vector<mutation>> node_excess_mutations(total_nodes);
            std::vector<std::vector<mutation>> node_imputed_mutations(total_nodes);

            std::vector<int> node_set_difference;

            if (print_parsimony_scores) {
                node_set_difference.resize(total_nodes);
            }

            size_t best_node_num_leaves = 0;
            int best_set_difference = 1e9;
            size_t best_j = 0;
            size_t num_best = 1;
            bool best_node_has_unique = false;
            Node* best_node = T.root;
#if DEBUG == 1
            std::vector<bool> node_has_unique(total_nodes, false);
            std::vector<size_t> best_j_vec;
            best_j_vec.emplace_back(0);
#endif

            auto grain_size = 400; 
            tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes, grain_size),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k){
                        mapper2_input inp;
                        inp.T = &T;
                        inp.node = dfs[k];
                        inp.node_mutations = &node_mutations;
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
#if DEBUG == 1
                        inp.best_j_vec = &best_j_vec;
                        inp.node_has_unique = &(node_has_unique);
#endif

                        mapper2_body(inp, print_parsimony_scores);
                    }       
            }); 

            if (!print_parsimony_scores) {
                fprintf(stderr, "Current tree size (#nodes): %zu\tMissing sample: %s\tParsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", total_nodes, sample.c_str(), \
                        best_set_difference, num_best);
                if (num_best > 3) {
                    low_confidence_samples.push_back(sample);
                    fprintf(stderr, "WARNING: Too many parsimony-optimal placements found. Placement done without high confidence.\n");
                }
            }
            else {
                fprintf(stderr, "Missing sample: %s\t Best parsimony score: %d\tNumber of parsimony-optimal placements: %zu\n", sample.c_str(), \
                        best_set_difference, num_best);
            }

#if DEBUG == 1
            
            fprintf (stderr, "Sample mutations:\t");
            if (missing_sample_mutations[s].size() > 0) {
                for (auto m: missing_sample_mutations[s]) {
                    if (m.is_missing) {
                        continue;
                    }
                    fprintf(stderr, "|%s", (get_nuc_char(m.par_nuc) + std::to_string(m.position)).c_str());
                    for (size_t c_size =0; c_size < m.mut_nuc.size(); c_size++) {
                        fprintf(stderr, "%c", get_nuc_char(m.mut_nuc[c_size]));
                        if (c_size + 1 < m.mut_nuc.size()) {
                            fprintf(stderr, ",");
                        }
                    }
                    fprintf(stderr, "| ");
                }
            }
            fprintf (stderr, "\n");

            assert((best_j_vec.size() == num_best) && (num_best > 0));

            //best_node_vec.emplace_back(best_node);
            if (num_best > 0) {
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
                    for (auto m: node_mutations[node]) {
                        s += get_nuc_char(m.par_nuc) + std::to_string(m.position) + get_nuc_char(m.mut_nuc[0]) + '|';
                    }
                    if (node_mutations[node].size() > 0) {
                        muts.push_back(std::move(s));
                    }
                    
                    for (auto anc: T.rsearch(node->identifier)) {
                        s = "|";
                        for (auto m: node_mutations[anc]) {
                            s += get_nuc_char(m.par_nuc) + std::to_string(m.position) + get_nuc_char(m.mut_nuc[0]) + '|';
                        }
                        if (node_mutations[anc].size() > 0) {
                            muts.push_back(std::move(s));
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
            
            assert(std::find(best_j_vec.begin(), best_j_vec.end(), best_j) != best_j_vec.end());
#endif

            if (print_parsimony_scores) {
                fprintf (stdout, "#Sample\tTree node\tParsimony score\tParsimony-increasing mutations (for optimal nodes)\n"); 
                for (size_t k = 0; k < total_nodes; k++) {
                    fprintf (stdout, "%s\t%s\t%d\t", sample.c_str(), dfs[k]->identifier.c_str(), node_set_difference[k]); 
                    if (node_set_difference[k] == best_set_difference) {
                        assert (node_excess_mutations[k].size() == node_set_difference[k]);
                        if (node_set_difference[k] == 0) {
                            fprintf(stdout, "*");
                        }
                        for (auto m: node_excess_mutations[k]) {
                            fprintf(stdout, "|%s", (get_nuc_char(m.par_nuc) + std::to_string(m.position)).c_str());
                            for (size_t c_size =0; c_size < m.mut_nuc.size(); c_size++) {
                                fprintf(stdout, "%c", get_nuc_char(m.mut_nuc[c_size]));
                                if (c_size + 1 < m.mut_nuc.size()) {
                                    fprintf(stdout, ",");
                                }
                            }
                            fprintf(stdout, "| ");
                        }
                    }
                    fprintf(stdout, "\n");
                }
            }
            else {
                if (T.get_node(sample) == NULL) {
                    if (best_node->is_leaf() || best_node_has_unique) {
                        std::string nid = std::to_string(++T.curr_internal_node);
                        T.create_node(nid, best_node->parent->identifier);
                        T.create_node(sample, nid);
                        T.move_node(best_node->identifier, nid);
                        std::vector<mutation> common_mut, l1_mut, l2_mut;
                        std::vector<mutation> curr_l1_mut;

                        if (node_mutations.find(best_node) != node_mutations.end()) {
                            for (auto m1: node_mutations[best_node]) {
                                mutation m;
                                m.position = m1.position;
                                m.ref_nuc = m1.ref_nuc;
                                m.par_nuc = m1.par_nuc;
                                m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                curr_l1_mut.emplace_back(m);
                            }
                            node_mutations.erase(best_node);
                        }

                        for (auto m1: curr_l1_mut) {
                            bool found = false;
                            for (auto m2: node_excess_mutations[best_j]) {
                                if (m1.position == m2.position) {
                                    if (m1.mut_nuc[0] == m2.mut_nuc[0]) {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (!found) {
                                mutation m;
                                m.position = m1.position;
                                m.ref_nuc = m1.ref_nuc;
                                m.par_nuc = m1.par_nuc;
                                m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                l1_mut.emplace_back(m);
                            }
                        }
                        for (auto m1: node_excess_mutations[best_j]) {
                            bool found = false;
                            for (auto m2: curr_l1_mut) {
                                if (m1.position == m2.position) {
                                    if (m1.mut_nuc[0] == m2.mut_nuc[0]) {
                                        found = true;
                                        mutation m;
                                        m.position = m1.position;
                                        m.ref_nuc = m1.ref_nuc;
                                        m.par_nuc = m1.par_nuc;
                                        m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                        common_mut.emplace_back(m);
                                        break;
                                    }
                                }
                            }
                            if (!found) {
                                mutation m;
                                m.position = m1.position;
                                m.ref_nuc = m1.ref_nuc;
                                m.par_nuc = m1.par_nuc;
                                m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                l2_mut.emplace_back(m);
                            }
                        }

                        if (common_mut.size() > 0) {
                            if (!std::is_sorted(common_mut.begin(), common_mut.end(), compare_by_position)) {
                                std::sort(common_mut.begin(), common_mut.end(), compare_by_position);
                            }
                            node_mutations[T.get_node(nid)] = common_mut;
                        }
                        if (l1_mut.size() > 0) {
                            if (!std::is_sorted(l1_mut.begin(), l1_mut.end(),compare_by_position)) {
                                std::sort(l1_mut.begin(), l1_mut.end(),compare_by_position); 
                            }
                            node_mutations[T.get_node(best_node->identifier)] = l1_mut;
                        }
                        if (l2_mut.size() > 0) {
                            if (!std::is_sorted(l2_mut.begin(), l2_mut.end(), compare_by_position)) {
                                std::sort(l2_mut.begin(), l2_mut.end(), compare_by_position); 
                            }
                            node_mutations[T.get_node(sample)] = l2_mut;
                        }
                    }
                    else {
                        T.create_node(sample, best_node->identifier);
                        Node* node = T.get_node(sample);
                        std::vector<mutation> node_mut;

                        std::vector<mutation> curr_l1_mut;

                        if (node_mutations.find(best_node) != node_mutations.end()) {
                            for (auto m1: node_mutations[best_node]) {
                                mutation m;
                                m.position = m1.position;
                                m.ref_nuc = m1.ref_nuc;
                                m.par_nuc = m1.par_nuc;
                                m.mut_nuc.emplace_back(m1.mut_nuc[0]);
                                curr_l1_mut.emplace_back(m);
                            }
                        }

                        for (auto m1: node_excess_mutations[best_j]) {
                            bool found = false;
                            for (auto m2: curr_l1_mut) {
                                if (m1.position == m2.position) {
                                    if (m1.mut_nuc[0] == m2.mut_nuc[0]) {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (!found) {
                                mutation m;
                                m.position = m1.position;
                                m.ref_nuc = m1.ref_nuc;
                                m.par_nuc = m1.par_nuc;
                                for (auto nuc: m1.mut_nuc) {
                                    m.mut_nuc.emplace_back(nuc);
                                }
                                node_mut.emplace_back(m);
                            }
                        }
                        std::sort(node_mut.begin(), node_mut.end(), compare_by_position); 
                        node_mutations[node] = node_mut;
                    }

                    if (node_imputed_mutations[best_j].size() > 0) {
                        fprintf (stderr, "Imputed mutations:\t");
                        size_t tot = node_imputed_mutations[best_j].size();
                        for (size_t curr = 0; curr < tot; curr++) {
                            if (curr < tot-1) {
                                fprintf (stderr, "%i:%c;", node_imputed_mutations[best_j][curr].position, get_nuc_char(node_imputed_mutations[best_j][curr].mut_nuc[0]));
                            }
                            else {
                                fprintf (stderr, "%i:%c", node_imputed_mutations[best_j][curr].position, get_nuc_char(node_imputed_mutations[best_j][curr].mut_nuc[0]));
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                }
            }

            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
        
    }
            
    if (print_parsimony_scores) {
        return 0;
    }

    if (collapse_tree) {
        fprintf(stderr, "Collapsing final tree. \n");
        timer.Start();

        bfs.clear();
        bfs = T.breadth_first_expansion();
        
        for (size_t idx = 1; idx < bfs.size(); idx++) {
            auto mutations = node_mutations[bfs[idx]];
            if (mutations.size() == 0) {
                auto node = bfs[idx];
                auto parent = node->parent;
                auto children = node->children;
                for (auto child: children) {
                    T.move_node(child->identifier, parent->identifier);
                }
            }
        }
        
        condensed_T = TreeLib::create_tree_from_newick_string(TreeLib::get_newick_string(T, false, true));
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        fprintf(stderr, "Condensing identical sequences. \n");
        timer.Start();

        bfs.clear();
        bfs = T.breadth_first_expansion();
        auto condensed_bfs = condensed_T.breadth_first_expansion();

        assert(condensed_bfs.size() == bfs.size());

        condensed_nodes.clear();

        for (size_t it = 0; it < condensed_bfs.size(); it++) {
            auto condensed_node = condensed_bfs[it];
            condensed_node_mutations.insert(std::pair<Node*, std::vector<mutation>>(condensed_node, std::vector<mutation>(node_mutations[bfs[it]].size())));
            for (size_t k = 0; k < node_mutations[bfs[it]].size(); k++) {
                condensed_node_mutations[condensed_node][k] = node_mutations[bfs[it]][k];
            }
        }

        auto tree_leaves = T.get_leaves();
        for (auto l1: tree_leaves) {
            std::vector<std::string> polytomy_nodes;

            if (std::find(missing_samples.begin(), missing_samples.end(), l1->identifier) != missing_samples.end()) {
                continue;
            }
            if (node_mutations[l1].size() > 0) {
                continue;
            }
            if (condensed_T.get_node(l1->identifier) == NULL) {
                continue;
            }

            for (auto l2: l1->parent->children) {
                if (std::find(missing_samples.begin(), missing_samples.end(), l2->identifier) != missing_samples.end()) {
                    continue;
                }
                if (l2->is_leaf() && (condensed_T.get_node(l2->identifier) != NULL) && (node_mutations[l2].size() == 0)) {
                    polytomy_nodes.push_back(l2->identifier);
                }
            }

            if (polytomy_nodes.size() > 1) {
                std::string new_node_name = "node_" + std::to_string(1+condensed_nodes.size()) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
                auto curr_node = condensed_T.get_node(l1->identifier);
                condensed_T.create_node(new_node_name, curr_node->parent->identifier, l1->branch_length);
                auto new_node = condensed_T.get_node(new_node_name);
                condensed_node_mutations.insert(std::pair<Node*, std::vector<mutation>>(new_node, std::vector<mutation>(0)));
                condensed_nodes[new_node_name] = std::vector<std::string>(polytomy_nodes.size());

                for (size_t it = 0; it < polytomy_nodes.size(); it++) {
                    condensed_nodes[new_node_name][it] = polytomy_nodes[it];
                    condensed_T.remove_node(polytomy_nodes[it], false);
                }
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        
        fprintf(stderr, "Printing condensed tree. \n");
        fprintf(stdout, "%s\n", TreeLib::get_newick_string(condensed_T, true, true).c_str());
    }

    fprintf(stderr, "Printing final tree. \n");
    fprintf(stdout, "%s\n", TreeLib::get_newick_string(T, true, true).c_str());

    if (print_uncondensed_tree) {
        fprintf(stderr, "Printing uncondensed final tree. \n");
        
        timer.Start();
        
        if (!collapse_tree && (condensed_nodes.size() > 0)) {
            Tree T_to_print = TreeLib::create_tree_from_newick_string(TreeLib::get_newick_string(T, false, true)); 
            for (size_t it = 0; it < condensed_nodes.size(); it++) {
                auto cn = condensed_nodes.begin();
                std::advance(cn, it);

                auto n = T_to_print.get_node(cn->first);
                auto par = (n->parent != NULL) ? n->parent : n;

                size_t num_samples = cn->second.size();

                if (num_samples > 0) {
                    T_to_print.rename_node(n->identifier, cn->second[0]);
                }
                
                for (size_t s = 1; s < num_samples; s++) {
                    T_to_print.create_node(cn->second[s], par->identifier, n->branch_length);
                }
            }
            fprintf(stdout, "%s\n", TreeLib::get_newick_string(T_to_print, true, true).c_str());
        }
        else {
            fprintf(stdout, "%s\n", TreeLib::get_newick_string(T, true, true).c_str());
        }
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    
    if (missing_samples.size() > 0) {
        fprintf(stderr, "Printing mutation paths in the new tree for missing samples.\n");  

        timer.Start();
        
        for (size_t s=0; s<missing_samples.size(); s++) {
            auto sample = missing_samples[s];
            auto sample_node = T.get_node(sample);
            std::stack<std::string> mutation_stack;
            std::string curr_node_mutation_string;

            auto curr_node_mutations = node_mutations[sample_node];
            if (curr_node_mutations.size() > 0) {
                curr_node_mutation_string = sample + ":";
                size_t num_mutations = curr_node_mutations.size();
                for (size_t k = 0; k < num_mutations; k++) {
                    curr_node_mutation_string += get_nuc_char(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + get_nuc_char(curr_node_mutations[k].mut_nuc[0]); 
                    if (k < num_mutations-1) {
                        curr_node_mutation_string += ',';
                    }
                    else {
                        curr_node_mutation_string += ' ';    
                    }
                }
                mutation_stack.push(curr_node_mutation_string);
            }

            for (auto anc_node: T.rsearch(sample)) {
                curr_node_mutations = node_mutations[anc_node];
                if (curr_node_mutations.size() > 0) {
                    curr_node_mutation_string = anc_node->identifier + ":";
                    size_t num_mutations = curr_node_mutations.size();
                    for (size_t k = 0; k < num_mutations; k++) {
                        curr_node_mutation_string += get_nuc_char(curr_node_mutations[k].par_nuc) + std::to_string(curr_node_mutations[k].position) + get_nuc_char(curr_node_mutations[k].mut_nuc[0]); 
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
            
            fprintf(stderr, "%s\t", sample.c_str()); 
            while (mutation_stack.size()) {
                fprintf(stderr, "%s", mutation_stack.top().c_str()); 
                mutation_stack.pop();
            }
            fprintf(stderr, "\n"); 
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    if ((print_subtrees_size > 1) && (missing_samples.size() > 0)) {
        fprintf(stderr, "Printing subtrees for added samples. \n");

        timer.Start();
        
        std::vector<bool> displayed_mising_sample (missing_samples.size(), false);
        
        int num_subtrees = 0;
        for (size_t i = 0; i < missing_samples.size(); i++) {
            if (displayed_mising_sample[i]) {
                continue;
            }
            
            for (auto anc: T.rsearch(missing_samples[i])) {
                size_t num_leaves = T.get_num_leaves(anc);
                if (num_leaves < print_subtrees_size) {
                    continue;
                }

                std::string newick = TreeLib::get_newick_string(T, anc, false, true);
                Tree new_T = TreeLib::create_tree_from_newick_string(newick);

                std::unordered_map<Node*, std::vector<mutation>> subtree_node_mutations;
                std::vector<mutation> subtree_root_mutations;
                auto dfs1 = T.depth_first_expansion(anc);
                auto dfs2 = new_T.depth_first_expansion();

                assert(dfs1.size() == dfs2.size());

                for (size_t k = 0; k < dfs1.size(); k++) {
                    auto n1 = dfs1[k];
                    auto n2 = dfs2[k];
                    subtree_node_mutations[n2] = node_mutations[n1];

                    if (k == 0) {
                        for (auto p: T.rsearch(n1->identifier)) {
                            for (auto m: node_mutations[p]) {
                                subtree_root_mutations.push_back(m);
                            }
                        }
                        std::reverse(subtree_root_mutations.begin(), subtree_root_mutations.end());
                    }
                }

                if (num_leaves > print_subtrees_size) {
                    auto last_anc = new_T.get_node(missing_samples[i]);
                    auto ancestors = new_T.rsearch(missing_samples[i]);
                    if (ancestors.size() > 1) {
                        last_anc = ancestors[ancestors.size()-2];
                    }
                    std::vector<Node*> siblings;
                    for (auto child: new_T.root->children) {
                        if (child->identifier != last_anc->identifier) {
                            siblings.push_back(child);
                        }
                    }
                    
                    for (size_t k=0; k<siblings.size(); k++) {
                        auto curr_sibling = siblings[k];
                        auto sibling_leaves = new_T.get_leaves(curr_sibling->identifier);
                        if (num_leaves-sibling_leaves.size() < print_subtrees_size) {
                            for (auto child: curr_sibling->children) {
                                siblings.push_back(child);
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

                    newick = TreeLib::get_newick_string(new_T, true, true);
                }

#pragma omp parallel for
                for (size_t j = i+1; j < missing_samples.size(); j++) {
                    if (!displayed_mising_sample[j]) {
                        if (new_T.get_node(missing_samples[j]) != NULL) {
                            displayed_mising_sample[j] = true;
                        }
                    }
                }
                
                fprintf(stderr, "Subtree %d.\n", ++num_subtrees);
                fprintf(stdout, "%s\n", newick.c_str());

                fprintf(stderr, "Mutations at the nodes of subtree %d.\n", num_subtrees);
                if (subtree_root_mutations.size() > 0) {
                    fprintf(stderr, "ROOT->%s: ", new_T.root->identifier.c_str());
                    for (auto m: subtree_root_mutations) {
                        fprintf(stderr, "|%s", (get_nuc_char(m.par_nuc) + std::to_string(m.position)).c_str());
                        for (size_t c_size =0; c_size < m.mut_nuc.size(); c_size++) {
                            fprintf(stderr, "%c", get_nuc_char(m.mut_nuc[c_size]));
                            if (c_size + 1 < m.mut_nuc.size()) {
                                fprintf(stderr, ",");
                            }
                        }
                        fprintf(stderr, "| ");
                    }
                    fprintf(stderr, "\n");
                }
                for (auto n: new_T.depth_first_expansion()) {
                    fprintf(stderr, "%s: ", n->identifier.c_str());
                    for (auto m: subtree_node_mutations[n]) {
                        fprintf(stderr, "|%s", (get_nuc_char(m.par_nuc) + std::to_string(m.position)).c_str());
                        for (size_t c_size =0; c_size < m.mut_nuc.size(); c_size++) {
                            fprintf(stderr, "%c", get_nuc_char(m.mut_nuc[c_size]));
                            if (c_size + 1 < m.mut_nuc.size()) {
                                fprintf(stderr, ",");
                            }
                        }
                        fprintf(stderr, "| ");
                    }
                    fprintf(stderr, "\n");
                }

                bool has_condensed = false;
                for (auto l: new_T.get_leaves()) {
                    if (condensed_nodes.find(l->identifier) != condensed_nodes.end()) {
                        if (!has_condensed) {
                            fprintf(stderr, "\nExpanded condensed nodes for subtree %d:\n", num_subtrees);
                            has_condensed = true;
                        }
                        fprintf(stderr, "%s: ", l->identifier.c_str());
                        for (auto n: condensed_nodes[l->identifier]) {
                            fprintf(stderr, "%s ", n.c_str());
                        }
                        fprintf(stderr, "\n");
                    }
                }
                if (has_condensed) {
                    fprintf(stderr, "\n");
                }

                break;
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    if (low_confidence_samples.size() > 0) {
        fprintf(stderr, "WARNING: Following samples had several possibilities of parsimony-optimal placements:\n");
        for (auto lcs: low_confidence_samples) { 
            fprintf(stderr, "%s\n", lcs.c_str());
        }
    }

    if (dout_filename != "") {
        fprintf(stderr, "Saving assignments. \n");

        timer.Start();

        Parsimony::data data;

        if (!collapse_tree) {
            data.set_newick(TreeLib::get_newick_string(T, false, true));

            auto dfs = T.depth_first_expansion();

            for (size_t idx = 0; idx < dfs.size(); idx++) {
                auto mutation_list = data.add_node_mutations();
                auto mutations = node_mutations[dfs[idx]];
                for (auto m: mutations) {
                    auto mut = mutation_list->add_mutation();
                    mut->set_chromosome(m.chrom);
                    mut->set_position(m.position);
                    mut->set_ref_nuc(m.ref_nuc);
                    mut->set_par_nuc(m.par_nuc);
                    mut->clear_mut_nuc();
                    for (auto nuc: m.mut_nuc) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            }
        }
        else {
            data.set_newick(TreeLib::get_newick_string(condensed_T, false, true));

            auto dfs = condensed_T.depth_first_expansion();

            for (size_t idx = 0; idx < dfs.size(); idx++) {
                auto mutation_list = data.add_node_mutations();
                auto mutations = condensed_node_mutations[dfs[idx]];
                for (auto m: mutations) {
                    auto mut = mutation_list->add_mutation();
                    mut->set_chromosome(m.chrom);
                    mut->set_position(m.position);
                    mut->set_ref_nuc(m.ref_nuc);
                    mut->set_par_nuc(m.par_nuc);
                    mut->clear_mut_nuc();
                    for (auto nuc: m.mut_nuc) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            }

            // Add condensed nodes
            for (auto cn: condensed_nodes) {
                auto cn_ptr = data.add_condensed_nodes();
                cn_ptr->set_node_name(cn.first);
                for (auto l: cn.second) {
                    cn_ptr->add_condensed_leaves(l);
                }
            }
        }

        std::ofstream outfile(dout_filename, std::ios::out | std::ios::binary);
        data.SerializeToOstream(&outfile);
        outfile.close();
        
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    
    google::protobuf::ShutdownProtobufLibrary();

#if SAVE_PROFILE == 1
    Instrumentor::Get().EndSession();
#endif

    return 0;
}

