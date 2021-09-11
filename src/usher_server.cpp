#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <memory>
#include <limits>
#include <filesystem>
#include <utility>
#include "boost/filesystem.hpp"
#include "usher_graph.hpp"
#include "parsimony.pb.h"
#include "version.hpp"
#include "usher_common.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
namespace fs = boost::filesystem;




int main(int argc, char** argv) {

    //Variables to load command-line options using Boost program_options
    std::string arg_dirname;
    std::string MAT_list_filename;
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    po::options_description desc{"Options"};
    uint32_t sleep_length;
    uint32_t termination_character;


    desc.add_options()

    ("arguments,a", po::value<std::string>(&arg_dirname)->required(), "Input argument directory that will contain argument files with arguments for usher [REQUIRED]")
    ("list-mutation-annonated-trees,i", po::value<std::string>(&MAT_list_filename)->default_value(""), "File containing list of mutation-annotated tree objects")
    ("sleep-length,s", po::value<uint32_t>(&sleep_length)->default_value(100), "Time in milliseconds that the program waits until checking for input in argument file")
    ("termination-char,c", po::value<uint32_t>(&termination_character)->default_value(94), "Character that will ddetermine if the argument file is ready to be read. Default is '^'")
    ("help,h", "Print help messages");


    po::options_description all_options;
    all_options.add(desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        //return with error code 1 unless the user specifies help
        std::cerr << desc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    //slot for a MAT not in the MAT_list file 
    MAT::Tree loaded_MAT;
    bool loaded_MAT_avail = false;//keep track if loaded_MAT is a new copy of MAT and can be used
    std::string loaded_MAT_name = "";
    
    fs::path p = fs::current_path();
    p/=arg_dirname;//get the path to the argument directory

    if (!fs::is_directory(p)) {
        fprintf(stderr, "ERROR: Argument directory provided is not a directory: %s!\n", arg_dirname.c_str());
        exit(1);
    }
    MAT::Tree *curr_tree; //MAT that is used in the iteration
    Timer timer;
    std::unordered_map<std::string, MAT::Tree> MAT_list; //store list of trees
    std::unordered_map<std::string, bool> MAT_list_avail; //stores information on whether the MATs in the list are available for use

    //if there's a specified MAT_list file, it would go through the file and load the MATs in
    if(MAT_list_filename != "") {
        if(!fs::exists(MAT_list_filename)) {
            std::cout << "MAT list file not found" <<std::endl;
            return 1;
        }
        std::string MAT_filename;
        MAT::Tree temp_MAT;
        std::ifstream MAT_list_file(MAT_list_filename);
        while(std::getline(MAT_list_file, MAT_filename)) { //set up MATs in the list
            // Load mutation-annotated tree and store it
            timer.Start();
            fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", MAT_filename.c_str());
            MAT_list[MAT_filename] = MAT::load_mutation_annotated_tree(MAT_filename);
            MAT_list_avail[MAT_filename] = true;
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
        MAT_list_file.close();
    }



    while(true) {
        //if there's a MAT loaded in and it's not available to be used, load it in
        if((loaded_MAT_name != "") && (!loaded_MAT_avail)) {
            timer.Start();
            fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", loaded_MAT_name.c_str());
            // Load mutation-annotated tree and store it
            MAT::clear_tree(loaded_MAT);
            loaded_MAT = MAT::load_mutation_annotated_tree(loaded_MAT_name);
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            loaded_MAT_avail = true;
        }
        //iterate through MATs specified by MAT_list file if there are used trees then remove them and load them back in
        for(auto itr = MAT_list_avail.begin(); itr != MAT_list_avail.end(); itr++) {
            if(!(itr->second)) { //if a MAT pointed by this iterator is not available, load
                timer.Start();
                fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", (itr->first).c_str());
                MAT::clear_tree(MAT_list[itr->first]);
                MAT_list[itr->first] = MAT::load_mutation_annotated_tree(itr->first);
                itr->second = true;
                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            }
        }

        //if there are no files in the argument directory, wait.
        if(fs::is_empty(p)) {
            fprintf(stderr, "Waiting for more arguments\n\n");
            while(fs::is_empty(p)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(sleep_length));
            }
        }

        //sort files in the directory by when they were last modified

        std::vector<std::pair<fs::path, time_t>> list_argfiles;
        
        
        for(auto ent: fs::directory_iterator(p)){
            fs::path curr_path = ent.path();
            list_argfiles.emplace_back(std::make_pair(curr_path, fs::last_write_time(curr_path)));
        }
        std::sort(list_argfiles.begin(), list_argfiles.end());
        
        //go through each file in the argument directory
        for(auto arg_file_pair: list_argfiles){
            std::ifstream arguments_file(arg_file_pair.first.string());
            if (!arguments_file) {
                fprintf(stderr, "ERROR: Could not open the arguments file: %s!\n", arg_file_pair.first.string().c_str());
                exit(1);
            }

            arguments_file.seekg(-2, arguments_file.end);
            char last_char = arguments_file.get();//last char could be at the very end or one before
            if(last_char != ((char) termination_character)) {
                last_char = arguments_file.get();
                if(last_char != ((char) termination_character)) {//if the termination character does not exist, skip the file
                    continue;
                }
            }
            
            arguments_file.seekg(0, arguments_file.beg);
            std::string argument;

            //get a line of argument and feed it into usher
            while(std::getline(arguments_file, argument)){
                argument.erase(std::remove(argument.begin(), argument.end(), ((char) termination_character)), argument.end());
                fprintf(stderr, "Argument: %s \n\n", argument.c_str());//print this line's argument
                std::istringstream arg(argument);
                std::vector<std::string> arg_vector; //to store each word from arg
                arg_vector.emplace_back("./usher"); //to replicate commandline argument
                std::string tempStr; //hold the word from argument before adding to the vector
                while(arg >> tempStr) {
                    arg_vector.emplace_back(tempStr);
                }

                int argc_line = arg_vector.size();
                const char* argv_line[argc_line];
                for(int i = 0; i < argc_line; i++) {
                    argv_line[i] = arg_vector[i].c_str();
                }

                //Variables to load command-line options using Boost program_options
                std::string din_filename;
                std::string dout_filename;
                std::string outdir;
                std::string vcf_filename;
                uint32_t num_threads;
                uint32_t max_trees = 1; //only one tree for usher_server
                uint32_t max_uncertainty;
                uint32_t max_parsimony;
                bool sort_before_placement_1 = false;
                bool sort_before_placement_2 = false;
                bool sort_before_placement_3 = false;
                bool reverse_sort = false;
                bool collapse_tree=false;
                bool collapse_output_tree=false;
                bool print_uncondensed_tree = false;
                bool print_parsimony_scores = false;
                bool retain_original_branch_len = false;
                bool no_add = false;
                bool detailed_clades = false;
                size_t print_subtrees_size=0;
                size_t print_subtrees_single=0;
                po::options_description desc{"Options"};

                std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
                desc.add_options()
                ("vcf,v", po::value<std::string>(&vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
                ("load-mutation-annotated-tree,i", po::value<std::string>(&din_filename)->required(), "Load mutation-annotated tree object [REQUIRED]")
                ("outdir,d", po::value<std::string>(&outdir)->default_value("."), "Output directory to dump output and log files [DEFAULT uses current directory]")
                ("save-mutation-annotated-tree,o", po::value<std::string>(&dout_filename)->default_value(""), "Save output mutation-annotated tree object to the specified filename")
                ("sort-before-placement-1,s", po::bool_switch(&sort_before_placement_1), \
                "Sort new samples based on computed parsimony score and then number of optimal placements before the actual placement [EXPERIMENTAL].")
                ("sort-before-placement-2,S", po::bool_switch(&sort_before_placement_2), \
                "Sort new samples based on the number of optimal placements and then the parsimony score before the actual placement [EXPERIMENTAL].")
                ("sort-before-placement-3,A", po::bool_switch(&sort_before_placement_3), \
                "Sort new samples based on the number of ambiguous bases [EXPERIMENTAL].")
                ("reverse-sort,r", po::bool_switch(&reverse_sort), \
                "Reverse the sorting order of sorting options (sort-before-placement-1 or sort-before-placement-2) [EXPERIMENTAL]")
                ("collapse-tree,c", po::bool_switch(&collapse_tree), \
                "Collapse internal nodes of the input tree with no mutations and condense identical sequences in polytomies into a single node and the save the tree to file condensed-tree.nh in outdir")
                ("collapse-output-tree,C", po::bool_switch(&collapse_output_tree), \
                "Collapse internal nodes of the output tree with no mutations before the saving the tree to file final-tree.nh in outdir")
                ("max-uncertainty-per-sample,e", po::value<uint32_t>(&max_uncertainty)->default_value(1e6), \
                "Maximum number of equally parsimonious placements allowed per sample beyond which the sample is ignored")
                ("max-parsimony-per-sample,E", po::value<uint32_t>(&max_parsimony)->default_value(1e6), \
                "Maximum parsimony score of the most parsimonious placement(s) allowed per sample beyond which the sample is ignored")
                ("write-uncondensed-final-tree,u", po::bool_switch(&print_uncondensed_tree), "Write the final tree in uncondensed format and save to file uncondensed-final-tree.nh in outdir")
                ("write-subtrees-size,k", po::value<size_t>(&print_subtrees_size)->default_value(0), \
                "Write minimum set of subtrees covering the newly added samples of size equal to this value")
                ("write-single-subtree,K", po::value<size_t>(&print_subtrees_single)->default_value(0), \
                "Similar to write-subtrees-size but produces a single subtree with all newly added samples along with random samples up to the value specified by this argument")
                ("write-parsimony-scores-per-node,p", po::bool_switch(&print_parsimony_scores), \
                "Write the parsimony scores for adding new samples at each existing node in the tree without modifying the tree in a file names parsimony-scores.tsv in outdir")
                ("retain-input-branch-lengths,l", po::bool_switch(&retain_original_branch_len), \
                "Retain the branch lengths from the input tree in out newick files instead of using number of mutations for the branch lengths.")
                ("no-add,n", po::bool_switch(&no_add), \
                "Do not add new samples to the tree")
                ("detailed-clades,D", po::bool_switch(&detailed_clades), \
                "In clades.txt, write a histogram of annotated clades and counts across all equally parsimonious placements")
                ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
                ("version", "Print version number")
                ("reload", "Reload the MAT_list")
                ("help,h", "Print help messages");

                po::options_description all_options;
                all_options.add(desc);

                po::variables_map vm;
                try {
                    po::store(po::command_line_parser(argc_line, argv_line).options(all_options).run(), vm);
                    po::notify(vm);
                } catch(std::exception &e) {
                    if (vm.count("version")) {
                        std::cout << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
                    } else if(vm.count("reload")) { //reload option
                        if(MAT_list_filename != "") {
                            if(!boost::filesystem::exists(MAT_list_filename)) {
                                std::cout << "MAT list file not found" <<std::endl;
                                return 1;
                            }
                            //if the trees specified in the MAT_list file is loaded in already, 
                            //discard them so new version could be loaded in
                            for(auto itr = MAT_list.begin(); itr != MAT_list.end(); itr++) {
                                MAT::clear_tree(itr->second);//delete all the trees in the list
                            }
                            MAT_list.clear();
                            MAT_list_avail.clear();
                            std::string MAT_filename;
                            MAT::Tree temp_MAT;
                            std::ifstream MAT_list_file(MAT_list_filename);
                            while(std::getline(MAT_list_file, MAT_filename)) { //set up MATs in the list
                                // Load mutation-annotated tree and store it
                                timer.Start();
                                fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", MAT_filename.c_str());
                                MAT_list[MAT_filename] = MAT::load_mutation_annotated_tree(MAT_filename);
                                MAT_list_avail[MAT_filename] = true;
                                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
                            }
                            MAT_list_file.close();
                        }
                    } else {
                        std::cerr << "UShER (v" << PROJECT_VERSION << ")" << std::endl;
                        std::cerr << desc << std::endl;
                    }
                    // Return with error code 1 unless the user specifies help or version
                    if(vm.count("help") || vm.count("version") || vm.count("reload"))
                        continue;//if help or version then go to next line
                    else
                        break;//if error encountered then stop reading this file
                }
                
                if(MAT_list.count(din_filename) != 0) { //if the MAT is in the list
                    if(!(MAT_list_avail[din_filename])) { //if the MAT is not available, load
                        timer.Start();
                        fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", din_filename.c_str());
                        MAT::clear_tree(MAT_list[din_filename]);
                        MAT_list[din_filename] = MAT::load_mutation_annotated_tree(din_filename);
                        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
                    }
                    curr_tree = &MAT_list[din_filename];
                    MAT_list_avail[din_filename] = false;

                } else if(din_filename != loaded_MAT_name) {//if the MAT is not in the loaded_MAT slot

                    timer.Start();
                    fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", din_filename.c_str());

                    if(loaded_MAT_name != "") { //if there is an existing trees, delete it
                        MAT::clear_tree(loaded_MAT);
                    }
                    // Load mutation-annotated tree and store it
                    loaded_MAT = MAT::load_mutation_annotated_tree(din_filename);
                    loaded_MAT_name = din_filename;
                    curr_tree = &loaded_MAT;
                    loaded_MAT_avail = false;
                    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

                //if the tree is in the loaded MAT slot, but not available,
                //load the tree that can be used
                } else if(!loaded_MAT_avail) {
                    timer.Start();
                    fprintf(stderr, "Loading existing mutation-annotated tree object from file %s\n", din_filename.c_str());
                    MAT::clear_tree(loaded_MAT);
                    loaded_MAT = MAT::load_mutation_annotated_tree(din_filename);
                    curr_tree = &loaded_MAT;
                    loaded_MAT_avail = false;
                    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
                } else { //loaded_MAT can be used
                    curr_tree = &loaded_MAT;
                    loaded_MAT_avail = false;
                }



                // Variables below used to store the different fields of the input VCF file
                bool header_found = false;
                std::vector<std::string> variant_ids;
                std::vector<Missing_Sample> missing_samples;

                // Vector used to store all tree nodes in breadth-first search (BFS) order
                std::vector<MAT::Node*> bfs;
                // Map the node identifier string to index in the BFS traversal
                std::unordered_map<std::string, size_t> bfs_idx;

                // Vectore to store the names of samples which have a high number of
                // parsimony-optimal placements
                std::vector<std::string> low_confidence_samples;
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
                } catch(const boost::iostreams::gzip_error& e) {
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
                                if ((curr_tree->get_node(words[j]) == NULL) && (curr_tree->condensed_leaves.find(words[j]) == curr_tree->condensed_leaves.end())) {
                                    missing_samples.emplace_back(Missing_Sample(words[j]));
                                    missing_idx.emplace_back(j);
                                } else {
                                    fprintf(stderr, "WARNING: Ignoring sample %s as it is already in the tree.\n", words[j].c_str());
                                }
                            }
                            header_found = true;
                        }
                    } else if (header_found) {
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
                                        } else {
                                            auto nuc = MAT::get_nuc_id(allele[0]);
                                            if (nuc == MAT::get_nuc_id('N')) {
                                                m.is_missing = true;
                                            } else {
                                                m.is_missing = false;
                                            }
                                            m.mut_nuc = nuc;
                                        }
                                        (*iter).mutations.emplace_back(m);
                                    }
                                } else {
                                    m.is_missing = true;
                                    m.mut_nuc = MAT::get_nuc_id('N');
                                    (*iter).mutations.emplace_back(m);
                                }
                                if ((m.mut_nuc & (m.mut_nuc-1)) !=0) {
                                    (*iter).num_ambiguous++;
                                }
                            }
                        }
                    }
                }
                fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

                //run usher on the argument
                int return_val = usher_common(dout_filename, outdir, num_threads, max_trees, max_uncertainty, max_parsimony,
                                            sort_before_placement_1, sort_before_placement_2, sort_before_placement_3, reverse_sort, collapse_tree,
                                            collapse_output_tree, print_uncondensed_tree, print_parsimony_scores, retain_original_branch_len, no_add,
                                            detailed_clades, print_subtrees_size, print_subtrees_single, missing_samples, low_confidence_samples, curr_tree);

                if(return_val != 0) {
                    break;//if error encountered then stop reading the file for now
                }
            }
	    fs::remove(arg_file_pair.first);//arguments were run, so delete the file
        }
    }
}
