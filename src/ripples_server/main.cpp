#include <src/mutation_annotated_tree.hpp>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include "src/usher_graph.hpp"
#include "src/usher_common.hpp"
#include "src/ripples/ripples_fast/ripples_runner.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

int main(int argc, char** argv) {
    std::string tree_filename;
    std::string sample_filename;
    std::string vcf_filename;
    uint32_t num_threads;

    po::options_description desc("Options");
    desc.add_options()("input-mat,i", po::value<std::string>(&tree_filename)->required(),
                       "Input mutation-annotated tree file [REQUIRED].")
                       ("samples,s", po::value<std::string>(&sample_filename)->required(),
                       "Select samples by explicitly naming them. One per line [REQUIRED].")
                       ("read-vcf,v", po::value<std::string>(&vcf_filename)->default_value(""),
                       "input VCF file containing the samples to be palced")
                       ("threads,T", po::value<uint32_t>(&num_threads)->default_value(tbb::task_scheduler_init::default_num_threads()),
                       "Number of Threads to use")
                       ("help,h", "Print help messages");


    po::options_description all_options;
    all_options.add(desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help or version
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }
        
    fprintf(stderr, "\nInitializing %u worker threads.\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);

    Timer timer;
    MAT::Tree T;

    // Load mutation-annotated tree
    timer.Start();
    fprintf(stderr, "Loading existing mutation-annotated tree from file %s\n", tree_filename.c_str());
    T = MAT::load_mutation_annotated_tree(tree_filename);

    if (T.root == NULL) {
        fprintf(stderr, "ERROR: Empty tree.\n");
        exit(1);
    }

    T.uncondense_leaves();
    fprintf(stderr, "Tree loaded in %ld msec \n\n", timer.Stop());
    fprintf(stderr, "LEAVES: %ld\n\n", T.get_num_leaves());

    // Loading Sequence List
    timer.Start();
    std::vector<std::string> sample_names;
    std::ifstream infile(sample_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", sample_filename.c_str());
        exit(1);
    }
    
    fprintf(stderr, "Loading the sequence names from file %s\n", sample_filename.c_str());
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() > 1) {
            fprintf(stderr, "WARNING: Sequence file %s contains excess columns; ignoring\n", sample_filename.c_str());
        } else if (words.size() == 0) {
            fprintf(stderr, "WARNING: Empty line in Sequence file %s; ignoring\n", sample_filename.c_str());
            continue;
        }
        //remove carriage returns from the input to handle windows os
        auto sname = words[0];
        if (sname[sname.size()-1] == '\r') {
            sname = sname.substr(0,sname.size()-1);
        }
        sample_names.push_back(std::move(sname));
    }
    infile.close();
    fprintf(stderr, "Sequences loaded in %ld msec \n\n", timer.Stop());

    // Loading Sequences
    fprintf(stderr, "Loading the sequences from file %s\n", vcf_filename.c_str());
    std::vector<Missing_Sample> vcf_samples;
    if (vcf_filename != "")
    {
        timer.Start();
        boost::filesystem::ifstream fileHandler(vcf_filename);
        std::string s;
        bool header_found = false;
        while (getline(fileHandler, s)) {
            std::vector<std::string> words;
            MAT::string_split(s, words);
            if (words.size() > 1) {
                //Checking for header
                if (words[1] == "POS") 
                {
                    header_found = true;
                    //Leave certain fields based on our VCF format
                    for (int j = 9; j < (int)words.size(); j++)
                    {
                        vcf_samples.emplace_back(Missing_Sample(words[j]));
                    }
                }
                else if (header_found) 
                {
                    std::vector<std::string> alleles;
                    //Checking for different alleles at a site
                    MAT::string_split(words[4], ',', alleles);
                    for (int j = 9; j < (int)words.size(); j++) 
                    {
                        int idx = j - 9;
                        auto iter = vcf_samples.begin();
                        std::advance(iter, idx);

                        MAT::Mutation m;
                        m.chrom = words[0];
                        m.position = std::stoi(words[1]);
                        //Checking the mutating allele value within the allele sizes
                        if (std::stoi(words[j]) > int(alleles.size())) 
                        {
                            fprintf(stderr, "\n\nVCF ERROR at Position: %d, idx = %d, Allele_id: %d, Alleles_size: %ld\n\n", m.position, idx, std::stoi(words[j]), alleles.size());
                        }
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
                                (*iter).mutations.emplace_back(m);
                            }
                        } 
                        else {
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

        fprintf(stderr, "VCF read in %ld msec \n", timer.Stop());
    }
    else
    {
        fprintf(stderr, "No VCF file provided! \n");
    }

    // Finding missing samples
    std::vector<Missing_Sample> missing_samples;
    for (auto s: sample_names)
    {
        if (T.get_node(s) != NULL)
        {
            fprintf(stderr, "WARNING: %s is already present in the tree.\n", s.c_str());
        }
        else
        {
            auto it = std::find_if(vcf_samples.begin(), vcf_samples.end(),
                       [&](const Missing_Sample& ms) {
                           return ms.name == s;
                       });
            if (it != vcf_samples.end()) {
                missing_samples.emplace_back(*it);
            }
        }
    }

    // Adding missing samples to tree
    if (missing_samples.size())
    {
        // Copying the Tree
        timer.Start();
        MAT::Tree T_new = T.fast_tree_copy();
        fprintf(stderr, "\nTree copied in %ld msec \n\n", timer.Stop());
            
        std::vector<std::string> low_confidence_samples;
        int return_val = usher_common("", ".", 1, 1e6, 1e6,
                            false, false, false, false, false,
                            false, false, false, false, false,
                            false, 0, 0, missing_samples, low_confidence_samples, &T_new);

        if(return_val != 0) {
            exit(1);
        }
        
        // FIX for new nodes NOT present in all_nodes
				/*
        std::vector<MAT::Node*> new_nodes;
        for (auto node: T_new.depth_first_expansion()) {
            if (std::find(sample_names.begin(), sample_names.end(), node->identifier) != sample_names.end()) {
                new_nodes.push_back(node);
            }
        }
				*/
        T_new.update_all_nodes(T_new.depth_first_expansion());

				// Set ripples-fast parameters TODO: Can be set from cli options with defaults
				uint32_t branch_len = 3;
        std::string outdir = ".";
        uint32_t num_desc = 10;
        ripples_runner runner(T_new, missing_samples, num_threads, branch_len,
                              num_desc, outdir);
        runner();
    }
}
