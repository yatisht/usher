#include <src/mutation_annotated_tree.hpp>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include "src/usher_graph.hpp"

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
    std::unordered_map<std::string, std::vector<MAT::Mutation>> samples;
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
                        samples.insert(std::make_pair(words[j], std::vector<MAT::Mutation>{}));
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
                                samples[words[j]].emplace_back(m);
                            }
                        } 
                        else {
                            m.is_missing = true;
                            m.mut_nuc = MAT::get_nuc_id('N');
                            samples[words[j]].emplace_back(m);
                        }
                    }
                }
            }
        }

        fprintf(stderr, "VCF read in %ld msec \n\n", timer.Stop());
    }
    else
        fprintf(stderr, "No VCF file provided! \n\n");

    // Placing the missing sequences on the tree
    for (auto s: sample_names)
    {
        if (T.get_node(s) == NULL)
        {
            fprintf(stderr, "%s is missing from tree.\n", s.c_str());
            
            // Copying the Tree
            timer.Start();
            MAT::Tree T_new;
            
            std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
            std::vector<MAT::Node*> dfs_orig = T.breadth_first_expansion();
            std::vector<MAT::Node> dfs_new(dfs_orig.size());
            int idx = 0;

            remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(T.root, NULL));
            while (!remaining_nodes.empty())
            {
                MAT::Node* new_node = &dfs_new[idx++];
                auto orig_node = remaining_nodes.front().first;
                auto new_node_parent = remaining_nodes.front().second;
                remaining_nodes.pop();
                
                // Create node
                if (new_node_parent == NULL) {
                    new_node->identifier = orig_node->identifier;
                    new_node->parent = NULL;
                    new_node->level = 1;
                    new_node->branch_length = orig_node->branch_length;
                    new_node->mutations.clear();
                    new_node->clade_annotations.clear();
                    T_new.root = new_node;
                }
                else {
                    new_node->identifier = orig_node->identifier;
                    new_node->parent = new_node_parent;
                    new_node->level = new_node_parent->level + 1;
                    new_node->branch_length = orig_node->branch_length;
                    new_node->mutations.clear();
                    new_node->clade_annotations.clear();
                    new_node_parent->children.emplace_back(new_node);
                }

                //Add children to remaining_nodes    
                for (auto child: orig_node->children) {
                    remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, new_node));
                }
            }
            
            // Adding annotations and mutations
            static tbb::affinity_partitioner ap;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs_orig.size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t k=r.begin(); k<r.end(); ++k) {
                    auto orig_node = dfs_orig[k];
                    auto new_node = &dfs_new[k];
                    // Add clade annotations
                    new_node->clade_annotations.resize(orig_node->clade_annotations.size()); 
                    for (size_t i=0; i < orig_node->clade_annotations.size(); i++) {
                        new_node->clade_annotations[i] = orig_node->clade_annotations[i];
                    }
                    // Add mutations
                    new_node->mutations.resize(orig_node->mutations.size());
                    for (size_t i=0; i < orig_node->mutations.size(); i++) {
                        new_node->mutations[i] = orig_node->mutations[i].copy();
                    }
                }
            }, ap);

            fprintf(stderr, "Tree copied in %ld msec \n\n", timer.Stop());
            fprintf(stderr, "Root children: %ld\n\n", T_new.root->children.size());
            fprintf(stderr, "LEAVES: %ld\n\n", T_new.get_num_leaves());

            ////CHECK
            //for (int i = 0; i < dfs_orig.size(); i++)
            //{
            //    auto orig_node = dfs_orig[i];
            //    auto new_node = &dfs_new[i];

            //    if ((orig_node->level == new_node->level) && (orig_node->branch_length == new_node->branch_length) && (orig_node->identifier == new_node->identifier) && (orig_node->clade_annotations[0] == new_node->clade_annotations[0]) && (orig_node->clade_annotations[1] == new_node->clade_annotations[1]) && ((orig_node->parent == NULL) || ((orig_node->parent != NULL) && (orig_node->parent->identifier == new_node->parent->identifier))))
            //    {
            //        for (int i = 0; i < orig_node->children.size(); i++)
            //        {
            //            if (orig_node->children[i]->identifier != new_node->children[i]->identifier)
            //                fprintf(stderr, "%s children are not equal", orig_node->identifier.c_str());
            //        }
            //        for (int i = 0; i < orig_node->mutations.size(); i++)
            //        {
            //            if (orig_node->mutations[i].mut_nuc != new_node->mutations[i].mut_nuc)
            //                fprintf(stderr, "%s mutations are not equal", orig_node->identifier.c_str());
            //        }
            //    }
            //    else 
            //    {
            //        fprintf(stderr, "%s is not equal", orig_node->identifier.c_str());
            //    }
            //}
        }
        else
            fprintf(stderr, "WARNING: %s is already present in the tree.\n", s.c_str());
    }

}
