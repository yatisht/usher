#include "select.hpp"

/*
Functions in this module take a variety of arguments, usually including a MAT
and return a set of samples as a std::vector<std::string>
*/

std::vector<std::string> read_sample_names (std::string sample_filename) {
    std::vector<std::string> sample_names;
    std::ifstream infile(sample_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", sample_filename.c_str());
        exit(1);
    }    
    std::string line;
    bool warned = false;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() != 1 && !warned) {
            fprintf(stderr, "WARNING: Input file %s contains excess columns; ignoring\n", sample_filename.c_str());
            warned = true;
        }
        //remove carriage returns from the input to handle windows os 
        auto sname = words[0];
        if (sname[sname.size()-1] == '\r') {
            sname = sname.substr(0,sname.size()-1);
        }
        sample_names.push_back(std::move(sname));
    }
    infile.close();
    return sample_names;
}

std::vector<std::string> get_clade_samples (MAT::Tree T, std::string clade_name) {
    //fetch the set of sample names associated with a clade name to pass downstream in lieu of reading in a sample file.

    std::vector<std::string> csamples;
    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        if (s->is_leaf()) {
            //check if the clade name is listed anywhere in the annotation vector.
            if (std::find(s->clade_annotations.begin(), s->clade_annotations.end(), clade_name) != s->clade_annotations.end()) {
                csamples.push_back(s->identifier);
            }
        }
    }
    return csamples;
}

std::vector<std::string> get_mutation_samples (MAT::Tree T, std::string mutation_id) {
    //fetch the set of sample names which contain a given mutation.
    //this is a naive implementation parallel to describe::mutation_paths
    std::vector<std::string> good_samples;

    for (auto node: T.get_leaves()) {
        bool assigned = false;
        //first, check if this specific sample has the mutation 
        for (auto m: node->mutations) {
            if (m.get_string() == mutation_id) {
                good_samples.push_back(node->identifier);
                assigned = true;
                break;
            }
        } 
        if (!assigned) {
            std::vector<MAT::Node*> path = T.rsearch(node->identifier);
            //for every ancestor up to the root, check if they have the mutation
            //if they do, break, save the name, move to the next sample
            for (auto anc_node: path) {
                if (!assigned) {
                    for (auto m: anc_node->mutations) {
                        good_samples.push_back(node->identifier);
                        assigned = true;
                        break;                
                    }
                } else {
                    break;
                }
            }
        }
    }
    return good_samples;
}
