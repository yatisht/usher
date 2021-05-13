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
        std::vector<std::string> canns = s->clade_annotations;
        if (canns.size() > 0) {
            if (canns.size() > 1 || canns[0] != "") {
                //the empty string is the default clade identifier attribute
                //skip entries which are annotated with 1 clade but that clade is empty
                //but don't skip entries which are annotated with 1 clade and its not empty
                //check if this clade root matches the clade name.
                for (auto c: canns) {
                    if (c == clade_name) {
                        //this is the root of the input clade (first one encountered in tree)
                        //get the set of samples descended from this clade root
                        csamples = T.get_leaves_ids(s->identifier);
                        //and break out by returning 
                        return csamples;
                    }
                }
            }
        }
    }
    //if it was never found, return an empty vector.
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
                        if (m.get_string() == mutation_id) {
                            good_samples.push_back(node->identifier);
                            assigned = true;
                            break;     
                        }
                    }
                } else {
                    break;
                }
            }
        }
    }
    // fprintf(stderr, "# of good samples %ld\n", good_samples.size());
    return good_samples;
}

std::vector<std::string> get_parsimony_samples (MAT::Tree T, int max_parsimony) {
    //simple selection- get samples which have less than X parsimony score (e.g. branch length)
    std::vector<std::string> good_samples;
    auto dfs = T.get_leaves();
    for (auto n: dfs) {
        if (n->mutations.size() <= static_cast<size_t>(max_parsimony)) {
            good_samples.push_back(n->identifier);
        }
    }
    return good_samples;
}

std::vector<std::string> get_clade_representatives(MAT::Tree T) {
    timer.Start();
    fprintf(stderr, "Selecting clade representative samples...");
    //get a pair of representative leaves for every clade currently annotated in the tree 
    //and return as a vector of samples
    //specifically, the leaves LCA must be the clade root node
    //to accomplish this, first identify the clade root,
    //then, since the clade root is guaranteed to have at least two children
    //depth first search down the first child until a sample is encountered, record the name
    //then depth first search down the second child until a sample is encountered, record that name
    //and add both samples to the representative output
    std::vector<std::string> rep_samples;
    //expand and identify clades seen
    std::unordered_set<std::string> clades_seen;
    
    auto dfs = T.breadth_first_expansion();
    for (auto n: dfs) {
        std::string curpath;
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                //if this is a new clade annotation, we need its info
                //the first time any new annotation is encountered in an expansion, its the root of that lineage
                if (clades_seen.find(ann) == clades_seen.end()) {
                    clades_seen.insert(ann);
                    //this should always be true
                    assert (n->children.size() > 1);
                    //search down the first child
                    auto first_dfs = T.depth_first_expansion(n->children[0]);
                    for (auto sn: first_dfs) {
                        //pick a sample which hasn't already been selected to represent some other clade
                        //since clades are nested (is this the correct way to handle this?)
                        if (sn->is_leaf() && std::find(rep_samples.begin(), rep_samples.end(), sn->identifier) == rep_samples.end()) {
                            //when a sample is found, grab it and break
                            rep_samples.push_back(sn->identifier);
                            break;
                        }
                    }
                    //and the second
                    auto second_dfs = T.depth_first_expansion(n->children[1]);
                    for (auto sn: second_dfs) {
                        if (sn->is_leaf() && std::find(rep_samples.begin(), rep_samples.end(), sn->identifier) == rep_samples.end()) {
                            rep_samples.push_back(sn->identifier);
                            break;
                        }
                    }
                    //there may be cases where a clade iterates all the way through and every sample which could represent it is already
                    //selected by some other clade, but in that case its perfectly represented anyways, so theres no reason to get more 
                    //samples just for it
                }
            }
        }
    }
    //this relationship should generally be 2-1, though nesting/overlap can make it a little less than 2 to 1 samples to clades
    fprintf(stderr, "%ld samples chosen to represent ", rep_samples.size());
    fprintf(stderr, "%ld unique clades\n", clades_seen.size());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    return rep_samples;
}

std::vector<std::string> sample_intersect (std::vector<std::string> samples, std::vector<std::string> nsamples) {
    //helper function to get the intersection of two sample identifier vectors
    //used when chaining together other select functions
    std::vector<std::string> inter_samples;
    for (auto s: samples) {
        if (std::find(nsamples.begin(), nsamples.end(), s) != nsamples.end()) {
            inter_samples.push_back(s);
        }
    }
    return inter_samples;
}

std::vector<std::string> get_nearby (MAT::Tree T, std::string sample_id, int number_to_get) {
    //get the nearest X neighbors to sample_id and return them as a vector
    //by far the simplest way to do this is to get_leaves_ids and subset out a distance around the index
    //increment the number to get by 1.
    number_to_get++;
    //trying the set of leaves in breadth-first search order
    std::vector<std::string> all_leaves;
    for (auto n: T.depth_first_expansion()) {
        if (n->is_leaf()) {
            all_leaves.push_back(n->identifier);
        }
    }
    //auto all_leaves = T.get_leaves_ids();
    auto target = std::find(all_leaves.begin(), all_leaves.end(), sample_id);
    if (target == all_leaves.cend()) {
        fprintf(stderr, "ERROR: Indicated sample does not exist in the tree!\n");
        exit(1);
    } else if (all_leaves.size() < static_cast<size_t>(number_to_get)) {
        fprintf(stderr, "ERROR: Not enough samples in tree to get neighborhood subtree of requested size!\n");
        exit(1);
    }
    int tindex = std::distance(all_leaves.begin(), target);
    int subset_start = tindex - (number_to_get/2);
    int subset_end = tindex + (number_to_get/2);
    if (subset_start < 0) {
        subset_start = 0;
        subset_end = number_to_get;
    } else if (static_cast<size_t>(subset_end) >= all_leaves.size()) {
        subset_start = all_leaves.size() - number_to_get - 1;
        subset_end = all_leaves.size() - 1;
    }
    //fprintf(stderr, "Start index is %d\n", subset_start);
    //fprintf(stderr, "Stop index is %d\n", subset_end);
    std::vector<std::string> neighborhood_leaves(all_leaves.begin() + subset_start, all_leaves.begin() + subset_end);
    return neighborhood_leaves;
}

std::vector<std::string> get_short_steppers(MAT::Tree T, std::vector<std::string> samples_to_check, int max_mutations) {
    //for each sample in samples_to_check, this function rsearches along that samples history in the tree
    //if any of the ancestors have greater than max_mutations mutations, then it breaks and marks that sample as a toss
    //including the sample itself. It takes a list of samples to check because rsearching is not a super fast process
    //and we may as well be efficient about it when time permits.
    std::vector<std::string> good_samples;
    if (samples_to_check.size() == 0) {
        //if nothing is passed in, then check the whole tree.
        samples_to_check = T.get_leaves_ids();
    }
    for (auto s: samples_to_check) {
        auto n = T.get_node(s);
        //check this sample immediately before spending cycles getting the ancestors
        if (n->mutations.size() > static_cast<size_t>(max_mutations)) {
            continue;
        }
        auto anc_nodes = T.rsearch(s);
        bool badanc = false;
        for (auto an: anc_nodes) {
            if (an->mutations.size() > static_cast<size_t>(max_mutations)) {
                badanc = true;
                continue;
            }
        }
        if (!badanc) {
            good_samples.push_back(s);
        }
    }
    return good_samples;
}

std::map<std::string,std::map<std::string,std::string>> read_metafile(std::string metainf) {
    std::ifstream infile(metainf);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", metainf.c_str());
        exit(1);
    }    
    std::string line;
    bool first = true;
    char delim = '\t';
    if (metainf.find(".csv\0") != std::string::npos) {
        delim = ',';
    }
    std::map<std::string,std::map<std::string,std::string>> metamap;
    std::vector<std::string> keys;

    while (std::getline(infile, line)) {
        std::cerr << "reading\n";
        std::vector<std::string> words;
        if (line[line.size()-1] == '\r') {
            line = line.substr(0, line.size()-1);
        }
        MAT::string_split(line, delim, words);
        if (first) {
            for (auto w: words) {
                keys.push_back(w);
            }
            first = false;
            std::cerr << keys[0].c_str() << ',' << keys[1].c_str() << "\n";
        } else {
            for (size_t i=1; i < words.size(); i++) {
                std::cerr << words[i];
                metamap[keys[i]][words[0]] = words[i];
            }
        }
    }
    infile.close();
    return metamap;
}