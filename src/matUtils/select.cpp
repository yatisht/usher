#include "select.hpp"
#include <random>
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

std::vector<std::string> get_clade_samples (MAT::Tree* T, std::string clade_name) {
    //fetch the set of sample names associated with a clade name to pass downstream in lieu of reading in a sample file.

    std::vector<std::string> csamples;
    auto dfs = T->depth_first_expansion();
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
                        csamples = T->get_leaves_ids(s->identifier);
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

std::vector<std::string> get_mutation_samples (MAT::Tree* T, std::string mutation_id) {
    //fetch the set of sample names which contain a given mutation.
    //this is a naive implementation parallel to describe::mutation_paths
    std::vector<std::string> good_samples;

    for (auto node: T->get_leaves()) {
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
            std::vector<MAT::Node*> path = T->rsearch(node->identifier);
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

std::vector<std::string> get_parsimony_samples (MAT::Tree* T, int max_parsimony) {
    //simple selection- get samples which have less than X parsimony score (e.g. branch length)
    std::vector<std::string> good_samples;
    auto dfs = T->get_leaves();
    for (auto n: dfs) {
        if (n->mutations.size() <= static_cast<size_t>(max_parsimony)) {
            good_samples.push_back(n->identifier);
        }
    }
    return good_samples;
}

std::vector<std::string> get_clade_representatives(MAT::Tree* T, size_t samples_per_clade) {
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
    std::unordered_set<std::string> rep_samples;
    //expand and identify clades seen
    std::unordered_set<std::string> clades_seen;
    static std::random_device rd;
    static std::mt19937 gen(rd());

    auto bfs = T->breadth_first_expansion();
    for (auto n: bfs) {
        std::string curpath;
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                //if this is a new clade annotation, we need its info
                //the first time any new annotation is encountered in an expansion, its the root of that lineage
                if (clades_seen.find(ann) == clades_seen.end()) {
                    clades_seen.insert(ann);
                    std::vector<std::string> leaf_ids = T->get_leaves_ids(n->identifier);
                    if (leaf_ids.size() <= samples_per_clade) {
                        // Add all leaves
                        rep_samples.insert(leaf_ids.begin(), leaf_ids.end());
                    } else {
                        // Randomly select leaves; keep trying if a leaf has already been selected,
                        // but don't keep trying forever in case there just aren't enough
                        // unselected leaves.
                        size_t added = 0, already_selected = 0;
                        std::uniform_int_distribution<> distrib(0, leaf_ids.size() - 1);
                        while (added < samples_per_clade && already_selected < samples_per_clade) {
                            int ix = distrib(gen);
                            if (rep_samples.find(leaf_ids[ix]) == rep_samples.end()) {
                                rep_samples.insert(leaf_ids[ix]);
                                added++;
                            } else {
                                already_selected++;
                            }
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
    std::vector<std::string> rep_sample_vec;
    rep_sample_vec.reserve(rep_samples.size());
    rep_sample_vec.insert(rep_sample_vec.end(), rep_samples.begin(), rep_samples.end());
    return rep_sample_vec;
}

std::vector<std::string> sample_intersect (std::vector<std::string> samples, std::vector<std::string> nsamples) {
    //helper function to get the intersection of two sample identifier vectors
    //used when chaining together other select functions
    assert (samples.size() > 0);
    assert (nsamples.size() > 0);
    std::vector<std::string> inter_samples;
    for (auto s: samples) {
        if (std::find(nsamples.begin(), nsamples.end(), s) != nsamples.end()) {
            inter_samples.push_back(s);
        }
    }
    return inter_samples;
}

std::vector<std::string> get_nearby (MAT::Tree* T, std::string sample_id, int number_to_get) {
    //get the nearest X neighbors to sample_id and return them as a vector
    //the simple indexing method is not guaranteed to get the very closest neighbors when the query sample is out near the edge of a large clade
    //unfortunately. so we have to brute force it.
    MAT::Node* last_anc = T->get_node(sample_id);
    assert (number_to_get > 0);
    if (last_anc == NULL) {
        fprintf(stderr, "ERROR: %s is not present in the tree!\n", sample_id.c_str() );
    }
    std::vector<std::string> leaves_to_keep;
    for (auto anc: T->rsearch(sample_id, true)) {
        int num_leaves = static_cast<int>(T->get_num_leaves(anc));
        if (num_leaves <= number_to_get) {
            last_anc = anc;
            continue;
        }

        if (num_leaves > number_to_get) {
            struct NodeDist {
                MAT::Node* node;
                uint32_t num_mut;

                NodeDist(MAT::Node* n, uint32_t d) {
                    node = n;
                    num_mut = d;
                }

                inline bool operator< (const NodeDist& n) const {
                    return ((*this).num_mut < n.num_mut);
                }
            };

            for (auto l: T->get_leaves(last_anc->identifier)) {
                leaves_to_keep.emplace_back(l->identifier);
            }

            std::vector<NodeDist> node_distances;
            for (auto l: T->get_leaves(anc->identifier)) {
                if (T->is_ancestor(last_anc->identifier, l->identifier)) {
                    continue;
                }

                uint32_t dist = 0;
                for (auto a: T->rsearch(l->identifier, true)) {
                    if (a == anc) {
                        break;
                    }
                    dist += a->mutations.size();
                }

                node_distances.emplace_back(NodeDist(l, dist));
            }

            std::sort(node_distances.begin(), node_distances.end());
            for (auto n: node_distances) {
                if (leaves_to_keep.size() == number_to_get) {
                    break;
                }
                leaves_to_keep.emplace_back(n.node->identifier);
            }
        } else {
            for (auto l: T->get_leaves(anc->identifier)) {
                leaves_to_keep.emplace_back(l->identifier);
            }
        }
        if (static_cast<int>(leaves_to_keep.size()) >= number_to_get) {
            break;
        }
    }
    return leaves_to_keep;
}

std::vector<std::string> get_short_steppers(MAT::Tree* T, std::vector<std::string> samples_to_check, int max_mutations) {
    //for each sample in samples_to_check, this function rsearches along that samples history in the tree
    //if any of the ancestors have greater than max_mutations mutations, then it breaks and marks that sample as a toss
    //including the sample itself. It takes a list of samples to check because rsearching is not a super fast process
    //and we may as well be efficient about it when time permits.
    std::vector<std::string> good_samples;
    if (samples_to_check.size() == 0) {
        //if nothing is passed in, then check the whole tree.
        samples_to_check = T->get_leaves_ids();
    }
    for (auto s: samples_to_check) {
        auto n = T->get_node(s);
        //check this sample immediately before spending cycles getting the ancestors
        if (n->mutations.size() > static_cast<size_t>(max_mutations)) {
            continue;
        }
        auto anc_nodes = T->rsearch(s);
        bool badanc = false;
        for (auto an: anc_nodes) {
            if (an->mutations.size() > static_cast<size_t>(max_mutations)) {
                badanc = true;
                break;
            }
        }
        if (!badanc) {
            good_samples.push_back(s);
        }
    }
    return good_samples;
}

std::unordered_map<std::string,std::unordered_map<std::string,std::string>> read_metafile(std::string metainf, std::set<std::string> samples_to_use) {
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
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> metamap;
    std::vector<std::string> keys;

    while (std::getline(infile, line)) {
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
        } else {
            for (size_t i=1; i < words.size(); i++) {
                metamap[keys[i]][words[0]] = words[i];
            }
        }
    }
    infile.close();
    return metamap;
}

std::vector<std::string> get_sample_match(MAT::Tree* T, std::string substring) {
    //get the set of samples which match the regular expression pattern and return them.
    //simple enough.
    std::regex pat (substring);
    std::vector<std::string> matchsamples;
    for (auto l: T->get_leaves_ids()) {
        if (std::regex_match(l, pat)) {
            matchsamples.emplace_back(l);
        }
    }
    return matchsamples;
}

std::vector<std::string> fill_random_samples(MAT::Tree* T, std::vector<std::string> current_samples, size_t target_size, bool lca_limit) {
    //expand the current sample selection with random samples until it is the indicated size.
    //alternatively, prune random samples from the selection until it is the indicated size, as necessary.
    std::set<std::string> choices;
    std::srand(std::time(nullptr));
    fprintf(stderr, "Selected sample set is %ld samples with %ld requested subtree size; ", current_samples.size(), target_size);
    if (current_samples.size() > target_size) {
        fprintf(stderr, "removing random samples\n");
        for (size_t i = 0; i < current_samples.size(); i++) {
            //technically, what this implementation is doing is selecting random samples to keep from among the current set.
            auto l = current_samples.begin();
            std::advance(l, std::rand() % current_samples.size());
            choices.insert(*l);
            if (choices.size() >= target_size) {
                break;
            }
        }
    } else if (current_samples.size() < target_size) {
        fprintf(stderr, "filling in with random samples\n");
        std::vector<std::string> all_leaves_ids;
        if (lca_limit) {
            std::string lca = MAT::get_subtree(*T, current_samples).root->identifier;
            std::cerr << "Selecting only samples below " << lca << "\n";
            all_leaves_ids = T->get_leaves_ids(lca);
        } else {
            all_leaves_ids = T->get_leaves_ids();
        }
        choices.insert(current_samples.begin(), current_samples.end());
        for (size_t i = 0; i < all_leaves_ids.size(); i++) {
            auto l = all_leaves_ids.begin();
            std::advance(l, std::rand() % all_leaves_ids.size());
            choices.insert(*l);
            if (choices.size() >= target_size) {
                break;
            }
        }
    } else {
        fprintf(stderr, "continuing\n");
        choices.insert(current_samples.begin(), current_samples.end());
    }
    std::vector<std::string> filled_samples (choices.begin(), choices.end());
    //assert (filled_samples.size() == target_size);
    if ((filled_samples.size() < target_size) && (lca_limit)) {
        fprintf(stderr, "WARNING: Less than the requested number of samples are available beneath the LCA; %ld samples included in output\n", filled_samples.size());
    }
    return filled_samples;
}
