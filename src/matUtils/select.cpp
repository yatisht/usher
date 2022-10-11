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
        if (words.size() > 1 && !warned) {
            fprintf(stderr, "WARNING: Input file %s contains excess columns; ignoring\n", sample_filename.c_str());
            warned = true;
        } else if (words.size() == 0) {
            fprintf(stderr, "WARNING: Empty line in input file %s; ignoring\n", sample_filename.c_str());
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
    MAT::Mutation* mutobj = MAT::mutation_from_string(mutation_id);
    for (auto node: T->get_leaves()) {
        bool assigned = false;
        //first, check if this specific sample has the mutation
        //"having the mutation" specifically means that it mutated to this base at this location
        //at the most recent time this location mutated
        for (auto m: node->mutations) {
            if (mutobj->position == m.position) {
                if (mutobj->mut_nuc == m.mut_nuc) {
                    good_samples.push_back(node->identifier);
                }
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
                        if (mutobj->position == m.position) {
                            //if there's a mutation at this position, break out of the loop regardless
                            //but only store it if its the correct base
                            if (mutobj->mut_nuc == m.mut_nuc) {
                                good_samples.push_back(node->identifier);
                            }
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

std::vector<std::string> get_parsimony_samples (MAT::Tree* T, std::vector<std::string> samples_to_check, int max_parsimony) {
    //simple selection- get samples which have less than X parsimony score (e.g. branch length)
    if (samples_to_check.size() == 0) {
        //if nothing is passed in, then check the whole tree.
        samples_to_check = T->get_leaves_ids();
    }
    std::vector<std::string> good_samples;
    for (auto sn: samples_to_check) {
        auto n = T->get_node(sn);
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

std::vector<std::string> sample_intersect (std::unordered_set<std::string> samples, std::vector<std::string> nsamples) {
    //helper function to get the intersection of two a sample identifier set and vector
    //used when chaining together other select functions
    assert (samples.size() > 0);
    assert (nsamples.size() > 0);
    std::vector<std::string> inter_samples;
    for (auto s: nsamples) {
        if (samples.find(s) != samples.end()) {
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
                if ((int)leaves_to_keep.size() == number_to_get) {
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

std::vector<std::string> get_short_paths (MAT::Tree* T, std::vector<std::string> samples_to_check, int max_path) {
    std::vector<std::string> good_samples;
    std::unordered_set<std::string> sampleset;
    if (samples_to_check.size() != 0) {
        sampleset.insert(samples_to_check.begin(), samples_to_check.end());
    }
    std::unordered_map<std::string, size_t> path_lengths;
    for (auto n: T->depth_first_expansion()) {
        if (!n->is_leaf()) {
            //if its an internal node, add it to the tracker. Path length is the length to its parent plus its mutations.
            //path length to the root is 0.
            if (n->is_root()) {
                path_lengths[n->identifier] = 0;
            } else {
                path_lengths[n->identifier] = path_lengths[n->parent->identifier] + n->mutations.size();
            }
        } else {
            //if its not an internal node, check to see if the length to its parent plus the length to the leaf is under the maximum.
            if (path_lengths[n->parent->identifier] + n->mutations.size() <= max_path) {
                if (samples_to_check.size() == 0 || sampleset.find(n->identifier) != sampleset.end()) {
                    good_samples.push_back(n->identifier);
                }
            }
        }
    }
    return good_samples;
}

void prune_except_clade_roots(Mutation_Annotated_Tree::Node* node, std::unordered_set<std::string>& samples_to_prune, std::unordered_set<std::string>& nodes_exempt, std::string& chopped_node, const size_t leaf_count_total, const size_t mut_count_total, const double mut_density) {
    // Recursively descend node, adding IDs of all leaf descendants of node to samples_to_prune,
    // unless node is a clade root (in that case print a message for now***) or is in nodes_exempt.
    bool has_clade_annotation = false;
    for (auto cann: node->clade_annotations) {
        if (cann != "") {
            has_clade_annotation = true;
            fprintf(stderr, "Sparing node %s (from %s) because it has clade annotation %s\n", node->identifier.c_str(), chopped_node.c_str(), cann.c_str());
            break;
        }
    }
    if (! has_clade_annotation &&
        nodes_exempt.find(node->identifier) == nodes_exempt.end()) {
        if (node->children.size() == 0) {
            samples_to_prune.insert(node->identifier);
        } else {
            for (auto child: node->children) {
              prune_except_clade_roots(child, samples_to_prune, nodes_exempt, chopped_node, leaf_count_total, mut_count_total, mut_density);
            }
        }
    }
}

void add_nodes_to_set(Mutation_Annotated_Tree::Node* node, std::unordered_set<std::string>& set) {
    // Recursively descend node, adding IDs of all internal nodes to set.
    if (node->children.size() > 0) {
        set.insert(node->identifier);
        for (auto child: node->children) {
            add_nodes_to_set(child, set);
        }
    }
}

void filter_mut_density_helper(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Tree* T, const double max_mut_density, const std::unordered_map<std::string, size_t>& node_leaf_counts, size_t& leaf_count, size_t& mut_count, std::unordered_set<std::string>& samples_to_prune, std::unordered_set<std::string>& nodes_exempt) {
    // Recursively descend node, identifying samples to prune because their branches have too-high
    // mutation density and no redeeming qualities such as being an annotated clade.
    if (node->children.size() > 0) {
        // First, exempt all descendants of clade root nodes with small leaf count from pruning
        for (auto cann: node->clade_annotations) {
            if (cann != "") {
                size_t pre_prune_count = node_leaf_counts.at(node->identifier);
                if (pre_prune_count < 150) {
                    add_nodes_to_set(node, nodes_exempt);
                }
                break;
            }
        }
        // Internal node: sum up counts of children to get mutation density
        bool this_node_exempt = (nodes_exempt.find(node->identifier) != nodes_exempt.end());
        size_t leaf_count_total = 0, mut_count_total = 0, own_mut_count = node->mutations.size();
        for (auto child: node->children) {
            if (!this_node_exempt && child->children.size() == 0) {
                // Special case for leaf child of node that's not already exempt: prune if
                // leaf mutation count exceeds max_mut_density.
                size_t child_mut_count = child->mutations.size();
                if (child_mut_count > max_mut_density) {
                    samples_to_prune.insert(child->identifier);
                } else {
                    leaf_count_total++;
                    mut_count_total += child_mut_count;
                }
            } else {
                size_t leaf_count_child, mut_count_child;
                filter_mut_density_helper(child, T, max_mut_density, node_leaf_counts, leaf_count_child, mut_count_child, samples_to_prune, nodes_exempt);
                leaf_count_total += leaf_count_child;
                mut_count_total += mut_count_child;
            }
        }
        // If all children have been pruned, then this node has effectively been pruned;
        // just return zero counts.
        if (leaf_count_total == 0) {
            leaf_count = 0;
            mut_count = 0;
        } else {
            // Treat mutations at this node specially for purposes of computing mutation density:
            // divide by leaf count to avoid penalizing a long branch to a tight cluster.
            double own_mut_contribution = own_mut_count / (double)leaf_count_total;
            double mut_density = (own_mut_contribution + mut_count_total) / (double)leaf_count_total;
            if (mut_density > max_mut_density) {
                // Remove all descendants of this node (unless exempt) and return zero counts
                prune_except_clade_roots(node, samples_to_prune, nodes_exempt, node->identifier, leaf_count_total, mut_count_total, mut_density);
                leaf_count = 0;
                mut_count = 0;
            } else {
                leaf_count = leaf_count_total;
                // Include own_mut_count in return value without scaling
                mut_count = mut_count_total + own_mut_count;
                // Exempt this node from pruning if it has a very low mut_density
                if (mut_density < 1 / max_mut_density) {
                    nodes_exempt.insert(node->identifier);
                }
            }
        }
    } else {
        // Leaf node: just return counts
        leaf_count = 1;
        mut_count = node->mutations.size();
    }
}

size_t get_leaf_counts(Mutation_Annotated_Tree::Node* node, std::unordered_map<std::string, size_t>& node_leaf_counts) {
    // Recursively descend the tree to find and store leaf counts for each node.
    size_t leaf_count = 0;
    if (node->children.size() == 0) {
        leaf_count = 1;
    } else {
        for (auto child: node->children) {
            leaf_count += get_leaf_counts(child, node_leaf_counts);
        }
    }
    node_leaf_counts.insert({node->identifier, leaf_count});
    return leaf_count;
}

std::unordered_set<std::string> filter_mut_density(Mutation_Annotated_Tree::Tree* T, std::unordered_set<std::string>& samples_to_check, double max_mut_density) {
    // Prune branches that have a too-high mutation density, i.e. ratio of total mutations to number
    // of leaves.  Don't prune annotated clades or branches with very low mutation density.
    std::unordered_set<std::string> samples_to_prune, nodes_exempt;
    std::unordered_map<std::string, size_t> node_leaf_counts;
    get_leaf_counts(T->root, node_leaf_counts);
    size_t leaf_count, mut_count;
    filter_mut_density_helper(T->root, T, max_mut_density, node_leaf_counts, leaf_count, mut_count, samples_to_prune, nodes_exempt);
    std::unordered_set<std::string> samples_to_keep;
    for (auto sample = samples_to_check.begin();  sample != samples_to_check.end();  sample++) {
        if (samples_to_prune.find(*sample) == samples_to_prune.end()) {
            samples_to_keep.insert(*sample);
        }
    }
    return samples_to_keep;
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
                if (samples_to_use.find(words[0]) != samples_to_use.end()) {
                    metamap[keys[i]][words[0]] = words[i];
                }
            }
        }
    }
    infile.close();
    return metamap;
}

std::vector<std::string> get_sample_match(MAT::Tree* T, std::vector<std::string> samples_to_check, std::string substring) {
    //get the set of samples which match the regular expression pattern and return them.
    //simple enough.
    if (samples_to_check.size() == 0) {
        samples_to_check = T->get_leaves_ids();
    }
    std::regex pat (substring);
    std::vector<std::string> matchsamples;
    for (auto l: samples_to_check) {
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

std::vector<std::string> get_mrca_samples(MAT::Tree* T, std::vector<std::string> current_samples) {
    //get the subtree, get its root, get all descendents from that node id in the original tree (node ids are retained until storage in a pb).
    auto mrca = MAT::get_subtree(*T, current_samples).root->identifier;
    std::vector<std::string> mrca_samples = T->get_leaves_ids(mrca);
    return mrca_samples;
}

void closest_samples_dfs(MAT::Node *node, MAT::Node *target, size_t path_length, size_t max_path_length, std::vector<std::pair<MAT::Node *, size_t>> &leaves, bool fixed_k) {

    if (path_length > max_path_length) {
        return;
    }
    for (auto child : node->children) {
        if (child->is_leaf()) {
            if (fixed_k && path_length + child->mutations.size() <= max_path_length) {

                leaves.push_back(std::make_pair(child, path_length + child->mutations.size()));
            } else if (!fixed_k) {
                leaves.push_back(std::make_pair(child, path_length + child->mutations.size()));
            }
        } else {
            closest_samples_dfs(child, target, path_length + child->mutations.size(), max_path_length, leaves, fixed_k);
        }
    }
}

std::pair<std::vector<std::string>, size_t> get_closest_samples(MAT::Tree* T, std::string nid, bool fixed_k, size_t k) {
    // Returns a pair with (1) a vector of closest nodes to a target and (2) the distance from the target node
    std::pair<std::vector<std::string>, size_t> closest_samples;

    MAT::Node *target = T->get_node(nid);
    MAT::Node *target_parent = target->parent;
    MAT::Node *curr_target = T->get_node(nid);

    if (!target) {
        fprintf(stderr, "WARNING: Node %s not found in tree\n", nid.c_str());
        return closest_samples;
    }
    MAT::Node *parent = target->parent;

    size_t min_dist = std::numeric_limits<size_t>::max();
    size_t dist_to_orig_parent = 0; // cumulative distance to the parent of the initial target


    bool go_up = true;
    while (go_up && parent) {
        size_t parent_branch_length = parent->mutations.size() + dist_to_orig_parent;
        // make a vector of siblings of the current target.
        // for siblings that are internal nodes, add leaves in the descendant subtree
        // as pseudo-children if they are close enough
        std::vector<std::pair<MAT::Node *, size_t>> children_and_distances;
        std::vector<std::pair<MAT::Node *, size_t>> children_and_distances_fixed_k;

        size_t min_of_sibling_leaves = std::numeric_limits<size_t>::max();
        //collect leaves
        for (auto child : parent->children) {

            if (child->is_leaf()) {

                if (child->identifier == curr_target->identifier) {

                    continue; // skip the target node
                }
                size_t child_branch_length = child->mutations.size();

                if (child_branch_length < min_of_sibling_leaves) {

                    min_of_sibling_leaves = child_branch_length;
                }
            }
        }

        for (auto child : parent->children) {
            if (child->identifier == curr_target->identifier) {
                continue; // skip the target node
            }
            if (child->identifier == target_parent->identifier) {
                continue; // don't go back down path
            }

            size_t dist_so_far = dist_to_orig_parent + target->mutations.size() + child->mutations.size();
            if (!child->is_leaf()) {
                // for internal nodes, descend the tree, adding leaves as they are
                // encountered, restricting path lengths to less than the minimum of
                // the sibling leaves at the current level
                if (fixed_k) {
                    closest_samples_dfs(child, target, dist_so_far, k, children_and_distances_fixed_k, true);
                } else {
                    size_t max_path;
                    if (min_of_sibling_leaves == std::numeric_limits<size_t>::max()) {
                        max_path = min_of_sibling_leaves;
                    } else {
                        max_path = min_of_sibling_leaves + dist_so_far;
                    }
                    closest_samples_dfs(child, target, dist_so_far, max_path, children_and_distances, false);
                }

            } else { // leaf node
                if(fixed_k && dist_so_far <= k) {
                    children_and_distances_fixed_k.push_back(std::make_pair(child, dist_so_far));
                } else if (!fixed_k) {
                    children_and_distances.push_back(std::make_pair(child, dist_so_far));
                }
            }
        }

        if (fixed_k) {
            if (parent_branch_length > k) {
                go_up = false;
            }
            for (std::pair child_and_dist : children_and_distances_fixed_k) {
                MAT::Node *child = child_and_dist.first;
                size_t child_branch_length = child_and_dist.second;
                closest_samples.first.push_back(child->identifier);
                closest_samples.second = 0;
            }

        } else {
            for (std::pair child_and_dist : children_and_distances) {
                // for the siblings of the target node, if any branch lengths
                // are shorter than the path up a level, we can stop
                MAT::Node *child = child_and_dist.first;
                size_t child_branch_length = child_and_dist.second;
                if (child_branch_length < parent_branch_length) {
                    go_up = false;
                }
                if (child_branch_length < min_dist) {
                    min_dist = child_branch_length;
                    closest_samples.first.clear();
                    closest_samples.first.push_back(child->identifier);
                    closest_samples.second = min_dist;

                } else if (child_branch_length == min_dist) {
                    closest_samples.first.push_back(child->identifier);
                }
            }
        }

        curr_target = parent;
        parent = curr_target->parent;
        dist_to_orig_parent = parent_branch_length;

    }
    return closest_samples;
}
