#include "mask.hpp"
#include <random>
#include <algorithm>

po::variables_map parse_mask_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description filt_desc("mask options");
    filt_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-mat,o", po::value<std::string>()->required(),
     "Path to output masked mutation-annotated tree file [REQUIRED]")
    ("simplify,S", po::bool_switch(),
     "Use to automatically remove identifying information from the tree, including all sample names and private mutations.")
    ("restricted-samples,s", po::value<std::string>()->default_value(""),
     "Sample names to restrict. Use to perform masking")
    ("rename-samples,r", po::value<std::string>()->default_value(""),
     "Name of the TSV file containing names of the samples to be renamed and their new names")
    ("mask-mutations,m", po::value<std::string>()->default_value(""),
     "Name of a TSV or CSV containing mutations to be masked in the first column and locations to mask downstream from in the second. If only one column is passed, all instances of that mutation on the tree are masked.")
    ("condense-tree,c", po::bool_switch(),
     "Use to recondense the tree before saving.")
    ("move-nodes,M", po::value<std::string>()->default_value(""),
     "Name of the TSV file containing names of the nodes to be moved and their new parents. Use to move nodes around the tree between paths containing identical sets of mutations.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << filt_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void mask_main(po::parsed_options parsed) {
    po::variables_map vm = parse_mask_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string samples_filename = vm["restricted-samples"].as<std::string>();
    std::string mutations_filename = vm["mask-mutations"].as<std::string>();
    std::string move_nodes_filename = vm["move-nodes"].as<std::string>();
    bool recondense = vm["condense-tree"].as<bool>();
    bool simplify = vm["simplify"].as<bool>();
    std::string rename_filename = vm["rename-samples"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    //check for mutually exclusive arguments
    if ((simplify) & (rename_filename != "")) {
        //doesn't make any sense to rename nodes after you just scrambled their names. Or to rename them, then scramble them.
        fprintf(stderr, "ERROR: Sample renaming and simplification are mutually exclusive operations. Review argument choices\n");
        exit(1);
    }
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
        T.uncondense_leaves();
    }

    // If a restricted samples file was provided, perform masking procedure
    if (samples_filename != "") {
        fprintf(stderr, "Performing Masking...\n");
        restrictSamples(samples_filename, T);
    }
    if (simplify) {
        fprintf(stderr, "Removing identifying information...\n");
        simplify_tree(&T);
    }
    if (mutations_filename != "") {
        fprintf(stderr, "Masking mutations...\n");
        restrictMutationsLocally(mutations_filename, &T);
    }

    // If a rename file was provided, perform renaming procedure
    if (rename_filename != "") {
        fprintf(stderr, "Performing Renaming\n");
        renameSamples(rename_filename, T);
    }

    if (move_nodes_filename != "") {
        moveNodes(move_nodes_filename, &T);
    }

    // Store final MAT to output file
    if (output_mat_filename != "") {
        fprintf(stderr, "Saving Final Tree\n");
        if (recondense) {
            T.collapse_tree();
            T.condense_leaves();
        }
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
    }
}

void simplify_tree(MAT::Tree* T) {
    /*
    This function is intended for the removal of potentially problematic information from the tree while keeping the core structure.
    This renames all samples to arbitrary numbers, similar to internal nodes, and removes sample mutations.
    */
    auto all_leaves = T->get_leaves();
    std::shuffle(all_leaves.begin(), all_leaves.end(), std::default_random_engine(0));
    int rid = 0;
    for (auto l: all_leaves) {
        //only leaves need to have their information altered.
        //remove the mutations first, then change the identifier.
        l->mutations.clear();
        std::stringstream nname;
        //add the l to distinguish them from internal node IDs
        nname << "l" << rid;
        T->rename_node(l->identifier, nname.str());
        rid++;
    }
    auto tree_leaves = T->get_leaves_ids();
    for (auto l1_id: tree_leaves) {
        std::vector<MAT::Node*> polytomy_nodes;
        // Use the same method as T.condense_leaves to find sets of leaves that are identical
        // to their parent nodes.
        auto l1 = T->get_node(l1_id);
        if (l1 == NULL) {
            continue;
        }
        if (l1->mutations.size() > 0) {
            continue;
        }
        for (auto l2: l1->parent->children) {
            if (l2->is_leaf() && (T->get_node(l2->identifier) != NULL) && (l2->mutations.size() == 0)) {
                polytomy_nodes.push_back(l2);
            }
        }
        if (polytomy_nodes.size() > 1) {
            // Leave the first node in the set in the tree, but remove all other identical nodes.
            for (size_t it = 1; it < polytomy_nodes.size(); it++) {
                T->remove_node(polytomy_nodes[it]->identifier, false);
            }
        }
    }
}

void renameSamples (std::string rename_filename, MAT::Tree& T) {

    std::ifstream infile(rename_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the renaming file: %s!\n", rename_filename.c_str());
        exit(1);
    }
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() != 2) {
            fprintf(stderr, "ERROR: Incorrect format for the renaming file: %s!\n", rename_filename.c_str());
            exit(1);
        }
        if (T.get_node(words[0]) == NULL) {
            fprintf(stderr, "WARNING: Node %s not found in the MAT.\n", words[0].c_str());
        } else {
            fprintf(stderr, "Renaming node %s to %s.\n", words[0].c_str(), words[1].c_str());
            T.rename_node(words[0], words[1]);
        }
    }
}

bool match_mutations(MAT::Mutation* target, MAT::Mutation* query) {
    /*
    This function is used to compare two mutations. An N in the target is treated as a match to any other IUPAC base.
    */
    if (target->position != query->position) {
        return false;
    }
    if (target->ref_nuc != 0b1111) {
        if (target->ref_nuc != query->ref_nuc) {
            return false;
        }
    }
    if (target->mut_nuc != 0b1111) {
        if (target->mut_nuc != query->mut_nuc) {
            return false;
        }
    }
    return true;
}
void restrictMutationsLocally (std::string mutations_filename, MAT::Tree* T, bool global) {
    std::ifstream infile(mutations_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", mutations_filename.c_str());
        exit(1);
    }
    std::string line;
    char delim = '\t';
    if (mutations_filename.find(".csv\0") != std::string::npos) {
        delim = ',';
    }
    std::string rootid = T->root->identifier;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        if (line[line.size()-1] == '\r') {
            line = line.substr(0, line.size()-1);
        }
        MAT::string_split(line, delim, words);
        std::string target_node;
        std::string target_mutation;
        if ((words.size() == 1) || (global)) {
            //std::cerr << "Masking mutations globally.\n";
            target_mutation = words[0];
            target_node = rootid;
        } else {
            target_mutation = words[0];
            target_node = words[1];
        }
        MAT::Mutation* mutobj = MAT::mutation_from_string(target_mutation);
        size_t instances_masked = 0;
        MAT::Node* rn = T->get_node(target_node);
        if (rn == NULL) {
            fprintf(stderr, "ERROR: Internal node %s requested for masking does not exist in the tree. Exiting\n", target_node.c_str());
            exit(1);
        }
        // fprintf(stderr, "Masking mutation %s below node %s\n", ml.first.c_str(), ml.second.c_str());
        for (auto n: T->depth_first_expansion(rn)) {
            std::vector<MAT::Mutation> nmuts;
            for (auto& mut: n->mutations) {
                if (match_mutations(mutobj, &mut)) {
                    instances_masked++;
                } else {
                    nmuts.push_back(mut);
                }
            }
            n->mutations = nmuts;
        }
    }
    infile.close();
    T->collapse_tree();
}

void restrictSamples (std::string samples_filename, MAT::Tree& T) {
    std::ifstream infile(samples_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the restricted samples file: %s!\n", samples_filename.c_str());
        exit(1);
    }
    std::unordered_set<std::string> restricted_samples;
    std::string sample;
    while (std::getline(infile, sample)) {
        fprintf(stderr, "Checking for Sample %s\n", sample.c_str());
        if (T.get_node(sample) == NULL) {
            fprintf(stderr, "ERROR: Sample missing in input MAT!\n");
            std::cerr << std::endl;
            exit(1);
        }
        restricted_samples.insert(std::move(sample));
    }
    assert (restricted_samples.size() > 0);
    // Set of nodes rooted at restricted samples
    std::unordered_set<MAT::Node*> restricted_roots;
    std::unordered_map<std::string, bool> visited;
    for (auto s: restricted_samples) {
        visited[s] = false;
    }
    for (auto cn: T.breadth_first_expansion()) {
        auto s = cn->identifier;
        if (restricted_samples.find(s) == restricted_samples.end()) {
            continue;
        }
        if (visited[s]) {
            continue;
        }
        auto curr_node = T.get_node(s);
        for (auto n: T.rsearch(s)) {
            bool found_unrestricted = false;
            for (auto l: T.get_leaves_ids(n->identifier)) {
                if (restricted_samples.find(l)  == restricted_samples.end()) {
                    found_unrestricted = true;
                    break;
                }
            }
            if (!found_unrestricted) {
                for (auto l: T.get_leaves_ids(n->identifier)) {
                    visited[l] = true;
                }
                curr_node = n;
                break;
            }
        }
        restricted_roots.insert(curr_node);
    }

    fprintf(stderr, "Restricted roots size: %zu\n\n", restricted_roots.size());

    // Map to store number of occurences of a mutation in the tree
    std::unordered_map<std::string, int> mutations_counts;
    for (auto n: T.depth_first_expansion()) {
        for (auto mut: n->mutations) {
            if (mut.is_masked()) {
                continue;
            }
            auto mut_string = mut.get_string();
            if (mutations_counts.find(mut_string) == mutations_counts.end()) {
                mutations_counts[mut_string] = 1;
            } else {
                mutations_counts[mut_string] += 1;
            }
        }
    }

    // Reduce mutation counts for mutations in subtrees rooted at
    // restricted_roots. Mutations specific to restricted samples
    // will now be set to 0.
    for (auto r: restricted_roots) {
        //fprintf(stdout, "At restricted root %s\n", r->identifier.c_str());
        for (auto n: T.depth_first_expansion(r)) {
            for (auto mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                mutations_counts[mut_string] -= 1;
            }
        }
    }

    for (auto r: restricted_roots) {
        for (auto n: T.depth_first_expansion(r)) {
            for (auto& mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                if (mutations_counts[mut_string] == 0) {
                    fprintf(stderr, "Masking mutation %s at node %s\n", mut_string.c_str(), n->identifier.c_str());
                    mut.position = -1;
                    mut.ref_nuc = 0;
                    mut.par_nuc = 0;
                    mut.mut_nuc = 0;
                }
            }
        }
    }
}

void moveNodes (std::string node_filename, MAT::Tree* T) {
    // Function to move nodes between two identical placement paths. That is, move the target node so that is a child of the indicated new parent node,
    // but the current placement and new placement must involve exactly the same set of mutations for the move to be allowed.
    // Takes the path to a two-column tsv containing the names of nodes to be moved and the parents to move them to.
    std::ifstream infile(node_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the moving file: %s!\n", node_filename.c_str());
        exit(1);
    }
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words.size() != 2) {
            fprintf(stderr, "ERROR: Incorrect format for the moving file: %s!\n", node_filename.c_str());
            exit(1);
        }

        MAT::Node* mn = T->get_node(words[0]);
        MAT::Node* np = T->get_node(words[1]);
        if (mn == NULL) {
            fprintf(stderr, "ERROR: Node %s does not exist in the tree. Exiting\n", words[0].c_str());
            exit(1);
        }
        if (np == NULL) {
            fprintf(stderr, "ERROR: Node %s does not exist in the tree. Exiting\n", words[1].c_str());
            exit(1);
        }
        if (np->is_leaf()) {
            fprintf(stderr, "ERROR: Node %s is a leaf and therefore cannot be a parent. Exiting\n", words[1].c_str());
            exit(1);
        }
        if (mn->parent == np) {
            fprintf(stderr, "ERROR: Node %s is already a child of %s. Exiting\n", words[0].c_str(), words[1].c_str());
            exit(1);
        }
        //accumulate the set of mutations belonging to the current and the new placement
        //not counting mutations belonging to the target, but counting ones to the putative new parent.
        std::unordered_set<std::string> curr_mutations;
        std::unordered_set<std::string> new_mutations;
        for (auto an: T->rsearch(mn->identifier, false)) {
            for (auto mut: an->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                MAT::Mutation mut_opposite = mut.copy();
                //we check for whether this mutation is going to negate with something in the set
                //by identifying its opposite and checking whether the opposite is already present on the traversal.
                mut_opposite.par_nuc = mut.mut_nuc;
                mut_opposite.mut_nuc = mut.par_nuc;
                auto cml = curr_mutations.find(mut_opposite.get_string());
                if (cml != curr_mutations.end()) {
                    curr_mutations.erase(cml);
                } else {
                    curr_mutations.insert(mut.get_string());
                }
            }
        }
        for (auto an: T->rsearch(np->identifier, true)) {
            for (auto mut: an->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                MAT::Mutation mut_opposite = mut.copy();
                mut_opposite.par_nuc = mut.mut_nuc;
                mut_opposite.mut_nuc = mut.par_nuc;
                auto cml = new_mutations.find(mut_opposite.get_string());
                if (cml != new_mutations.end()) {
                    new_mutations.erase(cml);
                } else {
                    new_mutations.insert(mut.get_string());
                }
            }
        }
        if (curr_mutations != new_mutations) {
            fprintf(stderr, "ERROR: The current and new placement paths do not involve the same set of mutations. Exiting\n");
            exit(1);
        }
        //we can now proceed with the move.
        T->move_node(mn->identifier, np->identifier);
        fprintf(stderr, "Move of node %s to node %s successful.\n", words[0].c_str(), words[1].c_str());
    }
    fprintf(stderr, "All requested moves complete.\n");
}