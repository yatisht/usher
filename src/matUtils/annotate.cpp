#include "annotate.hpp"


po::variables_map parse_annotate_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description ann_desc("annotate options");
    ann_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-mat,o", po::value<std::string>()->required(),
     "Path to output processed mutation-annotated tree file [REQUIRED]")
    ("clade-names,c", po::value<std::string>()->default_value(""),
     "Path to a tsv file containing clade assignments of samples. An algorithm automatically locates and annotates clade root nodes.  If a clade has been found already using --clade-paths or appears in --clade-mutations then the clade's assignments in this file are ignored.")
    ("clade-to-nid,C", po::value<std::string>()->default_value(""),
     "Path to a tsv file mapping clades to their respective internal node identifiers.  Cannot be used in conjunction with --clade-paths, --clade-mutations or --clade-names.")
    ("clade-paths,P", po::value<std::string>()->default_value(""),
     "Path to a tsv file mapping clades to mutation paths.  A warning is printed if a path is not found exactly as specified.  Format is the same as the first and third columns of the output of matUtils extract --clade-paths.")
    ("clade-mutations,M", po::value<std::string>()->default_value(""),
     "Path to a tsv file mapping clades to sets of mutations (separated by spaces, commas and/or >s) which will be used instead of extracting mutations from samples named in the --clade-names file.  This method is more robust to small changes in paths, but significantly slower than --clade-paths.  If used together with --clade-paths, --clade-paths takes precedence.  If used together with --clade-names, this takes precedence.")
    ("allele-frequency,f", po::value<float>()->default_value(0.8),
     "Minimum allele frequency in input samples for finding the best clade root. Used only with -c")
    ("mask-frequency,m", po::value<float>()->default_value(0.2),
     "Minimum allele frequency below -f in input samples that should be masked for finding the best clade root. Used only with -c")
    ("set-overlap,s", po::value<float>()->default_value(0.6),
     "Minimum fraction of the clade samples that should be desecendants of the assigned clade root. Used only with -c")
    ("clip-sample-frequency,p", po::value<float>()->default_value(0.1),
     "Maximum proportion of samples in a branch that are exemplars from -c to consider when sorting candidate clade root nodes")
    ("clear-current,l", po::bool_switch(),
     "Use to remove current annotations before applying new values.")
    ("output-directory,d", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("write-mutations,u", po::value<std::string>()->default_value(""),
     "Write a tsv listing each clade and the mutations found in at least [-f] of samples. Used only with -c")
    ("write-details,D", po::value<std::string>()->default_value(""),
     "Write a tsv with details about the nodes considered for each clade root. Used only with -c")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(ann_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << ann_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

std::string make_out_dir(const std::string& dir_prefix) {
    // Create output directory if necessary and return its canonicalized path
    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory %s.\n\n", dir_prefix.c_str());
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    std::string canonical_dir_prefix = path.generic_string();
    return canonical_dir_prefix;
}

std::string add_out_dir(const std::string& dir_prefix, const std::string& filename) {
    // Prepend dir_prefix to filename unless filename is an absolute path or empty
    std::string filepath = filename;
    if (filename != "" && filename[0] != '/') {
        filepath = dir_prefix + "/" + filename;
    }
    return filepath;
}
FILE* must_open(const std::string& filename, const char* mode) {
    FILE* file = fopen(filename.c_str(), mode);
    if (! file) {
        fprintf(stderr, "ERROR: Could not open file '%s' with mode '%s'\n", filename.c_str(), mode);
        exit(1);
    }
    return file;
}

void annotate_main(po::parsed_options parsed) {
    //the annotate subcommand assigns samples to lineages and saves it as MAT metadata
    po::variables_map vm = parse_annotate_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = make_out_dir(vm["output-directory"].as<std::string>());
    std::string output_mat_filename = add_out_dir(dir_prefix, vm["output-mat"].as<std::string>());
    std::string clade_filename = vm["clade-names"].as<std::string>();
    std::string clade_to_nid_filename = vm["clade-to-nid"].as<std::string>();
    std::string clade_paths_filename = vm["clade-paths"].as<std::string>();
    std::string clade_mutations_filename = vm["clade-mutations"].as<std::string>();
    std::string mutations_filename = add_out_dir(dir_prefix, vm["write-mutations"].as<std::string>());
    std::string details_filename = add_out_dir(dir_prefix, vm["write-details"].as<std::string>());
    bool clear_current = vm["clear-current"].as<bool>();
    float allele_frequency = vm["allele-frequency"].as<float>();
    float mask_frequency = vm["mask-frequency"].as<float>();
    float set_overlap = vm["set-overlap"].as<float>();
    float clip_sample_frequency = vm["clip-sample-frequency"].as<float>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    int clade_spec_count = 0;
    if (clade_filename != "" || clade_mutations_filename != "" || clade_paths_filename != "") {
        clade_spec_count++;
    }
    if (clade_to_nid_filename != "") {
        clade_spec_count++;
    }
    if (clade_spec_count != 1) {
        fprintf(stderr, "ERROR: must specify either --clade-to-nid or [--clade-names and/or --clade-mutations and/or --clade-paths]!\n");
        exit(1);
    }

    tbb::task_scheduler_init init(num_threads);

    // Load input MAT and uncondense tree
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    timer.Start();
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
        T.uncondense_leaves();
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    fprintf(stderr, "Annotating Lineage Root Nodes\n");
    if (clade_filename != "" || clade_mutations_filename != "" || clade_paths_filename != "") {
        assignLineages(T, clade_filename, clade_mutations_filename, clade_paths_filename, allele_frequency, mask_frequency, set_overlap, clip_sample_frequency, clear_current, mutations_filename, details_filename);
    } else {
        assignLineages(T, clade_to_nid_filename, clear_current);
    }

    //condense_leaves() expects some samples to ignore. We don't have any such samples
    //this would be space to add an additional argument containing samples to not recondense
    //for now, just recondense everything
    fprintf(stderr, "Recondensing leaves\n");
    timer.Start();
    T.condense_leaves(std::vector<std::string>());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    // Store final MAT to output file
    if (output_mat_filename != "") {
        fprintf(stderr, "Saving Final Tree to %s\n", output_mat_filename.c_str());
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
    }
}

void init_annotations(std::vector<MAT::Node*>& dfs, bool clear_current) {
    size_t total_nodes = dfs.size();
    for (size_t idx = 0; idx < total_nodes; idx++) {
        //remove the current annotations if unwanted.
        if (clear_current) {
            dfs[idx]->clade_annotations.clear();
        }
        // Add a new entry in annotations
        dfs[idx]->clade_annotations.emplace_back("");
    }
}

void assignLineages (MAT::Tree& T, const std::string& clade_to_nid_filename, bool clear_current) {
    auto dfs = T.depth_first_expansion();
    init_annotations(dfs, clear_current);
    size_t num_annotations = T.get_num_annotations();

    std::ifstream infile(clade_to_nid_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the clade to node id assignment file: %s!\n", clade_to_nid_filename.c_str());
        exit(1);
    }
    std::string line;

    fprintf(stderr, "Reading clade to node id assignment file and completing assignments.\n");
    timer.Start();
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, '\t', words);
        if ((words.size() > 2) || (words.size() == 1)) {
            fprintf(stderr, "ERROR: Incorrect format for clade to node id assignment file: %s!\n", clade_to_nid_filename.c_str());
            exit(1);
        }
        auto n = T.get_node(words[1]);
        if (n == NULL) {
            fprintf(stderr, "ERROR: Node id %s not found!\n", words[1].c_str());
            exit(1);
        } else if (n->clade_annotations[num_annotations-1] != "") {
            fprintf(stderr, "WARNING: Assigning clade %s to node %s failed as the node is already assigned to clade %s!\n",
                    words[0].c_str(), words[1].c_str(), n->clade_annotations[num_annotations-1].c_str());

        } else {
            n->clade_annotations[num_annotations-1] = words[0];
        }
    }
    infile.close();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void parse_clade_mutations(const std::string& clade_mutations_filename,
                           std::map<std::string, std::vector<MAT::Mutation>>& clade_mutations,
                           std::unordered_set<std::string> clades_already_assigned) {
    std::ifstream infile(clade_mutations_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the clade mutations file: %s!\n", clade_mutations_filename.c_str());
        exit(1);
    }
    std::string line;

    fprintf(stderr, "Reading clade mutations file %s.\n", clade_mutations_filename.c_str());
    bool got_error = false;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, '\t', words);
        // Empty second word (for root node, no mutations) is ignored by string_split; add it back.
        if (words.size() == 1 && line.back() == '\t') {
            words.push_back("");
        }
        if (words.size() != 2) {
            fprintf(stderr, "ERROR: Incorrect format for clade mutations file: %s!  Expected 2 tab-separated words, got %ld (%s)\n", clade_mutations_filename.c_str(), words.size(), line.c_str());
            got_error = true;
            continue;
        }
        std::string clade = words[0];
        // If clade has already been assigned by a method with higher precedence, move on to the next one.
        if (clades_already_assigned.find(clade) != clades_already_assigned.end()) {
            continue;
        }
        // Parse mutations from words[1] and store in clade_mutations[words[0]]
        std::vector<MAT::Mutation> mutations;
        std::vector<int> mut_positions;
        std::vector<std::string> mut_words;
        MAT::string_split(words[1], mut_words);
        for (std::string path_el: mut_words) {
            // Ignore empty string or ">"
            if (path_el == "" || path_el == ">") {
                continue;
            }
            // Might be comma-sep list of mutations
            std::vector<std::string> mut_strings;
            MAT::string_split(path_el, ',', mut_strings);
            for (auto mut_string: mut_strings) {
                MAT::Mutation *mut = MAT::mutation_from_string(mut_string);
                if (mut == NULL) {
                    got_error = true;
                } else {
                    if (std::find(mut_positions.begin(), mut_positions.end(), mut->position) != mut_positions.end()) {
                        fprintf(stderr, "ERROR: Clade %s: position %d used multiple times, must appear only once per clade.\n",
                                clade.c_str(), mut->position);
                        got_error = true;
                        continue;
                    }
                    mut_positions.emplace_back(mut->position);
                    mutations.emplace_back(std::move(*mut));
                }
            }
        }
        clade_mutations[clade] = mutations;
    }
    infile.close();
    if (got_error) {
        fprintf(stderr, "Encountered errors -- exiting.\n");
        exit(1);
    }
}

void parse_clade_names(const std::string& clade_filename,
                       std::map<std::string, std::vector<MAT::Mutation>>& clade_mutations,
                       std::unordered_set<std::string>& clades_already_assigned,
                       std::map<std::string, std::vector<MAT::Node*>>& clade_map,
                       MAT::Tree& uncondensed_T, float min_freq, float mask_freq) {
    std::ifstream infile(clade_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the clade assignment file: %s!\n", clade_filename.c_str());
        exit(1);
    }
    std::string line;

    //Read the clade assignment file line by line and fill up map values
    fprintf(stderr, "Reading clade assignment file %s.\n", clade_filename.c_str());
    timer.Start();
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, '\t', words);
        if ((words.size() > 2) || (words.size() == 1)) {
            fprintf(stderr, "ERROR: Incorrect format for clade assignment file: %s!  "
                    "Expected 2 tab-separated words, got %ld\n", clade_filename.c_str(), words.size());
            exit(1);
        } else if (words.size() == 2) {
            std::string clade = words[0];
            // Skip node lookup & clade_map building if clade mutations file already
            // specified the mutations we want to use for placement, or if the clade was
            // already assigned by a method with higher precedence.
            if (clade_mutations.find(clade) == clade_mutations.end() &&
                clades_already_assigned.find(clade) == clades_already_assigned.end()) {
                auto n = uncondensed_T.get_node(words[1]);
                if (n != NULL) {
                    if (clade_map.find(clade) == clade_map.end()) {
                        clade_map[clade] = std::vector<MAT::Node*>();;
                    }
                    clade_map[clade].emplace_back(n);
                } else {
                    fprintf(stderr, "WARNING: Sample %s not found in input MAT!\n", words[1].c_str());
                }
            }
        }
    }
    infile.close();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    for (auto it: clade_map) {
        const std::string clade = it.first;
        const char* clade_c_str = clade.c_str();
        fprintf(stderr, "Finding mutations in clade %s samples.\n", clade_c_str);

        std::map<std::string, int> mutation_counts;

        static tbb::affinity_partitioner ap;
        tbb::mutex tbb_lock;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, it.second.size()),
        [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i) {
                auto n = it.second[i];
                std::vector<int> anc_positions;
                std::vector<MAT::Mutation> ancestral_mutations;

                // Add ancestral mutations to ancestral mutations. When multiple mutations
                // at same position are found in the path leading from the root to the
                // current node, add only the most recent mutation to the vector
                for (auto anc: uncondensed_T.rsearch(n->identifier, true)) {
                    for (auto m: anc->mutations) {
                        if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                            ancestral_mutations.emplace_back(m);
                            if (!m.is_masked()) {
                                anc_positions.emplace_back(m.position);
                            }
                        }
                    }
                }

                for (auto m: ancestral_mutations) {
                    if (m.ref_nuc != m.mut_nuc) {
                        std::string mut_string = m.chrom + "\t" + std::to_string(m.ref_nuc) + "\t" +
                                                 std::to_string(m.position) + "\t" + std::to_string(m.mut_nuc);
                        tbb_lock.lock();
                        if (mutation_counts.find(mut_string) == mutation_counts.end()) {
                            mutation_counts[mut_string] = 1;
                        } else {
                            mutation_counts[mut_string] += 1;
                        }
                        tbb_lock.unlock();
                    }
                }
            }
        }, ap);

        clade_mutations[clade] = std::vector<MAT::Mutation>();
        for (auto mc: mutation_counts) {
            if (static_cast<float>(mc.second)/it.second.size() >= min_freq) {
                std::vector<std::string> words;
                MAT::string_split(mc.first, words);
                MAT::Mutation m;
                m.chrom = words[0];
                m.ref_nuc = static_cast<int8_t>(std::stoi(words[1]));
                m.par_nuc = m.ref_nuc;
                m.position = std::stoi(words[2]);
                m.mut_nuc = static_cast<int8_t>(std::stoi(words[3]));
                clade_mutations[clade].emplace_back(m);
            } else if (static_cast<float>(mc.second)/it.second.size() >= mask_freq) {
                std::vector<std::string> words;
                MAT::string_split(mc.first, words);
                MAT::Mutation m;
                m.chrom = words[0];
                m.ref_nuc = static_cast<int8_t>(std::stoi(words[1]));
                m.par_nuc = m.ref_nuc;
                m.position = std::stoi(words[2]);
                m.mut_nuc = MAT::get_nuc_id('N');
                clade_mutations[clade].emplace_back(m);
            }

        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

std::string node_mutations_string(MAT::Node *node) {
    std::string s = node->identifier + ":";
    for (size_t ix = 0;  ix < node->mutations.size();  ix++) {
        if (ix > 0) {
            s += ',';
        }
        s += node->mutations[ix].get_string();
    }
    return s;
}

std::string path_to_node(MAT::Tree& T, MAT::Node* node) {
    std::vector<std::string> muts;
    std::string s = node_mutations_string(node);
    muts.emplace_back(std::move(s));
    for (auto anc: T.rsearch(node->identifier)) {
        s = node_mutations_string(anc);
        muts.emplace_back(std::move(s));
    }
    std::reverse(muts.begin(), muts.end());
    std::string path = "";
    for (size_t m = 0; m < muts.size(); m++) {
        if (m > 0) {
            path += " > ";
        }
        path += muts[m];
    }
    return path;
}

void write_mutations(FILE* f, const std::vector<MAT::Mutation>& mutations, const bool write_unmasked = true, const bool write_masked = true) {
    bool gotOne = false;
    for (size_t i = 0; i < mutations.size(); i++) {
        std::string mut = mutations[i].get_string();
        bool is_masked = (mut[mut.size()-1] == 'N');
        if ((write_unmasked && !is_masked) || (write_masked && is_masked)) {
            if (gotOne) {
                fprintf(f, ", ");
            }
            fprintf(f, "%s", mut.c_str());
            gotOne = true;
        }
    }
}

void get_freq_overlap(MAT::Tree& T, MAT::Node* node, std::vector<MAT::Node*>clade_samples, float* freq, float* overlap) {
    static tbb::affinity_partitioner ap;
    size_t num_desc = 0;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, clade_samples.size()),
    [&](const tbb::blocked_range<size_t> r) {
        for (size_t i=r.begin(); i<r.end(); ++i) {
            if (T.is_ancestor(node->identifier, clade_samples[i]->identifier)) {
                __sync_fetch_and_add(&num_desc, 1.0);
            }
        }
    }, ap);

    *freq = (static_cast<float>(num_desc) / T.get_num_leaves(node));
    *overlap = (static_cast<float>(num_desc) / clade_samples.size());
}

void assignLineages (MAT::Tree& T, const std::string& clade_filename,
                     const std::string& clade_mutations_filename,
                     const std::string& clade_paths_filename, float min_freq, float mask_freq,
                     float set_overlap, float clip_sample_frequency, bool clear_current,
                     const std::string& mutations_filename, const std::string& details_filename) {
    static tbb::affinity_partitioner ap;

    fprintf(stderr, "Copying tree with uncondensed leaves.\n");
    timer.Start();
    auto uncondensed_T = MAT::get_tree_copy(T);
    uncondensed_T.uncondense_leaves();

    auto dfs = T.depth_first_expansion();
    size_t total_nodes = dfs.size();
    FILE *mutations_file = NULL, *details_file = NULL;
    if (mutations_filename != "") {
        fprintf(stderr, "Writing clade root node mutations to file %s\n", mutations_filename.c_str());
        mutations_file = must_open(mutations_filename, "w");
        fprintf(mutations_file, "clade\tmutations\n");
    }
    if (details_filename != "") {
        fprintf(stderr, "Writing details to file %s\n", details_filename.c_str());
        details_file = must_open(details_filename.c_str(), "w");
        fprintf(details_file, "clade\tmutations\tmasked_mutations\tnode:freq:overlap\t"
                "already_assigned\tfinal_overlap\texemplar_count\tbest_node_path\n");
    }


    init_annotations(dfs, clear_current);
    std::unordered_map<std::string, size_t> dfs_idx;
    for (size_t idx = 0; idx < total_nodes; idx++) {
        dfs_idx[dfs[idx]->identifier] = idx;
    }
    size_t num_annotations = T.get_num_annotations();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    std::unordered_set<std::string> clades_already_assigned;
    if (clade_paths_filename != "") {
        assignLineagesFromPaths(T, clade_paths_filename, clades_already_assigned);
    }

    std::map<std::string, std::vector<MAT::Mutation>> clade_mutations_map;
    std::map<std::string, std::vector<MAT::Node*>> clade_map;
    if (clade_mutations_filename != "") {
        parse_clade_mutations(clade_mutations_filename, clade_mutations_map, clades_already_assigned);
    }
    if (clade_filename != "") {
        parse_clade_names(clade_filename, clade_mutations_map, clades_already_assigned, clade_map,
                          uncondensed_T, min_freq, mask_freq);

    }

    struct Node_freq {
        size_t best_j;
        float freq;
        float overlap;
        float freq_clipped;
        Node_freq (size_t a, float b, float c, float clip_freq) {
            best_j = a;
            freq = b;
            overlap = c;
            freq_clipped = (freq > clip_freq) ? clip_freq : freq;
        }
        // To sort by highest product of clipped freq and overlap*overlap
        inline bool operator< (const Node_freq& n) const {
            return ((this->freq_clipped * this->overlap * this->overlap) > (n.freq_clipped * n.overlap * n.overlap));
        }
    };

    struct Clade_Assignments {
        std::string clade_name;
        size_t clade_size;
        std::vector<Node_freq> best_node_frequencies;
        std::vector<MAT::Mutation> mutations;
        Clade_Assignments(std::string name, size_t sz, std::vector<MAT::Mutation> muts) {
            clade_name = name;
            clade_size = sz;
            mutations = muts;
        }
        // To sort with fewest best nodes and largest clade size first when samples have been given,
        // but prioritize clades whose mutations have been specified directly (--clade-mutations)
        // have been specified (--clade-names).  Those have clade_size=0 (no samples).
        inline bool operator< (const Clade_Assignments& c) const {
            if (this->clade_size == 0 && c.clade_size > 0) {
                return true;
            } else if (this->clade_size > 0 && c.clade_size == 0) {
                return false;
            } else {
                return ((this->best_node_frequencies.size() < c.best_node_frequencies.size()) ||
                        ((this->best_node_frequencies.size() == c.best_node_frequencies.size()) && (this->clade_size > c.clade_size)));
            }
        }
    };

    std::vector<Clade_Assignments> clade_assignments;

    size_t curr_idx = 0;
    for (auto it: clade_mutations_map) {
        const std::string clade = it.first;
        const char* clade_c_str = clade.c_str();
        std::vector<MAT::Mutation>& clade_mutations = it.second;
        // Mutations need to be sorted by position before placement
        std::sort(clade_mutations.begin(), clade_mutations.end());
        fprintf(stderr, "Mutations above the specified frequency in clade %s: ", clade_c_str);
        write_mutations(stderr, clade_mutations);
        fputc('\n', stderr);
        if (mutations_file) {
            fprintf(mutations_file, "%s\t", clade_c_str);
            write_mutations(mutations_file, clade_mutations, true, false);
            fputc('\n', mutations_file);
        }

        fprintf(stderr, "Finding best node for clade %s.\n", clade_c_str);
        timer.Start();
        size_t best_node_num_leaves = 0;
        int best_set_difference = 1e9;
        size_t best_j = 0;
        bool best_node_has_unique = false;

        std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
        std::vector<std::vector<MAT::Mutation>> node_imputed_mutations(total_nodes);

        std::vector<bool> node_has_unique(total_nodes, false);
        std::vector<size_t> best_j_vec;

        size_t num_best = 1;
        MAT::Node* best_node = T.root;
        best_j_vec.emplace_back(0);

        // Find the best placement(s)
        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t k=r.begin(); k<r.end(); ++k) {
                mapper2_input inp;
                inp.T = &T;
                inp.node = dfs[k];
                inp.missing_sample_mutations = &clade_mutations;
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
        }, ap);

        fprintf(stderr, "Parsimony score at the best node: %d\n", best_set_difference);
        // Show the paths to placement node(s)
        for (auto j: best_j_vec) {
            auto node = dfs[j];
            fprintf(stderr, "%s\t%d\t", clade_c_str, best_set_difference);
            fprintf(stderr, "%c\t%s\t", (node == best_node) ? '*' : '-', node->identifier.c_str());
            fprintf(stderr, "%s\n", path_to_node(T, node).c_str());
        }

        if (num_best > 1) {
            fprintf(stderr, "WARNING: found %zu possible assignments\n", num_best);
        }

        auto cm = clade_map.find(clade);
        if (cm != clade_map.end()) {
            std::vector<MAT::Node*>& clade_samples = cm->second;
            // Mutations for this clade were not specified with --clade-mutations, but instead by
            // representative samples; sort candidate clade_assignments by freq & overlap.
            clade_assignments.push_back(Clade_Assignments(clade, clade_samples.size(), std::move(clade_mutations)));

            float freq = 0;
            float overlap = 0;
            float best_freq = -1.0;

            std::stable_sort(best_j_vec.begin(), best_j_vec.end());

            for (auto j: best_j_vec) {
                // for each placement, this for loop traverses the ancestors of the
                // placement node all the way up the the root as long as the
                // frequency at which the clade sammples occur in each descendants
                // monotonically increases, and adds each of those ancestors to clade
                // assignments for consideration.
                for (auto anc: T.rsearch(dfs[j]->identifier, true)) {
                    get_freq_overlap(T, anc, clade_samples, &freq, &overlap);
                    if ((freq >= best_freq) && (overlap >= set_overlap)) {
                        clade_assignments[curr_idx].best_node_frequencies.push_back(Node_freq(dfs_idx[anc->identifier], freq, overlap, clip_sample_frequency));
                        best_freq = freq;
                    } else {
                        break;
                    }
                }
            }
            if (clade_assignments[curr_idx].best_node_frequencies.size() == 0) {
                fprintf(stderr, "WARNING: %s: no placement node or ancestor passed thresholds.\n", clade_c_str);
                for (auto j: best_j_vec) {
                    auto node = dfs[j];
                    get_freq_overlap(T, node, clade_samples, &freq, &overlap);
                    fprintf(stderr, "fail node\t%s\t%s\t%f\t%f\n", clade_c_str, node->identifier.c_str(), freq, overlap);
                }
            }

            // prioritize the assignments with the highest frequency and overlap
            std::stable_sort(clade_assignments[curr_idx].best_node_frequencies.begin(), clade_assignments[curr_idx].best_node_frequencies.end());

        } else {
            // Mutations were specified directly with --clade-mutations; choose the node with the
            // most leaf descendants.
            clade_assignments.push_back(Clade_Assignments(clade, 0, std::move(clade_mutations)));
            // However, if parsimony score > 0, don't assign the clade.  For example, the
            // SARS-CoV-2 public tree doesn't have any samples for some lineages; taking the
            // closest available node is the wrong thing to do, especially because it might belong
            // to some other lineage.
            if (best_set_difference == 0) {
                size_t best_node_j = 0;
                size_t max_num_leaves = 0;
                for (auto j: best_j_vec) {
                    MAT::Node *node = dfs[j];
                    size_t num_leaves = T.get_num_leaves(node);
                    if (num_leaves > max_num_leaves) {
                        best_node_j = j;
                        max_num_leaves = num_leaves;
                    }
                }
                clade_assignments[curr_idx].best_node_frequencies.push_back(Node_freq(best_node_j, 1, 1, clip_sample_frequency));
            } else {
                fprintf(stderr, "WARNING: skipping clade %s because the parsimony score is %d\n",
                        clade_c_str, best_set_difference);
            }
        }
        curr_idx++;
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    timer.Start();
    fprintf(stderr, "Sorting clades by the number of best nodes \n");
    std::sort(clade_assignments.begin(), clade_assignments.end());
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Now assigning clades to nodes \n");
    for (auto c: clade_assignments) {
        bool assigned = false;
        Node_freq *assigned_n = NULL;
        std::string already_assigned = "";
        for (auto n: c.best_node_frequencies) {
            auto j = n.best_j;
            if (dfs[j]->clade_annotations[num_annotations-1] == "") {
                fprintf(stderr, "\nAssigning %s to node %s\n", c.clade_name.c_str(), dfs[j]->identifier.c_str());
                if (c.clade_size > 0) {
                    fprintf(stderr, "%f fraction of %zu clade %s samples are descendants of the assigned node %s\n", n.overlap, c.clade_size, c.clade_name.c_str(), dfs[j]->identifier.c_str());
                }
                dfs[j]->clade_annotations[num_annotations-1] = c.clade_name;
                assigned = true;
                assigned_n = &n;

                break;
            } else {
                std::string node_name = dfs[j]->identifier;
                std::string assigned_clade = dfs[j]->clade_annotations[num_annotations-1];
                fprintf(stderr, "\nNode %s already assigned to %s, cannot assign to %s.\n",
                        node_name.c_str(), assigned_clade.c_str(), c.clade_name.c_str());
                if (already_assigned.size() > 0) {
                    already_assigned += ", ";
                }
                already_assigned += (node_name + ":" + assigned_clade);
            }
        }

        if (!assigned) {
            if (c.best_node_frequencies.size() > 0) {
                fprintf(stderr, "\nWARNING: Could not assign a node to clade %s with %zu samples since all possible nodes were already assigned some other clade!\n", c.clade_name.c_str(), c.clade_size);
            } else {
                fprintf(stderr, "\nWARNING: Could not assign a node to clade %s with %zu samples since placement node(s) did not overlap with enough clade samples!\n", c.clade_name.c_str(), c.clade_size);
            }
        }
        if (details_file) {
            fprintf(details_file, "%s\t", c.clade_name.c_str());
            write_mutations(details_file, c.mutations, true, false);
            fputc('\t', details_file);
            write_mutations(details_file, c.mutations, false, true);
            fputc('\t', details_file);
            if (c.best_node_frequencies.size() == 0) {
                fprintf(details_file, "n/a");
            } else {
                for (size_t i = 0;  i < c.best_node_frequencies.size();  i++) {
                    if (i > 0) {
                        fprintf(details_file, ", ");
                    }
                    Node_freq bnf = c.best_node_frequencies[i];
                    MAT::Node *node = dfs[bnf.best_j];
                    fprintf(details_file, "%s:%f:%f", node->identifier.c_str(), bnf.freq, bnf.overlap);
                }
            }
            if (already_assigned == "") {
                fprintf(details_file, "\tn/a");
            } else {
                fprintf(details_file, "\t%s", already_assigned.c_str());
            }
            fprintf(details_file, "\t%f", assigned ? assigned_n->overlap : 0.0);
            fprintf(details_file, "\t%lu", c.clade_size);
            std::string path_string = "n/a";
            if (assigned) {
                path_string = path_to_node(T, dfs[assigned_n->best_j]);
            } else if (c.best_node_frequencies.size() > 0) {
                path_string = path_to_node(T, dfs[c.best_node_frequencies[0].best_j]);
            }
            fprintf(details_file, "\t%s\n", path_string.c_str());
        }
    }

    if (mutations_file) {
        fclose(mutations_file);
    }
    if (details_file) {
        fclose(details_file);
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

bool node_has_muts(MAT::Node *node, std::vector<std::string>& muts) {
    // Return true if node's mutations are the same as muts, where each string in muts is
    // of the form <par_nuc><position><mut_nuc>.  Mutations can appear in different orders.
    // The length of muts is expected to be very small (typically 1, usually < 5).
    if (node->mutations.size() != muts.size()) {
        return false;
    }
    std::vector<std::string> node_muts(node->mutations.size());
    for (MAT::Mutation node_mut: node->mutations) {
        bool got_it = false;
        std::string node_mut_string = node_mut.get_string();
        for (std::string mut: muts) {
            if (mut == node_mut_string) {
                got_it = true;
                break;
            }
        }
        if (! got_it) {
            return false;
        }
    }
    return true;
}

void assignLineagesFromPaths(MAT::Tree& T, const std::string& clade_paths_filename, std::unordered_set<std::string>& clades_already_assigned) {
    // Read in a tsv file of clade names and mutation paths like first and third columns of the
    // output of matUtils extract --clade-paths e.g. "B\t" (B, no mutations) or
    // "A\t> T28144C > C8782T" or B.52\tC8090T > A1904G,C10789T".
    // Follow the path of mutations from root and assign clade to terminal node of path;
    // warn and move on to next clade & path if no child has the specified mutations.
    // Return a set of clades whose paths were found.
    auto dfs = T.depth_first_expansion();
    size_t num_annotations = T.get_num_annotations();

    std::ifstream infile(clade_paths_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the clade paths file: %s!\n", clade_paths_filename.c_str());
        exit(1);
    }
    std::string line;

    fprintf(stderr, "Reading clade paths file and making assignments.\n");
    int found_count = 0, fail_count = 0;
    timer.Start();
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, '\t', words);
        // Empty second word (for root node, no mutations) is ignored by string_split; add it back.
        if (words.size() == 1 && line.back() == '\t') {
            words.push_back("");
        }
        if (words.size() != 2) {
            fprintf(stderr, "ERROR: Incorrect format for clade paths file: %s!  Expected 2 tab-separated words, got %ld (%s)\n", clade_paths_filename.c_str(), words.size(), line.c_str());
            exit(1);
        }
        std::string clade = words[0];
        if (clades_already_assigned.find(clade) != clades_already_assigned.end()) {
            fprintf(stderr, "Clade %s has already been assigned, ignoring path.\n", clade.c_str());
            continue;
        }
        std::string path = words[1];
        // Parse path from words[1] and find terminal node of path, or exit with error.
        std::vector<std::string> path_words;
        MAT::string_split(path, path_words);
        MAT::Node* node = T.root;
        bool failed = false;
        for (std::string path_el: path_words) {
            // Ignore empty string or ">"
            if (path_el == "" || path_el == ">") {
                continue;
            }
            // Might be comma-sep list of mutations
            std::vector<std::string> muts;
            MAT::string_split(path_el, ',', muts);
            MAT::Node *child_with_muts = NULL;
            for (MAT::Node *kid: node->children) {
                if (node_has_muts(kid, muts)) {
                    child_with_muts = kid;
                    break;
                }
            }
            if (child_with_muts == NULL) {
                fprintf(stderr, "WARNING: path for %s not found: no child of node %s has mutations %s (path %s)\n",
                        clade.c_str(), node->identifier.c_str(), path_el.c_str(), line.c_str());
                failed = true;
                fail_count++;
                break;
            }
            node = child_with_muts;
        }
        if (!failed) {
            // assign clade to node
            if (node->clade_annotations[num_annotations-1] != "") {
                fprintf(stderr, "WARNING: Assigning clade %s to node %s for path %s failed as the node is already assigned to clade %s!\n",
                        clade.c_str(), node->identifier.c_str(), path.c_str(), node->clade_annotations[num_annotations-1].c_str());
            } else {
                node->clade_annotations[num_annotations-1] = clade;
            }
            clades_already_assigned.insert(clade);
            found_count++;
        }
    }
    infile.close();
    fprintf(stderr, "\nAnnotated %d clades; failed to find paths for %d clades\n", found_count, fail_count);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
