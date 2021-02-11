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
        ("lineage-names,l", po::value<std::string>()->required(),
         "Path to a file containing lineage asssignments of samples. Use to locate and annotate clade root nodes")
        ("allele-frequency,f", po::value<float>()->default_value(0.8),
         "Minimum allele frequency in input samples for finding the best clade root. Used only with -l")
        ("set-overlap,s", po::value<float>()->default_value(0.9),
        "Minimum fraction of the lineage samples that should be desecendants of the assigned clade root")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(ann_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << ann_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void annotate_main(po::parsed_options parsed) {
    //the annotate subcommand calculates and saves information about the tree, returning a protobuf file that is larger than the input
    po::variables_map vm = parse_annotate_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string lineage_filename = vm["lineage-names"].as<std::string>();
    float allele_frequency = vm["allele-frequency"].as<float>();
    float set_overlap = vm["set-overlap"].as<float>();
    //    bool get_parsimony = vm["get-parsimony"].as<bool>();
    //    bool fepps = vm["find-epps"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    // If the argument to calculate equally parsimonious placements was used, perform this operation
    //    if (fepps) {
    //        fprintf(stderr, "Calculating EPPs\n");
    //        T = findEPPs(T);
    //    }
    //    if (get_parsimony){
    //        //fprintf(stderr, "Calculating Total Parsimony\n");
    //        //T.total_parsimony = T.get_parsimony_score();
    //    }
    if (lineage_filename != "") {
        fprintf(stderr, "Annotating Lineage Root Nodes\n");
        assignLineages(T, lineage_filename, allele_frequency, set_overlap);
    }
    else {
        fprintf(stderr, "ERROR: must specifiy lineage-names!\n");
        exit(1);
    }

    //condense_leaves() expects some samples to ignore. We don't have any such samples
    //this would be space to add an additional argument containing samples to not recondense
    //for now, just recondense everything
    T.condense_leaves(std::vector<std::string>());

    // Store final MAT to output file
    if (output_mat_filename != "") {
        fprintf(stderr, "Saving Final Tree\n");
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
    }
}

void assignLineages (MAT::Tree& T, const std::string& lineage_filename, float min_freq, float set_overlap) {
    static tbb::affinity_partitioner ap;
    
    fprintf(stderr, "Copying tree with uncondensed leaves.\n"); 
    timer.Start();
    auto uncondensed_T = MAT::get_tree_copy(T);
    uncondensed_T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    
    std::map<std::string, std::vector<MAT::Node*>> lineage_map;
    std::ifstream infile(lineage_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the lineage asssignment file: %s!\n", lineage_filename.c_str());
        exit(1);
    }    
    std::string line;

    //Read the lineage assignment file line by line and fill up map values
    fprintf(stderr, "Reading lineage assignment file.\n"); 
    timer.Start();
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if ((words.size() > 2) || (words.size() == 1)) {
            fprintf(stderr, "ERROR: Incorrect format for lineage asssignment file: %s!\n", lineage_filename.c_str());
            exit(1);
        }
        else if (words.size() == 2) {
            auto n = uncondensed_T.get_node(words[1]);
            if (n != NULL) {
                if (lineage_map.find(words[0]) == lineage_map.end()) {
                    lineage_map[words[0]] = std::vector<MAT::Node*>();;
                }
                lineage_map[words[0]].emplace_back(n);
            }
            else {
                fprintf(stderr, "WARNING: Sample %s not found in input MAT!\n", words[1].c_str());
            }
        }
    }
    infile.close();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());


    for (auto it: lineage_map) {
        fprintf(stderr, "Finding best node for lineage %s.\n", it.first.c_str()); 
        timer.Start();
        
        std::map<std::string, int> mutation_counts;

        tbb::mutex tbb_lock;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, it.second.size()),
        [&](const tbb::blocked_range<size_t> r) {
            for (size_t i=r.begin(); i<r.end(); ++i){
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
                        }
                        else {
                            mutation_counts[mut_string] += 1;
                        }
                        tbb_lock.unlock();
                    }
                }
            }
        }, ap);

        std::vector<MAT::Mutation> lineage_mutations;
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
                lineage_mutations.emplace_back(m);
            }
        }
        // Mutations need to be sorted by position before placement
        std::sort(lineage_mutations.begin(), lineage_mutations.end());
        fprintf(stderr, "Number of mutations above the specified frequency in this clade: %zu \n", lineage_mutations.size());
        
        auto dfs = T.depth_first_expansion();
        size_t total_nodes = dfs.size();
        
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

        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes),
                [&](tbb::blocked_range<size_t> r) {
                for (size_t k=r.begin(); k<r.end(); ++k){
                   mapper2_input inp;
                   inp.T = &T;
                   inp.node = dfs[k];
                   inp.missing_sample_mutations = &lineage_mutations;
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
        
        struct Node_freq {
            size_t best_j;
            float freq;
            float overlap;
            Node_freq (size_t a, float b, float c) 
            {
                best_j = a;
                freq = b;
                overlap = c;
            }
            // To sort with highest frequency first
            inline bool operator< (const Node_freq& n) const {
                return ((this->freq > n.freq) || ((this->freq == n.freq) && (this->best_j < n.best_j)));
            }
        };

        std::vector<Node_freq> best_node_frequencies;
        if (num_best > 1) {
            fprintf(stderr, "WARNING: found %zu possible assignments\n", num_best);
        }
        
        float freq = 0;
        float overlap = 0;
            
        for (auto j: best_j_vec) {
            size_t num_desc = 0;

            tbb::parallel_for(tbb::blocked_range<size_t>(0, it.second.size()),
                    [&](const tbb::blocked_range<size_t> r) {
                    for (size_t i=r.begin(); i<r.end(); ++i){
                    if (T.is_ancestor(dfs[j]->identifier, it.second[i]->identifier)) {
                    __sync_fetch_and_add(&num_desc, 1.0);
                    }
                    }
                    }, ap);

            overlap = (static_cast<float>(num_desc) / it.second.size());
            if (overlap >= set_overlap) {
                freq = (static_cast<float>(num_desc) / T.get_num_leaves(dfs[j])); 
                best_node_frequencies.push_back(Node_freq(j, freq, overlap));
            }
        }
        
        std::sort(best_node_frequencies.begin(), best_node_frequencies.end());

        bool assigned = false;
        for (auto n: best_node_frequencies) {
            auto j = n.best_j;
            if (dfs[j]->clade == "") {
                fprintf(stderr, "Assigning %s to node %s\n", it.first.c_str(), dfs[j]->identifier.c_str());
                fprintf(stderr, "%f fraction of the lineage samples are descendants of the assigned node %s\n", n.overlap, dfs[j]->identifier.c_str());
                dfs[j]->clade = it.first;
                assigned = true;
                break;
            }
        }

        if (!assigned) {
            if (best_node_frequencies.size() > 0) {
                fprintf(stderr, "WARNING: Could not assign a node to clade %s since all possible nodes were already assigned some other clade!\n", it.first.c_str());
            }
            else {
                fprintf(stderr, "WARNING: Could not assign a node to clade %s since placement node(s) did not overlap with enough lineage samples!\n", it.first.c_str());
            }
        }

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}

//The below is commented out because it was accidental that I committed it to my master branch when updating the PR; its not ready and I'm still debugging it on a dedicated branch 
//do feel free to read through it and comment on the method/concept if you want though
//Basically, I'm using the below code to calculate the maximum pairwise distance between the most recent common ancestor of an arbitrary group of nodes
//this is an uncertainty metric related to equally parsimonious placements- big values mean bad samples
//a sample that has several EPPs that are very nearby on the tree will have a small neighborhood size value
//while a sample that only has two or three EPPs but they are on completely different parts of the tree will have a substantially larger one.

std::vector<MAT::Node*> get_common_nodes (std::vector<std::vector<MAT::Node*>> nodepaths) {
    //to identify common nodes, perform pairwise set intersections repeatedly for all path node vectors
    std::vector<MAT::Node*> common_nodes = nodepaths[0];
    std::sort(common_nodes.begin(), common_nodes.end()); //needs to be sorted for intersection. These are actually vectors of node POINTERS, so this should be fine.
    for (size_t s=1; s<nodepaths.size(); s++) {
        std::vector<MAT::Node*> nextint;
        std::vector<MAT::Node*> next = nodepaths[s];
        std::sort(next.begin(), next.end()); //sort each path
        std::set_intersection(common_nodes.begin(), common_nodes.end(), next.begin(), next.end(), std::back_inserter(nextint)); //intersect the values
        common_nodes = nextint; //store the intersected vector and move to the next node path vector
    }
    return common_nodes;
}

std::vector<float> get_all_distances(MAT::Node* target, std::vector<std::vector<MAT::Node*>> paths){
    std::vector<float> distvs;
    for (size_t p=0; p<paths.size(); p++) {
        //for this path to the target common ancestor, count up distances from the start until it is reached
        float tdist = 0;
        assert (paths[p].size() > 0);
        for (size_t i=0;i<paths[p].size();i++) {
            if (paths[p][i]->identifier == target->identifier) {
                break; //stop iterating when its reached this common ancestor (remember, path is sorted nearest to root)
            }
            tdist += paths[p][i]->branch_length;
            
        //then record tdist in distvs
        distvs.emplace_back(tdist);
        }
    }
    return distvs;
}

size_t get_neighborhood_size(std::vector<MAT::Node*> nodes, MAT::Tree* T) {
    //first step for this is collecting the full paths back to the root for all nodes
    assert (nodes.size() > 1); //doesn't make sense if there's only one best placement.
    std::vector<std::vector<MAT::Node*>> parentvecs;
    for (size_t s=0; s<nodes.size(); s++) {
        if (!nodes[s]->is_root()){ //if one of the epps sites is directly off the root, then this doesn't make much sense
            std::vector<MAT::Node*> npath;
            npath.emplace_back(nodes[s]); //include the node itself on the path so branch length is correctly accessed
            for (auto& it: T->rsearch(nodes[s]->identifier)) {
                npath.emplace_back(it);
            }
            parentvecs.emplace_back(npath);
        } else { //instead, just construct a path of length 1 that contains the root node only
            std::vector<MAT::Node*> npath;
            npath.emplace_back(nodes[s]);
            parentvecs.emplace_back(npath);
        }
    }
    //then we need to identify all common elements to all node vectors
    std::vector<MAT::Node*> common_nodes = get_common_nodes(parentvecs);
    assert (common_nodes.size() > 0); //bare minimum this will always include the root. therefore it is always > 0
    //then for all common nodes, we need to calculate the largest sum of paired distances for all samples to that specific common ancestor
    //the smallest of these largest sums is the best neighborhood size value
    size_t best_size = T->get_parsimony_score(); //bigger than the biggest maximum limit on neighborhood size. basically a big number
    for (size_t s=0; s<common_nodes.size(); s++) {
        //get the set of distances between each placement to this common ancestor with the path vectors
        std::vector<float> distances = get_all_distances(common_nodes[s], parentvecs);
        //now find the biggest sum of shared values for this specific common node
        float widest = 0.0;
        for (size_t i=0; i<distances.size(); i++) {
            for (size_t j=0; j<distances.size(); j++) {
                if (i!=j){
                    float spairdist = distances[i] + distances[j];
                    if (spairdist > widest){
                        widest = spairdist;
                    }
                }
            }
        }
        //after that oddness, I now have a value which is the longest path going between any two nodes in the placement set
        //which goes through this specific common ancestor
        //admittedly this is probably not the fastest way to do this, I should be able to eliminate common ancestors which are directly ancestral to common ancestors between the same complete set without counting them up
        //but I'm focusing on results first here
        //anyways, assign this longest pair path value to best_size if its smaller than any we've seen for other common ancestors
        size_t size_widest = static_cast<size_t>(widest);
        if (size_widest < best_size){
            best_size = size_widest;
        }
    }
    //at the end, we should be left with the proper neighborhood size value. 
    return best_size;
}

MAT::Tree findEPPs (MAT::Tree Tobj) {
    TIMEIT()
    //all comments on function JDM
    //first, need to iterate through all leaf nodes, including condensed nodes and single leaves
    //internal nodes have metadata objects but those are just gonna be default values 0 for epps for now which should indicate high confidence anyways
    //the simplest way to do this is to iterate through all nodes and check whether each one is a leaf

    MAT::Tree* T = &Tobj; //mapper wants a pointer.
    std::vector<MAT::Node*> fdfs = Tobj.depth_first_expansion();

    for (size_t s=0; s<fdfs.size(); s++){ //this loop is not a parallel for because its going to contain a parallel for.
        //get the node object.
        auto node = fdfs[s];
        if (node->is_leaf()) { 
            assert (node->children.size() == 0); //leaf nodes need to be leaves.
            //fprintf(stderr, "Node- ID %s ", node->identifier.c_str()); //write the identifier to stderr
            //retrieve the full set of mutations associated with this Node object from root to it
            //to do this, get the full set of ancestral nodes and their mutations
            //code copied from the usher mapper.
            std::vector<int> anc_positions; //tracking positions is required to account for backmutation/overwriting along the path
            std::vector<MAT::Mutation> ancestral_mutations;
            //first load in the current mutations
            for (auto m: node->mutations){ //I don't fully understand this code block from Angie's vcf generator; likely that at least some of it is unnecessary. Review?
                if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                    ancestral_mutations.emplace_back(m);
                    if (!m.is_masked()) {
                        anc_positions.emplace_back(m.position);
                    }            
                }
            }
            //then load in ancestral mutations
            for (auto n: Tobj.rsearch(node->identifier)) {
                for (auto m: n->mutations) {
                    if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                        ancestral_mutations.emplace_back(m);
                        if (!m.is_masked()) {
                            anc_positions.emplace_back(m.position);
                        }
                    }
                }
            }
            //fprintf(stderr, "Mutations %ld ", ancestral_mutations.size());
            //if there are any mutations in the set.
            if (ancestral_mutations.size()>0) {
                //the ancestral_mutations vector, plus the mutations assigned to this specific node, constitute the "missing_sample" equivalents for calling the mapper
                //COPIED FROM SOME PART OF USHER 
                auto dfs = T->depth_first_expansion();
                size_t total_nodes = dfs.size();

                // Stores the excess mutations to place the sample at each
                // node of the tree in DFS order. When placement is as a
                // child, it only contains parsimony-increasing mutations in
                // the sample. When placement is as a sibling, it contains 
                // parsimony-increasing mutations as well as the mutations
                // on the placed node in common with the new sample. Note
                // guaranteed to be corrrect only for optimal nodes since
                // the mapper can terminate the search early for non-optimal
                // nodes
                std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
                // Stores the imputed mutations for ambiguous bases in the
                // sampled in order to place the sample at each node of the 
                // tree in DFS order. Again, guaranteed to be corrrect only 
                // for pasrimony-optimal nodes 
                std::vector<std::vector<MAT::Mutation>> node_imputed_mutations(total_nodes);

                // Stores the parsimony score to place the sample at each
                // node of the tree in DFS order.
                std::vector<int> node_set_difference;
                size_t best_node_num_leaves = 0;
                // The maximum number of mutations is bound by the number
                // of mutations in the missing sample (place at root)
                //int best_set_difference = 1e9;
                // TODO: currently number of root mutations is also added to
                // this value since it forces placement as child but this
                // could be changed later 
                int best_set_difference = ancestral_mutations.size() + T->root->mutations.size() + 1;

                size_t best_j = 0;
                size_t num_best = 1;
                bool best_node_has_unique = false;
                MAT::Node* best_node = T->root;
                std::vector<bool> node_has_unique(total_nodes, false);
                std::vector<size_t> best_j_vec;
                best_j_vec.emplace_back(0);

                // Parallel for loop to search for most parsimonious
                // placements. Real action happens within mapper2_body
                auto grain_size = 400; 
                tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes, grain_size),
                        [&](tbb::blocked_range<size_t> r) {
                        for (size_t k=r.begin(); k<r.end(); ++k){
                            if (dfs[k]->identifier != node->identifier) { //do not allow self-mapping (e.g. can't remap leaf as child of itself)
                                mapper2_input inp;
                                inp.T = T;
                                inp.node = dfs[k];
                                inp.missing_sample_mutations = &ancestral_mutations;
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
                        }       
                        }); 
                //node->epps = num_best; //save the value.
                //additional metadata attribute assigning related to remapping would go here.

                //additional metadata value (not even starting to assign it to an attribute yet) is maximum pairwise distance between equally parsimonious placement cluster members
                //commented out code that was accidentally committed to master
                //to find this, first we need the nodes.
                if (num_best > 1){ //only worth calculating if there's more than one best placement. Otherwise its just 0.
                    std::vector<MAT::Node*> best_placements;
                    //for every index in best_j_vec, find the corresponding node from dfs
                    for (size_t z=0; z<best_j_vec.size(); z++) {
                        auto nobj = dfs[best_j_vec[z]];
                        best_placements.emplace_back(nobj);
                    }
                    //size_t neighborhood_size = get_neighborhood_size(best_placements, T);
                    //node->neighborhood_size = neighborhood_size;
                    //fprintf(stderr, "Neighborhood Size: %ld\n", neighborhood_size);
                } 
                //else {
                    //node->neighborhood_size = 0;
                    //fprintf(stderr, "Neighborhood Size: 0\n");
                //}

            } 
            //else {
                //node->epps = 1;
                //node->neighborhood_size = 0;
                //fprintf(stderr, "Neighborhood Size: 0\n");
                //no mutations for this sample compared to the reference. This means it's leaf off the root/identical to the reference
                //there's just one place for that, ofc.
            //}
        }
    }
    return Tobj; //return the actual object.
}

