#include "uncertainty.hpp"

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
        // for (auto i = distances.begin(); i != distances.end(); ++i) {
        //     std::cerr << *i << ',';
        // }
        // //now find the biggest sum of shared values for this specific common node
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

void findEPPs (MAT::Tree Tobj, std::string sample_file, std::string fepps, std::string fneigh) {
    timer.Start();
    std::ofstream eppfile;
    std::ofstream neighfile;
    if (fepps != ""){
        eppfile.open(fepps);
        //add column names
        //auspice doesn't care what's actually in column 1 name, but it does for column 2
        eppfile << "sample\tepps\n";
    }
    if (fneigh != ""){
        neighfile.open(fneigh);
        neighfile << "sample\tneighborhood_size\n";
    }
    
    MAT::Tree* T = &Tobj; //mapper wants a pointer.

    std::vector<MAT::Node*> fdfs;

    if (sample_file != "") { //read in the samples files and get the nodes corresponding to each sample.
        std::ifstream infile(sample_file);
        if (!infile) {
            fprintf(stderr, "ERROR: could not open the indicated sample file\n");
            exit(1);
        }
        std::string psample;
        std::string sample;
        while (std::getline(infile, psample)) {
            //adding code to handle carriage returns in the sample text file
            //in case a text file edited on Windows is uploaded or passed into matUtils
            if (psample.size() && psample[psample.size()-1] == '\r') {
                sample = psample.substr(0,psample.size()-1);
            } else {
                sample = psample;
            }

            if (T->get_node(sample) == NULL) {
                fprintf(stderr, "ERROR: Sample missing in input MAT!\n");
                std::cerr << sample;
                std::cerr << std::endl;
                exit(1);
            }
            fdfs.emplace_back(T->get_node(sample));
        }
    } else {
        //if unspecified, it does it for all samples in the MAT
        //which can take several hours on the full reference tree
        //this is useful default behavior for future subtree MATs though
        fprintf(stderr, "No sample file specified; calculating metrics for all samples in the tree\n");
        fdfs = Tobj.depth_first_expansion();
    }
    fprintf(stderr, "Processing %ld samples\n", fdfs.size());

    for (size_t s=0; s<fdfs.size(); s++){ //this loop is not a parallel for because its going to contain a parallel for.
        //get the node object.
        auto node = fdfs[s];
        if (node->is_leaf()) { 
            assert (node->children.size() == 0); //leaf nodes need to be leaves.
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
                
                if (fepps != "") {
                    eppfile << node->identifier << "\t" << num_best << "\n";
                }

                if (fneigh != "") { 
                    if (num_best > 1) { //only worth calculating if there's more than one best placement. Otherwise its just 0.
                        std::vector<MAT::Node*> best_placements;
                        //for every index in best_j_vec, find the corresponding node from dfs
                        for (size_t z=0; z<best_j_vec.size(); z++) {
                            auto nobj = dfs[best_j_vec[z]];
                            best_placements.emplace_back(nobj);
                        }
                        size_t neighborhood_size = get_neighborhood_size(best_placements, T);
                        neighfile << node->identifier << "\t" << neighborhood_size << "\n";
                    } else {
                        //one best placement, total distance is 0
                        neighfile << node->identifier << "\t0\n";
                    }
                }
            } 
            else {
                //save default values 
                if (fepps != "") {
                    eppfile << node->identifier << "\t1\n";
                }
                if (fneigh != "") {
                    neighfile << node->identifier << "\t0\n";
                }

            }
        }
    }
    if (fepps != ""){
        eppfile.close();
    }
    if (fneigh != ""){
        neighfile.close();
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

po::variables_map parse_uncertainty_command(po::parsed_options parsed) {    

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description ann_desc("annotate options");
    ann_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("samples,s", po::value<std::string>()->default_value(""),
        "Path to a simple text file of sample names to calculate uncertainty for.")
        ("get-parsimony,g", po::bool_switch(),
        "Calculate and print the total tree parsimony score.")
        ("find-epps,e", po::value<std::string>()->default_value(""),
        "Name for an Auspice-compatible tsv file output of equally parsimonious placements for input samples.")
        ("find-neighborhood,n", po::value<std::string>()->default_value(""),
        "Name for an Auspice-compatible tsv file output of neighborhood score values for the equally parsimonious placements for input samples.")
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

void uncertainty_main(po::parsed_options parsed) {
    //the uncertainty command, for now, produces Auspice-compatible TSV files for visualization
    //in the future I would like it to be able to produce a metadata-annotated JSON file if indicated
    po::variables_map vm = parse_uncertainty_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string sample_file = vm["samples"].as<std::string>();
    bool get_parsimony = vm["get-parsimony"].as<bool>();
    std::string fepps = vm["find-epps"].as<std::string>();
    std::string fneigh = vm["find-neighborhood"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    if (fepps != "" || fneigh != "") {
        fprintf(stderr, "Calculating placement uncertainty\n");
        findEPPs(T, sample_file, fepps, fneigh);
    }
    if (get_parsimony){
        fprintf(stderr, "Total Tree Parsimony %ld", T.get_parsimony_score());
    }
}
