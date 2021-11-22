#include "uncertainty.hpp"

std::vector<MAT::Node*> get_common_nodes (std::vector<std::vector<MAT::Node*>> nodepaths) {
    //to identify common nodes, perform pairwise set intersections repeatedly for all path node vectors
    std::vector<MAT::Node*> common_nodes = nodepaths[0];
    //needs to be sorted for intersection. These are vectors of pointers, so there's not a problem with how to order them.
    std::sort(common_nodes.begin(), common_nodes.end());
    for (size_t s=1; s<nodepaths.size(); s++) {
        std::vector<MAT::Node*> nextint;
        std::vector<MAT::Node*> next = nodepaths[s];
        std::sort(next.begin(), next.end());
        std::set_intersection(common_nodes.begin(), common_nodes.end(), next.begin(), next.end(), std::back_inserter(nextint));
        //store the intersected vector and move to the next node path vector
        common_nodes = nextint;
    }
    //the final vector is the set of nodes which existed in every path vector
    return common_nodes;
}

std::vector<float> get_all_distances(MAT::Node* target, std::vector<std::vector<MAT::Node*>> paths) {
    std::vector<float> distvs;
    for (size_t p=0; p<paths.size(); p++) {
        //for this path to the target common ancestor, count up distances from the start until it is reached
        float tdist = 0;
        assert (paths[p].size() > 0);
        for (size_t i=0; i<paths[p].size(); i++) {
            if (paths[p][i]->identifier == target->identifier) {
                //stop iterating when its reached the indicated common ancestor (remember, path is sorted nearest to root)
                break;
            }
            tdist += paths[p][i]->mutations.size();

            //then record the total distance in distvs
            distvs.emplace_back(tdist);
        }
    }
    return distvs;
}

size_t get_neighborhood_size(std::vector<MAT::Node*> nodes, MAT::Tree* T) {
    /*
    The basic concept behind neighborhood size is that it is the longest direct path
    (direct meaning without passing over the same connection twice)
    between any two nodes in the input set. This is parallel to the concept of the
    complete-linkage cluster density value.

    Obtaining this value requires two basic steps. First, the most recent common ancestor,
    which the longest traversible direct path must pass through, must be identified.
    To do this, the set of all common ancestors for all nodes in the input set must be identified.
    Then the total distances for each sample to each common ancestor must be tabulated; the
    most recent common ancestor is the one with the shortest total distance to all nodes in the input.

    After the most recent common ancestor is identified, the distances from each input node
    to that ancestor are calculated. The longest path is the sum of the two largest distance values
    from an input node to the common ancestor.

    This metric is useful because it is indicative of the distribution of placements across the tree
    when ran on the set of nodes for which a sample could be equally parsimoniously placed.
    For example, let's say we have two samples with five equally parsimonious placements.
    The first sample has all five equally parsimonious placements in a small clade with a very recent common
    ancestor. The second has two of these placements on a very different parts of the tree from the other three,
    with the common ancestor back at the root. If you only look at the raw number of equally parsimonious placements,
    they look equally trustworthy, but its likely that you could believe that the first sample is properly placed
    somewhere in that small clade, while the second sample should be removed or not trusted at all.
    The first of these two samples will have a small neighborhood size value
    while the second will have a large neighborhood size value. This metric thus complements the
    number of equally parsimonious placements when evaluating sample placement quality.
    */

    //first step for this is collecting the full paths back to the root for all nodes
    assert (nodes.size() > 1);
    std::vector<std::vector<MAT::Node*>> parentvecs;
    for (size_t s=0; s<nodes.size(); s++) {
        if (!nodes[s]->is_root()) {
            std::vector<MAT::Node*> npath;
            npath.emplace_back(nodes[s]); //include the node itself on the path so branch length is correctly accessed
            for (auto& it: T->rsearch(nodes[s]->identifier)) {
                npath.emplace_back(it);
            }
            parentvecs.emplace_back(npath);
        } else {
            //if one of the epps sites is directly off the root, then this doesn't make much sense
            //instead, just construct a path of length 1 that contains the root node only
            std::vector<MAT::Node*> npath;
            npath.emplace_back(nodes[s]);
            parentvecs.emplace_back(npath);
        }
    }
    //then we need to identify all common elements to all node vectors
    std::vector<MAT::Node*> common_nodes = get_common_nodes(parentvecs);
    //bare minimum this will always include the root. therefore it is always > 0
    assert (common_nodes.size() > 0);
    //then for all common nodes, we need to calculate the largest sum of paired distances for all samples to that specific common ancestor
    //the smallest of these largest sums is the best neighborhood size value

    //the parsimony score is bigger than the biggest maximum limit on neighborhood size.
    //This is basically just a big number to compare to
    //TODO: There's definitely a better way to bound this that doesn't involve as much adding.
    size_t best_size = T->get_parsimony_score();

    for (size_t s=0; s<common_nodes.size(); s++) {
        //get the set of distances between each placement to this common ancestor with the path vectors
        std::vector<float> distances = get_all_distances(common_nodes[s], parentvecs);
        float widest = 0.0;
        for (size_t i=0; i<distances.size(); i++) {
            for (size_t j=0; j<distances.size(); j++) {
                if (i!=j) {
                    float spairdist = distances[i] + distances[j];
                    if (spairdist > widest) {
                        widest = spairdist;
                    }
                }
            }
        }
        //I now have a value which is the longest path going between any two nodes in the placement set
        //which goes through this specific common ancestor
        //admittedly this is probably not the fastest way to do this
        //I should be able to eliminate common ancestors which are directly ancestral
        //to common ancestors between the same complete set without counting them up
        //but I'm focusing on results first here
        //assign this longest pair path value to best_size if its smaller than any we've seen for other common ancestors
        size_t size_widest = static_cast<size_t>(widest);
        if (size_widest < best_size) {
            best_size = size_widest;
        }
    }
    //at the end, we should be left with the proper neighborhood size value.
    return best_size;
}

std::vector<MAT::Node*> findEPPs (MAT::Tree* T, MAT::Node* node, size_t* nbest, size_t* nsize) {
    //calculating neighborhood size is optional.

    //retrieve the full set of mutations associated with this Node object from root to it
    //to do this, get the full set of ancestral nodes and their mutations

    //NOTE (TODO): The below code was copied from the Usher main.cpp
    //if that code is ever refactored or updated in a way that significantly
    //affects efficiency or accuracy outside of the usher_mapper.cpp,
    //this code needs to be manually updated

    //tracking positions is required to account for backmutation/overwriting along the path
    std::vector<int> anc_positions;
    std::vector<MAT::Mutation> ancestral_mutations;
    std::vector<MAT::Node*> best_placements;
    //first load in the current mutations
    for (auto m: node->mutations) {
        if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
            ancestral_mutations.emplace_back(m);
            if (!m.is_masked()) {
                anc_positions.emplace_back(m.position);
            }
        }
    }
    //then load in ancestral mutations
    for (auto n: T->rsearch(node->identifier)) {
        for (auto m: n->mutations) {
            if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                ancestral_mutations.emplace_back(m);
                if (!m.is_masked()) {
                    anc_positions.emplace_back(m.position);
                }
            }
        }
    }
    if (ancestral_mutations.size()>0) {
        //the ancestral_mutations vector, plus the mutations assigned to this specific node, constitute the "missing_sample" equivalents for calling the mapper
        //this is where code copied from usher_main.cpp begins
        //comments included.
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
            for (size_t k=r.begin(); k<r.end(); ++k) {
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
        //assign the result to the input pointer for epps
        *nbest = num_best;
        //if the num_best is big enough and if the bool is set, get the neighborhood size value and assign it
        if (num_best > 1) {
            //for every index in best_j_vec, find the corresponding node from dfs
            for (size_t z=0; z<best_j_vec.size(); z++) {
                auto nobj = dfs[best_j_vec[z]];
                best_placements.emplace_back(nobj);
            }
            size_t neighborhood_size = get_neighborhood_size(best_placements, T);
            *nsize = neighborhood_size;
        } else {
            //one best placement, total distance is 0
            *nsize = 0;
            //record the original parent of this node as the single best placement.
            best_placements.emplace_back(node->parent);
        }
    }
    return best_placements;
}

void findEPPs_wrapper (MAT::Tree Tobj, std::string sample_file, std::string fepps, std::string flocs) {
    /*
    The number of equally parsimonious placements (EPPs) is a placement uncertainty metric that
    indicates when a sample is ambiguous and could have been produced by more than one path
    from the root (e.g. multiple possible origins or explanations for the sample).
    About 85% of samples historically have a single unique placement (EPPs=1).

    This value is calculated as part of the standard Usher sample placement code. In order to get
    these values for a specific set of already-placed samples, we identify where each indicated sample
    is placed on the tree, construct a virtual sample containing its associated set of mutations
    from the root, and run Usher's placement step, specifically blocking mapping of
    the extracted sample to itself on the tree. The resulting values are recorded, and the set
    of EPP nodes can be used to calculate the complementary metric neighborhood size.
    */
    timer.Start();
    std::ofstream eppfile;
    if (fepps != "") {
        eppfile.open(fepps);
        eppfile << "sample\tequally_parsimonious_placements\tneighborhood_size\n";
    }
    std::ofstream locfile;
    if (flocs != "") {
        locfile.open(flocs);
        locfile << "placement\tsample\n";
    }

    //mapper code wants a pointer.
    MAT::Tree* T = &Tobj;

    std::vector<std::string> samples;
    //read in the samples files and get the nodes corresponding to each sample.
    if (sample_file != "") {
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
            //fdfs.emplace_back(T->get_node(sample));
            samples.emplace_back(sample);
        }
        fprintf(stderr, "Processing %ld samples\n", samples.size());
    } else {
        fprintf(stderr, "WARNING: No sample file indicated; calculating for full tree\n");
        samples = T->get_leaves_ids();
    }
    //this loop is not a parallel for because its going to contain a parallel for.
    //this specific function would probably be better optimized if the outer was parallelized and the inner was not
    //but having the inner loop parallelized lets me use parallelization when calculating EPPs for very few or single samples in the future

    for (size_t s=0; s<samples.size(); s++) {
        //get the node object.
        auto node = T->get_node(samples[s]);
        size_t num_best;
        size_t neighborhood_size;
        auto best_placements = findEPPs(&Tobj, node, &num_best, &neighborhood_size);
        if (fepps != "") {
            eppfile << node->identifier << "\t" << num_best << "\t" << neighborhood_size << "\n";
        }
        if (flocs != "") {
            locfile << node->identifier << "\t" << node->identifier << "\n";
            for (auto pn: best_placements) {
                locfile << pn->identifier << "\t" << node->identifier << "\n";
            }
        }
    }
    if (fepps != "") {
        eppfile.close();
    }
    if (flocs != "") {
        locfile.close();
    }

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

//a variation on sample selection specific to uncertainty metrics.
std::vector<std::string> get_samples_epps (MAT::Tree* T, size_t max_epps, std::vector<std::string> to_check) {
    //calculate uncertainty for all samples in the tree
    //and return the set of samples which have EPPs less than max_epps
    //default filter value is 1, which 85% of samples have
    std::vector<std::string> good_samples;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        //check every sample if the ones to check is unset, else only calculate for the input sample set to_check
        if (to_check.size() == 0 || std::find(to_check.begin(), to_check.end(), n->identifier) != to_check.end()) {
            size_t nb;
            size_t ns;
            auto placements = findEPPs(T, n, &nb, &ns);
            if (nb <= max_epps) {
                good_samples.push_back(n->identifier);
            }
        }
    }
    return good_samples;
}

po::variables_map parse_uncertainty_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description ann_desc("uncertainty options");
    ann_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("samples,s", po::value<std::string>()->default_value(""),
     "Path to a simple text file of sample names to calculate uncertainty for.")
    ("find-epps,e", po::value<std::string>()->default_value(""),
     "Name for an Auspice-compatible tsv file output of equally parsimonious placements and neighborhood sizes for input samples.")
    ("record-placements,o", po::value<std::string>()->default_value(""),
     "Name for an Auspice-compatible two-column tsv which records potential parents for each sample in the query set.")
    ("dropout-mutations,d", po::value<std::string>()->default_value(""),
     "Name a file to calculate and save mutations which may be associated with primer dropout [EXPERIMENTAL].")
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

double fisher_test(unsigned a, unsigned b, unsigned c, unsigned d) {
    //based on https://github.com/usuyama/fisher_exact_test/blob/master/fisher.hpp
    unsigned N = a + b + c + d;
    unsigned r = a + c;
    unsigned n = c + d;
    unsigned max_for_k;
    if (r < n) {
        max_for_k = r;
    } else {
        max_for_k = n;
    }
    unsigned min_for_k;
    if (0 >= int(r+n-N)) {
        min_for_k = unsigned(0);
    } else {
        min_for_k = unsigned(int(r+n-N));
    }
    boost::math::hypergeometric_distribution<> hgd(r, n, N);
    double cutoff = boost::math::pdf(hgd, c);
    double tmp_p = 0.0;
    for(int k = min_for_k; k < static_cast<int>(max_for_k + 1); k++) {
        double p = boost::math::pdf(hgd, k);
        if(p <= cutoff) tmp_p += p;
    }
    return tmp_p;
}

std::map<std::string,size_t> get_mutation_count(MAT::Tree* T, MAT::Node* A = NULL, bool by_location = false) {
    std::map<std::string,size_t> mcm;
    for (auto n: T->depth_first_expansion(A)) {
        for (auto m: n->mutations) {
            std::string id;
            if (by_location) {
                id = std::to_string(m.position);
            } else {
                id = m.get_string();
            }
            if (mcm.find(id) == mcm.end()) {
                mcm[id] = 1;
            } else {
                mcm[id]++;
            }
        }
    }
    return mcm;
}

void check_for_droppers(MAT::Tree* T, std::string outf) {
    timer.Start();
    std::ofstream outfile (outf);
    outfile << "mutation\tbranch\tpvalue\tcorrected_pvalue\toccurrences_in\toccurrences_out\tsplit_size\tlocation_pvalue\tlocation_corrected_pvalue\n";
    std::map<std::string,double> pvals;
    std::map<std::string,std::string> nodetrack;
    std::map<std::string,size_t> ocintrack;
    std::map<std::string,size_t> splitstrack;
    std::map<std::string,double> lpvals;
    //first, get the overall mutation map for the global tree.
    auto gmap = get_mutation_count(T, NULL, false);
    auto locmap = get_mutation_count(T, NULL, true);
    size_t global_parsimony_score = 0;
    for (auto kv: gmap) {
        global_parsimony_score += kv.second;
    }
    size_t tests_performed = 0;
    size_t loc_tests_performed = 0;
    for (auto n: T->depth_first_expansion()) {
        //timer.Start();
        //for each split point, get the subtree and the mutation count of that subtree.
        auto lmap = get_mutation_count(T, n, false);
        size_t local_parsimony_score = 0;
        for (auto kv: lmap) {
            local_parsimony_score += kv.second;
        }
        if (local_parsimony_score < 50) {
            //totally arbitrary cutoff. No clue how important this will be.
            continue;
        }
        auto mloc = get_mutation_count(T, n, true);
        //do a fisher's exact test on the counts of this mutation vs the total parsimony score of the outside vs the inside
        //for each mutation in lmap.
        for (auto kv: lmap) {
            if (kv.second < 10) {
                //more arbitrary cutoffs.
                continue;
            }
            auto pv = fisher_test(kv.second, local_parsimony_score, gmap[kv.first]-kv.second, global_parsimony_score-local_parsimony_score);
            tests_performed++;
            if (pv < 0.05) {
                //perform a secondary location-based test.
                //this will have less total tests performed as it is conditioned on the previous test.
                std::string locstr = kv.first.substr(1, kv.first.size()-2);
                auto lpv = fisher_test(mloc[locstr], local_parsimony_score, locmap[locstr]-mloc[locstr], global_parsimony_score-local_parsimony_score);
                //std::cerr << "DEBUG: loc fisher " << lpv << " performed with: " << mloc[locstr] << "," << local_parsimony_score << "," << locmap[locstr]-mloc[locstr] << "," << global_parsimony_score-local_parsimony_score << "\n";
                loc_tests_performed++;
                if (pvals.find(kv.first) == pvals.end()) {
                    pvals[kv.first] = pv;
                    lpvals[kv.first] = lpv;
                    nodetrack[kv.first] = n->identifier;
                    ocintrack[kv.first] = kv.second;
                    splitstrack[kv.first] = local_parsimony_score;
                } else if (pv < pvals[kv.first]) {
                    pvals[kv.first] = pv;
                    lpvals[kv.first] = lpv;
                    nodetrack[kv.first] = n->identifier;
                    ocintrack[kv.first] = kv.second;
                    splitstrack[kv.first] = local_parsimony_score;
                }
            }
        }
    }
    for (auto kv: pvals) {
        outfile << kv.first << "\t" << nodetrack[kv.first] << "\t" << pvals[kv.first] << "\t" << (pvals[kv.first] * tests_performed) << "\t" << ocintrack[kv.first] << "\t" << gmap[kv.first] - ocintrack[kv.first] << "\t" << splitstrack[kv.first] << "\t" << lpvals[kv.first] << "\t" << (lpvals[kv.first] * loc_tests_performed) << "\n";
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

void uncertainty_main(po::parsed_options parsed) {
    po::variables_map vm = parse_uncertainty_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string sample_file = vm["samples"].as<std::string>();
    std::string fepps = vm["find-epps"].as<std::string>();
    std::string flocs = vm["record-placements"].as<std::string>();
    std::string dropmuts = vm["dropout-mutations"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T in this scope is the actual object and not a pointer
    if (T.condensed_nodes.size() > 0) {
        T.uncondense_leaves();
    }
    if (dropmuts != "") {
        fprintf(stderr, "Identifying primer-dropout associated mutations.\n");
        check_for_droppers(&T, dropmuts);
    }
    if (sample_file != "") {
        fprintf(stderr, "Calculating placement uncertainty\n");
        findEPPs_wrapper(T, sample_file, fepps, flocs);
    }
}
