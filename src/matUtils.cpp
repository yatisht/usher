#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
#include "usher_graph.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

Timer timer; 

void assignLineages (MAT::Tree& T, const std::string& lineage_filename, float min_freq) {
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
            if (static_cast<float>(mc.second)/it.second.size() > min_freq) {
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
            size_t freq;
            Node_freq (size_t a, size_t b) 
            {
                best_j = a;
                freq = b;
            }
            // To sort with highest frequency first
            inline bool operator< (const Node_freq& n) const {
                return ((*this).freq > n.freq);
            }
        };

        std::vector<Node_freq> best_node_frequencies;
        if (num_best > 1) {
            fprintf(stderr, "WARNING: found %zu possible assignments\n", num_best);
            for (auto j: best_j_vec) {
                size_t freq = 0;

                tbb::parallel_for(tbb::blocked_range<size_t>(0, it.second.size()),
                        [&](const tbb::blocked_range<size_t> r) {
                        for (size_t i=r.begin(); i<r.end(); ++i){
                            if (T.is_ancestor(dfs[j]->identifier, it.second[i]->identifier)) {
                                __sync_fetch_and_add(&freq, 1);
                            }
                        }
                    }, ap);

                Node_freq n(j, freq);
                best_node_frequencies.push_back(n);

            }
        }
        else {
            Node_freq n(best_j, 0);
            best_node_frequencies.push_back(n);
        }
        
        std::sort(best_node_frequencies.begin(), best_node_frequencies.end());

        bool assigned = false;
        for (auto n: best_node_frequencies) {
            auto j = n.best_j;
            if (dfs[j]->clade == "") {
                fprintf(stderr, "Assigning %s to node %s\n", it.first.c_str(), dfs[j]->identifier.c_str());
                dfs[j]->clade = it.first;
                assigned = true;
                break;
            }
        }

        if (!assigned) {
            fprintf(stderr, "WARNING: Could not assign a node to clade %s since all possible nodes were already assigned some other clade!\n", it.first.c_str());
        }

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}

//Filter subcommands follow
MAT::Tree restrictSamples (std::string samples_filename, MAT::Tree T) {
    // Load restricted sample names from the input file and add it to the set
    //BUG WARNING
    //This does not actually work right now- all samples passed in are caught as missing. 
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
            }
            else {
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
    return T;
}

//Annotate subcommands follow

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

//Convert subcommands follow

void write_vcf_header(FILE *vcf_file, std::vector<Mutation_Annotated_Tree::Node*> &dfs,
                      bool print_genotypes) {
    // Write minimal VCF header with sample names in same order that genotypes
    // will be printed out (DFS).
    fprintf(vcf_file, "##fileformat=VCFv4.2\n");
    fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if (print_genotypes) {
        fprintf(vcf_file, "\tFORMAT");
        for (auto node: dfs) {
            if (node->is_leaf()) {
                fprintf(vcf_file, "\t%s", node->identifier.c_str());
            }
        }
    }
    fputc('\n', vcf_file);
}

uint count_leaves(std::vector<Mutation_Annotated_Tree::Node*> &dfs) {
    // Return the number of leaf nodes in dfs
    uint count = 0;
    for (auto node: dfs) {
        if (node->is_leaf()) {
            count++;
        }
    }
    return count;
}

int8_t *new_gt_array(int size, int8_t ref) {
    // Allocate and return an array of int8_t (encoded nucleotide values) initialized to ref.
    int8_t *gt_array = new int8_t[size];
    for (int i = 0;  i < size;  i++) {
        gt_array[i] = ref;
    }
    return gt_array;
}

uint r_add_genotypes(MAT::Node *node,
                     std::unordered_map<std::string, std::vector<int8_t *>> &chrom_pos_genotypes, 
                     std::unordered_map<std::string, std::vector<int8_t>> &chrom_pos_ref,
                     uint leaf_count, uint leaf_ix, std::vector<struct MAT::Mutation *> &mut_stack) {
    // Traverse tree, adding leaf/sample genotypes for mutations annotated on path from root to node
    // to chrom_pos_genotypes (and reference allele to chrom_pos_ref).
    for (auto &mut: node->mutations) {
      if (mut.is_masked()) {
          continue;
      }
      mut_stack.push_back(&mut);
    }
    if (node->is_leaf()) {
        // Store genotypes in this leaf's column for all mutations on the path from root to leaf
        for (auto mut: mut_stack) {
            std::string chrom = mut->chrom;
            uint pos = (uint)mut->position;
            if (chrom.empty()) {
              fprintf(stderr, "mut->chrom is empty std::string at node '%s', position %u\n",
                      node->identifier.c_str(), pos);
            }
            if (chrom_pos_genotypes.find(chrom) == chrom_pos_genotypes.end()) {
                // First variant on chrom: initialize a vector mapping position to genotype.
                // Assume a genome size similar to SARS-CoV-2, resize if necessary.
                uint initSize = 30000;
                chrom_pos_genotypes[chrom] = std::vector<int8_t *>(initSize);
                chrom_pos_ref[chrom] = std::vector<int8_t>(initSize);
            }
            if (pos >= chrom_pos_genotypes[chrom].size()) {
                // chrom has larger positions than we assumed; allocate a larger vector.
                uint newSize = chrom_pos_genotypes[chrom].size() * 2;
                chrom_pos_genotypes[chrom].resize(newSize);
                chrom_pos_ref[chrom].resize(newSize);
            }
            if (! chrom_pos_genotypes[chrom][pos]) {
                // First variant reported at this position; allocate genotype array and
                // store reference allele (which is stored in par_nuc not ref_nuc).
                chrom_pos_genotypes[chrom][pos] = new_gt_array(leaf_count, mut->par_nuc);
                chrom_pos_ref[chrom][pos] = mut->par_nuc;
            }
            // Store the allele/genotype for this chrom / pos / sample.
            chrom_pos_genotypes[chrom][pos][leaf_ix] = mut->mut_nuc;
        }
        leaf_ix++;
    }
    for (auto child: node->children) {
        leaf_ix = r_add_genotypes(child, chrom_pos_genotypes, chrom_pos_ref, leaf_count, leaf_ix,
                                  mut_stack);
    }
    for (auto mut: node->mutations) {
        mut_stack.pop_back();
    }
    return leaf_ix;
}

std::unordered_map<int8_t, uint>count_alleles(int8_t *gt_array, uint gtCount)  {
    // Tally up the count of each allele (both ref and alts) from sample genotypes.
    std::unordered_map<int8_t, uint> allele_counts;
    for (uint i = 0;  i < gtCount;  i++) {
        int8_t allele = gt_array[i];
        if (allele_counts.find(allele) == allele_counts.end()) {
            allele_counts.insert({allele, 1});
        } else {
            allele_counts[allele]++;
        }
    }
    return allele_counts;
}

bool cmp_allele_count_desc(const std::pair<int8_t, uint>& a, const std::pair<int8_t, uint>& b) {
    // Compare counts of two alleles, for sorting in descending order.
    return a.second > b.second;
}

std::map<int8_t, uint>make_alts(std::unordered_map<int8_t, uint> &allele_counts, int8_t ref) {
    // Map alternate alleles, ordered by count (highest first), to counts.
    std::vector<std::pair<int8_t, uint>> pairs;
    for (auto &itr : allele_counts) {
        if (itr.first != ref) {
            pairs.push_back(itr);
        }
    }
    std::sort(pairs.begin(), pairs.end(), cmp_allele_count_desc);
    std::map<int8_t, uint> alts;
    for (auto &itr : pairs) {
      alts.insert(itr);
    }
    return alts;
}

std::string make_id(int8_t ref, uint pos, std::map<int8_t, uint> &alts) {
    // Return a C std::string comma-sep list of the form <ref><pos><alt1>[,<ref><pos><alt2>[,...]].
    std::string id;
    for (auto &itr : alts) {
        if (! id.empty()) {
            id += ",";
        }
        id += MAT::get_nuc(ref) + std::to_string(pos) + MAT::get_nuc(itr.first);
    }
    return id;
}

std::string make_alt_str(std::map<int8_t, uint> &alts) {
    // Return a C std::string comma-sep list of alternate alleles.
    std::string alt_str;
    for (auto &itr : alts) {
        if (! alt_str.empty()) {
          alt_str += ",";
        }
        alt_str += MAT::get_nuc(itr.first);
    }
    return alt_str;
}

std::string make_info(std::map<int8_t, uint> &alts, uint leaf_count) {
    // Return a C std::string VCF INFO value with AC (comma-sep list of alternate allele counts)
    // and AN (total genotype count).
    std::string alt_count_str;
    for (auto &itr : alts) {
        if (! alt_count_str.empty()) {
            alt_count_str += ",";
        }
        alt_count_str += std::to_string(itr.second);
    }
    std::string info = "AC=" + alt_count_str + ";AN=" + std::to_string(leaf_count);
    return info;
}

int *make_allele_codes(int8_t ref, std::map<int8_t, uint> &alts) {
    // Return an array that maps binary-encoded nucleotide to VCF genotype encoding:
    // 0 for reference allele, 1 for first alternate allele, and so on.
    int *al_codes = new int[256];
    for (int i = 0;  i < 256;  i++) {
        al_codes[i] = 0;
    }
    al_codes[(uint8_t)ref] = 0;
    int altIx = 1;
    for (auto &itr : alts) {
        al_codes[itr.first] = altIx++;
    }
    return al_codes;
}

void write_vcf_rows(FILE *vcf_file, MAT::Tree T, std::vector<MAT::Node*> &dfs, bool print_genotypes) {
    // Fill in a matrix of genomic positions and sample genotypes in the same order as the
    // sample names in the header, compute allele counts, and output VCF rows.
    uint leaf_count = count_leaves(dfs);
    // The int8_t here is mutation_annotated_tree.hpp's binary encoding of IUPAC nucleotide bases.
    std::unordered_map<std::string, std::vector<int8_t *>> chrom_pos_genotypes;
    std::unordered_map<std::string, std::vector<int8_t>> chrom_pos_ref;
    std::vector<struct MAT::Mutation *> mut_stack;
    r_add_genotypes(T.root, chrom_pos_genotypes, chrom_pos_ref, leaf_count, 0, mut_stack);
    // Write row of VCF for each variant in chrom_pos_genotypes[chrom]
    for (auto itr = chrom_pos_genotypes.begin();  itr != chrom_pos_genotypes.end();  ++itr) {
        std::string chrom = itr->first;
        std::vector<int8_t *> pos_genotypes = itr->second;
        for (uint pos = 0;  pos < pos_genotypes.size();  pos++) {
            int8_t *gt_array = pos_genotypes[pos];
            if (gt_array) {
              int8_t ref = chrom_pos_ref[chrom][pos];
              std::unordered_map<int8_t, uint>allele_counts = count_alleles(gt_array, leaf_count);
              std::map<int8_t, uint>alts = make_alts(allele_counts, ref);
              std::string id = make_id(ref, pos, alts);
              std::string alt_str = make_alt_str(alts);
              std::string info = make_info(alts, leaf_count);
              fprintf(vcf_file, "%s\t%d\t%s\t%c\t%s\t.\t.\t%s",
                      chrom.c_str(), pos, id .c_str(), MAT::get_nuc(ref), alt_str.c_str(),
                      info.c_str());
              if (print_genotypes) {
                  int *allele_codes = make_allele_codes(ref, alts);
                  fprintf(vcf_file, "\tGT");
                  for (uint i = 0;  i < leaf_count;  i++) {
                      int8_t allele = gt_array[i];
                      fprintf(vcf_file, "\t%d", allele_codes[allele]);
                  }
              }
              fputc('\n', vcf_file);
            }
        }
    }
}

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes) {
    auto vcf_filepath = "./" + vcf_filename; //taking away a touch of functionality for simplicity. just write it to this directory.
    FILE *vcf_file = fopen(vcf_filepath.c_str(), "w");
    std::vector<Mutation_Annotated_Tree::Node*> dfs = T.depth_first_expansion();
    write_vcf_header(vcf_file, dfs, !no_genotypes);
    write_vcf_rows(vcf_file, T, dfs, !no_genotypes);
    fclose(vcf_file);
}

//main and parsing functions follow

po::variables_map parse_annotate_command(po::parsed_options parsed) {    

    po::variables_map vm;
    po::options_description ann_desc("annotate options");
    ann_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->required(),
         "Path to output processed mutation-annotated tree file [REQUIRED]")
        ("lineage-names,l", po::value<std::string>()->required(),
         "Path to a file containing lineage asssignments of samples. Use to locate and annotate clade root nodes")
        ("allele-frequency,f", po::value<float>()->default_value(0.9),
         "Minimum allele frequency in input samples for finding the best clade root. Used only with -l")
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
    //    bool get_parsimony = vm["get-parsimony"].as<bool>();
    //    bool fepps = vm["find-epps"].as<bool>();

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
        assignLineages(T, lineage_filename, allele_frequency);
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

po::variables_map parse_filter_command(po::parsed_options parsed) {

    po::variables_map vm;
    po::options_description filt_desc("filter options");
    filt_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->required(),
         "Path to output filtered mutation-annotated tree file [REQUIRED]")
        ("restricted-samples,s", po::value<std::string>()->default_value(""), 
         "Sample names to restrict. Use to perform masking") 
//        ("placement-confidence,c", po::value<int>()->default_value(0),
//        "Maximum number of equally parsimonious placements among nodes included in the tree (lower values is better, with 1 as highest confidence). Set to 0 to skip filtering. Default 0")
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << filt_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void filter_main(po::parsed_options parsed) {
    //the filter subcommand prunes data from the protobuf based on some threshold, returning a protobuf file that is smaller than the input
    po::variables_map vm = parse_filter_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string samples_filename = vm["restricted-samples"].as<std::string>();
    //int maxcon = vm["placement-confidence"].as<int>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    // If a restricted samples file was provided, perform masking procedure
    if (samples_filename != "") {
        fprintf(stderr, "Performing Masking\n");
        T = restrictSamples(samples_filename, T);
    }
    //there's a very simple filtering procedure which doesn't require a dedicated function
    // if (maxcon > 0) { //value of 0 means skip this procedure (default behavior)
        //fprintf(stderr, "Removing nodes with more than %d equally parsimonious placements", maxcon);
        //just get all nodes, and for each one with an EPPs greater than maxcon, remove it.
        //for the smallest possible maximum, 1, this should retain about 84% of samples. More for any other value
        // auto dfs = T.depth_first_expansion(); //technically for now I would want to get just all leaves, but I think the epps concept could be extended to internal nodes that aren't true samples with some thought in the future
        //for (auto it: dfs) {
            //if (it->epps > maxcon) {
            //    T.remove_node(it->identifier, false); //fairly sure I want this to be false for leaves
            //}
        //}
    //}

    // Store final MAT to output file
    if (output_mat_filename != "") {
        fprintf(stderr, "Saving Final Tree\n");
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
    }    
}

po::variables_map parse_convert_command(po::parsed_options parsed) {

    po::variables_map vm;
    po::options_description conv_desc("convert options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("write-vcf,v", po::value<std::string>()->default_value(""),
         "Output VCF file ")
        ("no-genotypes,n", po::bool_switch(),
        "Do not include sample genotype columns in VCF output. Used only with the vcf option")
        ("write-tree,t", po::value<std::string>()->default_value(""),
         "Use to write a newick tree to the indicated file.")
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void convert_main(po::parsed_options parsed) {
    //the convert subcommand converts the protobuf into another file format, returning some other file type that represents an equivalent structure
    po::variables_map vm = parse_convert_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string tree_filename = vm["write-tree"].as<std::string>();
    std::string vcf_filename = vm["write-vcf"].as<std::string>();
    bool no_genotypes = vm["no-genotypes"].as<bool>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    //if a vcf filename was given, write a vcf to it
    if (vcf_filename != "") {
        fprintf(stderr, "Generating VCF\n");
        make_vcf(T, vcf_filename, no_genotypes);
    }
    //if a newick tree filename was given, write a tree to it
    if (tree_filename != "") {
        fprintf(stderr, "Generating Newick file\n");
        auto tree_filepath = "./" + tree_filename; //for simplicity, write it to this directory
        FILE *tree_file = fopen(tree_filepath.c_str(), "w");
        fprintf(tree_file, "%s\n",
            MAT::get_newick_string(T, true, true, true).c_str());
        fclose(tree_file);        
    }
}

int main (int argc, char** argv) {
    /*
    The new design principle for organizing matUtils is to divide the overall structure into three options. All three take protobuf files as input.
    First there is annotate, which calculates and stores uncertainty metrics and other tree-based metadata into the protobuf file. This returns a protobuf that is larger than the input (in bytes).
    Second there is filter, which prunes the tree based on threshold arguments and the uncertainty metrics from annotate. This returns a protobuf that is smaller than the input (in bytes).
    Third there is convert, which produces a different file format than protobuf. This returns non-protobuf files.

    Generally a workflow will call annotate first, then filter, then convert. Or, if metadata is precalculated and saved on our updated global tree, they should be able to call filter then convert without an extended annotate command.
    Ideally these would be chainable on the command line, though I'm not sure if we can write a protobuf file to stdout/read from stdin. This is something we may want to investigate as a future option.
    For now, they can use && and intermediate file names to write a simple matUtils pipeline.
    */
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, filter, convert.")
        ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    
    try {
        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();

        po::store(parsed, vm);
        std::string cmd = vm["command"].as<std::string>();
        if (cmd == "annotate"){
            annotate_main(parsed);
        } else if (cmd == "convert"){
            convert_main(parsed);
        } else if (cmd == "filter"){
            filter_main(parsed); 
        } else if (cmd == "help" || cmd == "--help" || cmd == "-h") { //trying to catch some of the intuitive things people will try.
            fprintf(stderr, "matUtils has three major subcommands: annotate, filter, and convert. All three take a MAT .pb file as input.\nAnnotate adds information to the MAT. Use annotate when you have information you want to calculate or incorporate into the MAT .pb. The command 'matUtils annotate --help' will describe related options.\nFilter removes nodes or samples from the MAT. Use filter when you want to strip out low-quality samples or mask samples you want to avoid. The command 'matUtils filter --help' will describe related options.\nConvert produces files which are not MAT .pb format, such as .vcf or newick text files. Use convert when you want other file types. The command 'matUtils convert --help' will describe related options.\n"); //very open to alternative wording/formatting to make this nicer/clearer.
            exit(0);
        } else {
            fprintf(stderr, "Invalid command. Please choose from annotate, filter, convert, or help and try again.\n");
            exit(1);
        }
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Please choose from annotate, filter, convert, or help and try again.\n");
        exit(0);
    }

    return 0;
}
