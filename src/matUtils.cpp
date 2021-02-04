#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
//#include "usher_graph.hpp"
#include "usher_mapper.cpp"//it shouldn't need this, but if I just include the hpp it insists that the mapper2_body function isn't defined. 

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

po::variables_map check_options(int argc, char** argv) {
    // Check command line options and return variable map.
    po::options_description desc{"Options"};
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file to mask [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->required(),
         "Output masked mutation-annotated tree file [REQUIRED]")
        ("restricted-samples,s", po::value<std::string>()->default_value("none"), //this will buf if they name their restricted sample file "none" with no extension for some insane reason
         "Sample names to restrict. Use to perform masking") //this is now optional, as the Utils may be doing things other than masking. Should still perform the same given the same commands as previous.
        ("find-epps,e", po::bool_switch(),
        "Use to calculate and store the number of equally parsimonious placements for all nodes")
        ("help,h", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(p)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

/*
As a general principle, I intend to write this utility such that each function is modular and self-contained and any or all of them can be called based on command line usage.
Each relevant function will both take and return a MAT object
The main() function will only contain MAT and option read in, a series of if() then function calls, and saving the tree at the end
*/

MAT::Tree restrictSamples (std::string samples_filename, MAT::Tree T) {
    // Load restricted sampl0e names from the input file and add it to the set
    std::ifstream infile(samples_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the restricted samples file: %s!\n", samples_filename.c_str());
        exit(1);
    }    
    std::unordered_set<std::string> restricted_samples;
    std::string sample;
    while (std::getline(infile, sample)) {
        if (T.get_node(sample) == NULL) {
            fprintf(stderr, "ERROR: Sample %s missing in input MAT!\n", sample.c_str());
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

MAT::Tree findEPPs (MAT::Tree Tobj) {
    TIMEIT()
    //all comments on function JDM
    //first, need to iterate through all leaf nodes, including condensed nodes and single leaves
    //internal nodes have metadata objects but those are just gonna be default values 0 for epps for now which should indicate high confidence anyways
    //the simplest way to do this is to iterate through all nodes and check whether each one is a leaf
    //I'm having segfault problems, which I'm going to guess have to do with how I'm editing things while iterating over them?
    //maybe try making a copy of the tree and popping from the copy, then updating the node in the original without ever creating a node?

    
    //MAT::Tree* T = &Tobj; //T in this function is a pointer because of how the mapper is written.
    //fprintf(stderr, "Number of Leaves on Tree Pre-everything %d \n", Tobj.get_num_leaves());
    //fprintf(stderr, "Assigned tree pointer\n");
    MAT::Tree TCopy = MAT::get_tree_copy(Tobj);
    MAT::Tree* T = &TCopy; //wants the pointer for mapping.
    auto fdfs = Tobj.depth_first_expansion(); //the full tree expanded for outer loop iteration.
    for (size_t s=0; s<fdfs.size(); s++){ //this loop is not a parallel for because its going to contain a parallel for
        //get the node object.
        auto node = fdfs[s];
        if (node->is_leaf()) {
            //retrieve the full set of mutations associated with this Node object from root to it
            //to do this, get the full set of ancestral nodes and their mutations
            //code copied from the usher mapper.
            std::vector<int> anc_positions; //tracking positions is required to account for backmutation/overwriting along the path
            std::vector<MAT::Mutation> ancestral_mutations;
            //first load in the current mutations
            for (auto m: node->mutations){
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
            //fprintf(stderr, "Loaded mutations for leaf into array\n");

            //JDM-create a copy of the tree for pruning and mapping.

            //fprintf(stderr, "Number of Leaves on Clone Tree Premapping %d \n", T->get_num_leaves());

            //now that we have the mutation set, pop the current node from the Tree
            //save its identifier, parent, and branch length for replacement
            //fprintf(stderr, "Preremoval Parent ID %s \n", node->parent->identifier.c_str());
            //const std::string &nparid = node->parent->identifier;
            MAT::Node* nparn = node->parent;
            const std::string &nid = node->identifier; //needed to relocate the original node
            float blen = node->branch_length;

            T->remove_node(node->identifier, true); //this should modify in-place. pop it from the copy
            //fprintf(stderr, "Postremoval Parent ID %s \n", nparid.c_str());

            //the ancestral_mutations vector, plus the mutations assigned to this specific node, constitute the "missing_sample" equivalents for calling the mapper
            //PUT THE MAPPING HERE?? THERES NO CLEAR FUNCTION TO DO THIS :( WHY YATISH, WHY WRITE IT SO HORRIBLY
            //just super confusing repetitive ass code from the main Usher cpp, ugh
            //I think I'll have to loop over every node in the pruned tree and check placement???
            //fprintf(stderr, "Node removed from the tree\n");
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
                    }); 
            //BACK TO MY CODE (JDM)
            //put the node back.
            //fprintf(stderr, "Node %s remapped\n", nid.c_str()); //seems fine...
            //fprintf(stderr, "Node branch length: %f \n", blen);
            //fprintf(stderr, "Postmapping Parent ID %s \n", nparid.c_str());
            //fprintf(stderr, "Number of Leaves on Clone Tree Postmapping %d \n", T->get_num_leaves());
            MAT::Node* repnode = NULL;
            //fprintf(stderr, "Created blank node object\n");
            repnode = T->create_node(nid, nparn, blen); //replace the node back onto the tree copy
            //can't add error messages to the create node function because it's called a bajillion times in the mapper and it floods me out
            //delete the tree copy.
            //fprintf(stderr, "Attempting to delete tree copy\n");
            //delete T;
            //fprintf(stderr, "Tree copy deleted\n");
            //fprintf(stderr, "Node recreated\n");
            //give the original tree node, which hasn't moved, the metadata.
            auto cnode = Tobj.get_node(nid);
            cnode->epps = num_best;
            //fprintf(stderr, "Node metadata updated\n");

        }
    }
    return Tobj; //return the actual object.
}

int main(int argc, char** argv) {

    // Command line options
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string samples_filename = vm["restricted-samples"].as<std::string>();
    bool fepps = vm["find-epps"].as<bool>();
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    // If a restricted samples file was provided, perform masking procedure
    if (samples_filename != "none") {
        T = restrictSamples(samples_filename, T);
    }
    // If the argument to calculate equally parsimonious placements was used, perform this operation
    if (fepps) {
        fprintf(stderr, "Attempting to calculate EPPs\n");
        T = findEPPs(T);
    }

    // Store final MAT to output file 
    MAT::save_mutation_annotated_tree(T, output_mat_filename);

    return 0;
}

