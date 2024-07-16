#include "mask.hpp"
#include <random>
#include <algorithm>
#include "maskselect.hpp"
#include <omp.h>

using namespace std;

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
    //("snp-distance,d", po::value<uint32_t>()->default_value(0),
    // "SNP distance between a sample and the internal node which will have all descendents masked for missing data.")
    ("max-snp-distance,D", po::value<uint32_t>()->default_value(0),
     "Maximum distance past snp-distance that can be used to find a local ancestor, if local ancestor exists past the max-SNP-distance, program will look for local ancestor within bounds.")
    ("diff-file,f", po::value<std::string>()->default_value(""),
    "Diff files for samples contained in the tree. Samples not included will not be considered in masking.")
    //("ignore-positions-file,p", po::value<std::string>()->default_value(""),
    //"Diff files for samples contained in the tree. Samples not included will not be considered in masking.")
    ("move-nodes,M", po::value<std::string>()->default_value(""),
     "Name of the TSV file containing names of the nodes to be moved and their new parents. Use to move nodes around the tree between paths containing identical sets of mutations.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    //std::cout << "opts: " << opts << std::endl;
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        fprintf(stderr, "stuck here\n");
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
    //uint32_t snp_distance = vm["snp-distance"].as<uint32_t>();
    uint32_t max_snp_distance = vm["max-snp-distance"].as<uint32_t>();
    std::string diff_file = vm["diff-file"].as<std::string>();
    //std::string pos_file = vm["ignore-positions-file"].as<std::string>();
    tbb::task_scheduler_init init(num_threads);
    fprintf(stderr, "made it to main function");

    //check for mutually exclusive arguments
    //LILY: make sure you check for need for exclusivity of your function 
    if ((simplify) & (rename_filename != "")) {
        //doesn't make any sense to rename nodes after you just scrambled their names. Or to rename them, then scramble them.
        fprintf(stderr, "ERROR: Sample renaming and simplification are mutually exclusive operations. Review argument choices\n");
        exit(1);
    }
    if ((max_snp_distance > 0) & (diff_file == "")) {
        //doesn't make any sense to rename nodes after you just scrambled their names. Or to rename them, then scramble them.
        fprintf(stderr, "ERROR: Must provide diff file of samples for local masking. Review argument choices\n");
        exit(1);
    }
     
    // Load input MAT and uncondense tree
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    timer.Start();
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
        fprintf(stderr, "Uncondensing condensed nodes.\n");
        timer.Start();
        T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
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

    if (max_snp_distance > 0) {
        bool load_all = true;
        fprintf(stderr, "made it to here");
        localMask(max_snp_distance, T, diff_file, output_mat_filename, num_threads);
    }

    // Store final MAT to output file
    if (output_mat_filename != "") {
        if (recondense) {
            timer.Start();
            fprintf(stderr, "Collapsing tree...\n");
            T.collapse_tree();
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
            timer.Start();
            fprintf(stderr, "Condensing leaves...\n");
            T.condense_leaves();
            fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        }
        fprintf(stderr, "Saving Final Tree to %s\n", output_mat_filename.c_str());
        timer.Start();
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}

std::map<std::string, std::map<int, int>> readDiff (const std::string& diff_file) {
    fprintf(stderr, "Reading Variation Information\n");
    //only storing missing data, stored as position and length of missing
    std::map<std::string, std::map<int, int>> data;

    try {
        //open the file
        std::ifstream file(diff_file);
        //determine if there is a file error
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + diff_file);
        }

        std::string line; // initialize line variable to store each line from file
        std::string current_sample; // initialize sample name variable to track samples in file
        
        //iterate through all lines
        while (std::getline(file, line)) {
            
            std::vector<std::string> substrings; //store parsed line in a vector
            size_t startPos = 0;
            size_t endPos;

            //new sample, set up map
            if (line[0] == '>') {
                current_sample = line.erase(0,1); // format line to remove '>'
                auto it = data.find(current_sample);
                if (it == data.end()) {
                    // Inner map doesn't exist for the current sample
                    // Create a new inner map for the current sample
                    data[current_sample] = std::map<int, int>();
                } else {
                    // Inner map already exists for the current sample
                    // throw an error 
                    throw std::runtime_error("Duplicate samples detected, inspect diff file for sample: " + current_sample);
                }                
            }
            //Find lines with missing data
            else if (line[0] == '-'){
                //find tab separators to create substrings
                while ((endPos = line.find('\t', startPos)) != std::string::npos) {
                    
                    //Extract the substring between startPos and endPos and add it to the substrings vector
                    substrings.push_back(line.substr(startPos, endPos - startPos));
                    
                    // Update startPos to the position after the separator
                    startPos = endPos + 1;                    
                }
            
            // Extract the substring after the last occurrence of the separator
            substrings.push_back(line.substr(startPos));

            //convert substrings to ints and add them to map for current sample
            int position = std::stoi(substrings[1]);
            int length = std::stoi(substrings[2]);
            data[current_sample][position] = length;    
            }
        }
        
        file.close();

    //more file error handling
    } catch (const std::exception& e) {
        std::cerr << "Exception caught while reading file: " << e.what() << std::endl;
        throw;  // Re-throw the exception to be handled where the function is called
    }
    
    fprintf(stderr, "All missing data retrieved\n");

    return data;
}

//new version 6/24
void nodeComp(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data) {

    //we are assessing the node mutations and the leaf missing 
    //this should be consistent thorugh whole code
    int node_len = node->mutations.size();
    //std::cout <<   "node len " << node_len << std::endl;
    //std::cout <<   "MADE IT TO node comp" << std::endl;

    //std::cout <<   "MISSING DATA" << missing_data.size() << std::endl;

    /*
    std::cout <<   "MUTATIONS " << node->identifier << std::endl;
    for (const auto& mut: node->mutations) {
            std::cout << "pos" << mut.get_string() << endl;
            std::cout << "pos" << mut.mut_nuc << endl;
            std::cout << "pos" << mut.ref_nuc << endl;
            std::cout << "pos" << mut.is_masked() << endl;
            std::cout << "pos" << mut.is_missing << endl;
    }
    */
    // the number of missing data points in a leaf, in a map
    //delete? 
    //std::cout <<   "leaf??" << leaf->identifier << std::endl;

    //auto leaf_len = diff_data[leaf->identifier].size();
    auto leaf_len = missing_data.size();
    
    // this is my iterator through the nodes
    //node counter goes up
    //int node_counter = 0;

    //node counter goes down 
    //delete?
    //int node_counter = node_len;

    // this is my leaf iterator
    //leaf iterator counts up 
    //auto leaf_counter = diff_data[leaf->identifier].begin();

    //leaf iterator counts down
    //auto leaf_counter = diff_data[leaf->identifier].end();
    auto leaf_counter = missing_data.end();
    // a separate iterator for the
    //prob dont need this?  
    auto node_it = node->mutations.end();

    //decrement to start in the right place
    node_it --;
    leaf_counter --;
    //leaf_counter1 --;
    //std::cout << "leaf_counter1" << leaf_counter1->first << std::endl;
    /*
    
    for (const auto& muttys: node->mutations) {
            std::cout << "muttys" << muttys.get_string() << endl;
    }
    */
    //i want to make this 
    //will delete this at some point 
    std::vector<MAT::Mutation> del_muts;
    bool last_iteration = false;

    //std::cout <<   "starting leaf counter " << leaf_counter->first  << std::endl;
    //std::cout <<   "starting node counter " << node_it  << std::endl;
    //if (node_len != 0 ){
    //        std::cout <<   "starting node counter " << node_it->get_string()  << std::endl;
    //    }
    //else {
    //    std::cout <<   "no mutations " << std::endl;
    //}

    //std::cout <<   "current node " << node->identifier << " " << node_len << std::endl;
    //std::cout <<   "leaf " << leaf->identifier << " " << leaf_len << std::endl;
    
    //counting down not up
    //one of these lists is still going
    //while (node_it != node->mutations.begin() || leaf_counter != diff_data[leaf->identifier].begin()){
    while (true) {
    //this one counts up 
    //while (node_counter < node_len || leaf_counter != diff_data[leaf->identifier].end()){
        
        //no mutations or missing info, end now dont waste time
        if (node_len == 0 || leaf_len == 0) {
            //std::cout << "loop ended " << std::endl;
            break;
        }
        //std::cout <<   "node counter " << node_it->get_string()  << std::endl;
        //both lists still going
        //cant use bool below bc i need to use beginning of list
        //if (node_it != node->mutations.begin() && leaf_counter != diff_data[leaf->identifier].begin()){
        
        //retrieve mutation
        std::string mutation = node_it->get_string();

        //std::cout << "mutation " << mutation << std::endl;
        //std::cout << "mutation " << node_it << std::endl;

        int mut_pos = stoi(mutation.substr(1, mutation.length() - 2));
        //std::cout << "mutation " << mutation << std::endl;
        //std::cout << "mutation position " << mut_pos << std::endl;

        //get position for missing block
        //int missing_start1 = leaf_counter1->first;
        int missing_start = leaf_counter->first;
        int missing_end = missing_start + leaf_counter->second;
        //std::cout << "missing " << missing_start << " " << missing_end << std::endl;    
        
        //check if mutation is inside missing
        if (mut_pos >= missing_start && mut_pos <= missing_end) {
            // may need to make this fancier later
            //std::cout << "mut is fully inside missing " << std::endl;
            //delete mutation here
            //probably store in a list too
            //std::cout << "DELETE: " << mutation << " from" << node->identifier << " leaf is " << leaf->identifier << std::endl;
            //std::cout << "node id before delete" << node_it->get_string() << std::endl; // Prints the type of x

            //del_muts.push_back(node->mutations[node_counter]);
            node->mutations.erase(node_it);
            //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
            //std::cout << "node id AFTERDELETE " << node_it->get_string() << std::endl;


            //add one to mutation index (there might be more mutations in missing region so dont index down missing)
            //node_counter -= 1;

            //do i need to decrement if i delete node_it? 
            //POSSIBLE ERROR HERE
            if (node_it != node->mutations.begin()) {
                node_it --;
            }
            else {
                //need to do something first?
                //std::cout << "end of mutations!!!! break " << std::endl;
                break;
            }
        }  
        //figure out if mutation is before or after missing
        //if after i need to decrement node_it to get it closer to missing, do nothing to missing
        else if (mut_pos > missing_end) {
            // might have to make this fancier later 
            //std::cout << "mut is fully after missing: node: " << node->identifier << std::endl;
            //std::cout << "mut is fully after missing: leaf: " << leaf->identifier << std::endl;
            //node_counter --;
            //make sure youre not decrementing into nothing
            if (node_it != node->mutations.begin()) {
                node_it --;
            }
            else {
                //need to do something first?
                //std::cout << "end of mutations!!!! break " << std::endl;
                break;
            }

            //node_it --;
            //std::cout << "before missing beg " << missing_start << std::endl;
            //std::cout << "before missing end" << missing_end<< std::endl;
            //new_missing.push_back(std::make_tuple(temp_start, std::get<1>(temp_list[temp_counter])));
            //new_missing_idx += 1;
            //missing_counter += 1;
        }
             
        //if mut is before missing, missing needs to decrement to get closer to mut 
        else if (mut_pos < missing_start){
            //std::cout << "mut is fully before missing" << std::endl;
            //logging.debug('line less than mask')
            //new_missing.push_back(std::make_tuple(temp_start, std::get<1>(temp_list[temp_counter])));
            
            //leaf_counter --;
            //if (leaf_counter != diff_data[leaf->identifier].begin()) {
            if (leaf_counter != missing_data.begin()) {
                leaf_counter --;
            }
            else {
                //need to do something first?
                //std::cout << "end of missing!!!! break " << std::endl;
                break;
            }
            //std::cout << "made it here?" << std::endl;
        
        }
        else {
            //std::cout << "stuck here" << std::endl;
        }
        
        //if only one of them is still going 
        //still need to edit this

        /*
        else if (node_it == node->mutations.begin() || leaf_counter == diff_data[leaf->identifier].begin()) {
            std::cout << "work now???" << std::endl;
            node_it --;
            leaf_counter --;
            //std::cout << "work now???" << std::endl;

            last_iteration == true;
        }
        

        else {
            std::cout << "stuck here" << std::endl;
            std::cout << "node len " << node_len << std::endl;
            std::cout << "leaf len" << leaf_len << std::endl;
            //std::cout << "node " << node_counter << std::endl;
            std::cout << "leaf counter " << leaf_counter->first << " "  << leaf_counter->second << std::endl;
            std::cout << "what is this " << diff_data[leaf->identifier].end()->first <<  " " << diff_data[leaf->identifier].end()->second << std::endl;
            

            //std::cout << "leaf " << leaf_counter << std::endl;

            std::cout << "stuck here" << std::endl;

        } 
        */   
            //std::cout << "mutation " << node_it->get_string() << std::endl;
            //std::cout << "missing " << leaf_counter->first << " " << leaf_counter->second << std::endl;
 
    //}
        
 

   
   // this might be too slow but i have an idea? 
    }
            
}

bool prev_check(std::pair<int, int> prev, std::pair<int, int> line) {
    //std::cout <<   "PREV CHECK BB" << std::endl;
    //std::cout <<   "PREV " << prev.first << " " << prev.second << std::endl;
    //std::cout <<   "line " << line.first << " " << line.second << std::endl;
    if ( (prev.first + prev.second) >= line.first) {
        //std::cout <<   "YES OVERLAP PREV " <<  prev.first << " " << prev.first + prev.second << std::endl;
        //std::cout <<   "YES OVERLAP LINE " <<  line.first << " " << line.first + line.second << std::endl;

        return 1;
    }
    else {
        //std::cout <<   "No OVERLAP? " <<  prev.first + prev.second << " " << line.first << std::endl;
        return 0;
    }



}

void combine_missing(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, std::list<std::pair<int, int>>& missing_data, std::map<std::string, std::map<int, int>>& diff_data) {
    //get missing list lens for iterating 
    auto node_len = diff_data[node->identifier].size();
    auto node_iterator = diff_data[node->identifier].begin();
    auto leaf_len = diff_data[leaf->identifier].size();
    auto leaf_iterator = diff_data[leaf->identifier].begin();
    //std::cout <<   "node len" << node_len << std::endl;
    //std::cout <<   "node_iterator" << node_iterator->first << std::endl;
    //std::cout <<   "leaf len" << leaf_len << std::endl;
    //std::cout <<   "leaf_iterator" << leaf_iterator->first << std::endl;
    //while one list is still going
    std::pair<int, int> prev = {-1, -1};
    while (node_iterator != diff_data[node->identifier].end() || leaf_iterator != diff_data[leaf->identifier].end()) {
            //while two lists are still going 
            if (node_iterator != diff_data[node->identifier].end() && leaf_iterator != diff_data[leaf->identifier].end()) {
                //get positions for both lists
                int node_start = node_iterator->first;
                int node_end = node_start + node_iterator->second;
                int leaf_start = leaf_iterator->first;
                int leaf_end = leaf_start + leaf_iterator->second;
                //std::cout << "current node pos " << node_start << " " << node_end << std::endl;
                //std::cout << "current leaf pos " << leaf_start << " " << leaf_end << std::endl;
                if (node_start >= leaf_start && node_end <= leaf_end) {
                    //std::cout << "node missing is fully inside leaf missing " << std::endl;

                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    if (prev.first == -1) {
                        //missing_data.emplace_back(leaf_start, leaf_end-leaf_start);
                        prev = {leaf_start, leaf_end-leaf_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{leaf_start, leaf_end-leaf_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {leaf_start, leaf_end-leaf_start};

                        }

                    }
                    node_iterator ++;

                }
                else if (node_start <= leaf_start && node_end >= leaf_end) {
                    //std::cout << "leaf missing is fully inside node missing " << std::endl;

                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    if (prev.first == -1) {
                        //missing_data.emplace_back(node_start, node_end-node_start);
                        prev = {node_start, node_end-node_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{node_start, node_end-node_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {node_start, node_end-node_start};

                        }

                    }
                    //missing_data.emplace_back(node_start, node_end-node_start);

                    leaf_iterator ++;

                }
                else if (node_start <= leaf_end && node_end > leaf_end) {
                    //std::cout << "node missing overlaps leaf missing on the right" << std::endl;
                    if (prev.first == -1) {
                        //missing_data.emplace_back(leaf_start, node_end-leaf_start);
                        prev = {leaf_start, node_end-leaf_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>(leaf_start, node_end-leaf_start))) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {leaf_start, node_end-leaf_start};

                        }

                    }
                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    //missing_data.emplace_back(leaf_start, node_end-leaf_start);

                    leaf_iterator ++;

                }
                else if (node_start < leaf_start && node_end >= leaf_start) {
                    //std::cout << "leaf missing overlaps node missing on the right" << std::endl;
                    if (prev.first == -1) {
                        //missing_data.emplace_back(node_start, leaf_end-node_start);
                        prev = {node_start, leaf_end-node_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{node_start, leaf_end-node_start})) {
                           //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {node_start, leaf_end-node_start};

                        }

                    }
                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    //missing_data.emplace_back(node_start, leaf_end-node_start);
                    node_iterator ++;

                }
                else if (node_end <= leaf_start) {
                    //std::cout << "node missing is completely before leaf missing" << std::endl;
                    if (prev.first == -1) {
                        //missing_data.emplace_back(node_start, node_end-node_start);
                        prev = {node_start, node_end-node_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{node_start, node_end-node_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {node_start, node_end-node_start};

                        }

                    }
                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    //missing_data.emplace_back(node_start, node_end-node_start);
                    node_iterator ++;

                }
                else if (leaf_end <= node_start) {
                    //std::cout << "leaf missing is completely before node missing" << std::endl;
                    if (prev.first == -1) {
                       //missing_data.emplace_back(leaf_start, leaf_end-leaf_start);
                        prev = {leaf_start, leaf_end-leaf_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{leaf_start, leaf_end-leaf_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {leaf_start, leaf_end-leaf_start};


                        }

                    }
                    //del_muts.push_back(node->mutations[node_counter]);
                    //node->mutations.erase(node_it);
                    //std::cout << "CURRENTIT " << node_it->get_string() << std::endl;
                    //missing_data.emplace_back(leaf_start, leaf_end-leaf_start);
                    leaf_iterator ++;

                }
                else {
                    //std::cout << "WHAT IS HERE??? both loops" << std::endl;

                }
            }
            else if (node_iterator != diff_data[node->identifier].end()) {
                //std::cout << "node list is still going" << std::endl;
                int node_start = node_iterator->first;
                int node_end = node_start + node_iterator->second;
                if (prev.first == -1) {
                        //missing_data.emplace_back(node_start, node_end-node_start);
                        prev = {node_start, node_end-node_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{node_start, node_end-node_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {node_start, node_end-node_start};


                        }
                    }
                //missing_data.emplace_back(node_start, node_end-node_start);
                
                node_iterator ++;


            } 
            else if (leaf_iterator != diff_data[leaf->identifier].end()) {
                //std::cout << "leaf list is still going" << std::endl;
                int leaf_start = leaf_iterator->first;
                int leaf_end = leaf_start + leaf_iterator->second;
                //std::cout << "current node pos " << node_start << " " << node_end << std::endl;
                //std::cout << "current leaf pos " << leaf_start << " " << leaf_end << std::endl;
                if (prev.first == -1) {
                        //missing_data.emplace_back(leaf_start, leaf_end-leaf_start);
                        prev = {leaf_start, leaf_end-leaf_start};
                    }
                    else {
                        if (prev_check(prev,std::pair<int,int>{leaf_start, leaf_end-leaf_start})) {
                            //std::cout << "PREVE CHECK TRUE " << "there is overlap " << std::endl;
                            //std::cout << "PREVE " << prev.first << " " << prev.second << std::endl;
                            //                            std::cout << "PREVE " << prev.first << " " << prev.second << std::endl;


                        }
                        else {
                            //std::cout << "PREVE CHECK FALSE " << "there is no overlap " << std::endl;
                            missing_data.emplace_back(prev.first, prev.second);
                            prev = {leaf_start, leaf_end-leaf_start};


                        }
                    }
                //missing_data.emplace_back(leaf_start, leaf_end-leaf_start);
                leaf_iterator ++;
            }
            else {
                    std::cout << "WHAT IS HERE??? not both loops" << std::endl;

                }
            //add one to mutation index (there might be more mutations in missing region so dont index down missing)
            //node_counter -= 1;

            //do i need to decrement if i delete node_it? 
            //POSSIBLE ERROR HERE
            
            
        
            /*
            missing_iterator ++;
            data_iterator ++;

    }
    
    
*/
    //add last prev to the end 
    //cout << node->identifier << " " << leaf->identifier << std::endl;

    }
    //last one?
    missing_data.emplace_back(prev.first, prev.second);
    /*
    for (auto i: missing_data) {
        std::cout <<   "missing_data" << i.first << " " << i.second << std::endl;
        }
        */
//6/23 make sure combined missing data is complete, rework node comp to reduce number of node comparisons 
}


void traverseToMRCA(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data) {
    //
    auto current_node = node;
    while (current_node->identifier != mrca->identifier) {
        
            //std::lock_guard<std::mutex> guard(mtx); // Ensure thread-safe operation
            nodeComp(current_node, leaf, diff_data, missing_data);
        
        //std::cout << "current_node: " << current_node->identifier << std::endl;
        current_node = current_node->parent;
        //std::cout << "current_node: " << current_node->identifier << std::endl;
    }
}

// Function to traverse from leaf's parent to MRCA
void traverseFromLeafParentToMRCA(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::list<std::pair<int, int>>& missing_data) {
    auto current_node = leaf;
    //std::cout << "OTHER SIDE NODE: " << current_node->identifier << std::endl;
    while (current_node->identifier != mrca->identifier) {
        
        //std::lock_guard<std::mutex> guard(mtx); // Ensure thread-safe operation
        nodeComp(current_node, leaf, diff_data, missing_data);
        //std::cout << "missing data set after: " << missing_data.size() << std::endl;
        
        //std::cout << "current_node: " << current_node->identifier << std::endl;
        current_node = current_node->parent;
        //std::cout << "current_node: " << current_node->identifier << std::endl;
        //std::cout << "mrca: " << mrca->identifier << std::endl;
    }
}


//this switches leaf and node and i dont like it 
void getDistance(Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data, std::map<std::string, std::set<std::string>>& comparisons) {
    //std::map<int, int> missing = diff_data[node->identifier];
    
    
    //node is the neighbor, we are starting with deleting mutations from the neighbor
    //then we go up the path towards the mrca, deleting mutations that are in missing regions of leaf
    auto current_node = node; 

    //add leaf to comparisons if its not already there
    //which it shouldnt be? 
    //std::cout <<   "COMPARISONS BEFORE" << comparisons.size() << std::endl;


    if ( comparisons[leaf->identifier].find(current_node->identifier) == comparisons[leaf->identifier].end() ) {
        //std::cout <<   "NOT COMPARED BEFORE, Comparing now" << comparisons.size() << std::endl;
        std::list<std::pair<int, int>> missing_data;
        //std::mutex mtx; // Mutex for thread-safe operations

        comparisons[leaf->identifier].insert(node->identifier);
        comparisons[node->identifier].insert(leaf->identifier);
        //std::cout <<   "MISSING BEFORE" << missing_data.size() << std::endl;
        combine_missing(node, leaf, missing_data, diff_data);
        //std::cout <<   "MISSING AFTER" << missing_data.size() << std::endl;
        //std::cout <<   "LEAF MISSINg " << diff_data[leaf->identifier].size() << std::endl;
        //std::cout <<   "NODE MISSINg " << diff_data[node->identifier].size() << std::endl;

        //std::cout <<   "MADE IT TO GET DISTANCE" << std::endl;
        //std::thread t1(traverseToMRCA, leaf, node, mrca, std::ref(diff_data), std::ref(missing_data));
        //std::thread t2(traverseFromLeafParentToMRCA, leaf, node, mrca, std::ref(diff_data), std::ref(missing_data));
        //std::future<int> f = std::async(std::launch::async, traverseToMRCA, leaf, node, mrca, std::ref(diff_data), std::ref(missing_data), std::ref(mtx));
        //auto future1 = std::async(std::launch::async, traverseToMRCA, leaf, node, mrca, std::ref(diff_data), std::ref(missing_data));
        //auto future2 = std::async(std::launch::async, traverseFromLeafParentToMRCA, leaf, node, mrca, std::ref(diff_data), std::ref(missing_data));

        //future1.get();
        //future2.get();
        //int result1 = my_func(arg1, arg2, arg3);
        //int result2 = f.get(); 
        // Wait for both threads to complete
        //future1.join();
        //future2.join();
        
        //for the path on the opposite side of the mrca from the leaf
        
        while (current_node->identifier != mrca->identifier ) {
            //if (visited.find(current_node->identifier) != visited.end()){
            //std::cout << "missing data set before " << missing_data.size() << endl;
            

            nodeComp(current_node, leaf, diff_data, missing_data);
            //visited.insert(current_node->identifier);
            //std::cout << "missing data set after " << missing_data.size() << endl;

            //}
            //else {
            //    std::cout << "not SKIPPED YAY " << endl;
            //    nodeComp(current_node, leaf, diff_data,missing_data);


            //}
            //std::cout << "current_node" << current_node->identifier << endl;
            current_node = current_node->parent;
            //std::cout << "current_node" << current_node->identifier << endl;

        } 
        
        
        //when we hit the mrca, we move to the other side, the leaf side 
        //std::cout <<   "OTHER SIDE" << std::endl;
        //we dont need to include the leaf bc those mutations already dont exist
        
        //this might need to be leaf!!!!!!!!!
        current_node = leaf->parent; 
        //std::cout << "OTHER SIDE NODE" << current_node->identifier << std::endl;
        while (current_node->identifier != mrca->identifier ) {
            //std::cout << "OTHER SIDE WHILE LOOP" << std::endl;
            //if (visited.find(current_node->identifier) != visited.end()){
            //std::cout << "missing data set before " << missing_data.size() << endl;
            nodeComp(current_node, leaf, diff_data,missing_data);
            //visited.insert(current_node->identifier);
            //std::cout << "missing data set after " << missing_data.size() << endl;

            //}
            //else {
            //    std::cout << "Not SKIPPED YAY " << missing_data.size() << endl;
            //    nodeComp(current_node, leaf, diff_data,missing_data);

            //}
            //nodeComp(current_node, leaf, diff_data);
            //combine_missing(current_node, missing_data, diff_data);

            //std::cout << "current_node" << current_node->identifier << endl;
            current_node = current_node->parent;
            //std::cout << "current_node" << current_node->identifier << endl;
            //std::cout << "mrca" << mrca->identifier << endl;
            //if (current_node->identifier != mrca->identifier ) {
            //}
        } 
        
        //combine_missing(mrca, missing_data, diff_data);
        //do the mrca last
        nodeComp(mrca, leaf, diff_data,missing_data);
        //change comparisons here
    }
    else {
        //std::cout <<   "SKIPPED BC RENDUNDY" << leaf->identifier << node->identifier << std::endl;

    }
    //std::cout <<   "COMPARISONS AFTER" << comparisons.size() << std::endl;

    //this code assumes that all mutations entered are snps, will not work if mutations are longer than 1bp (which is currently usher's capability)
    
}


//new getDistance here
/*
void getDistance(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* leaf, Mutation_Annotated_Tree::Node* mrca, std::map<std::string, std::map<int, int>>& diff_data) {
    //std::map<int, int> missing = diff_data[node->identifier];
    std::cout <<   "MADE IT TO GET DISTANCE" << std::endl;
    
    //node is the starting leaf, we are deleting mutations from the starting leaf
    auto current_node = node; 
    
    std::cout <<   "MUTATIONS " << node->identifier << std::endl;
    for (const auto& mut: node->mutations) {
            std::cout << "pos" << mut.get_string() << endl;
            std::cout << "pos" << mut.mut_nuc << endl;
            std::cout << "pos" << mut.ref_nuc << endl;
            std::cout << "pos" << mut.is_masked() << endl;
            std::cout << "pos" << mut.is_missing << endl;

    }
    
    while (current_node->identifier != mrca->identifier ) {
        nodeComp(current_node, leaf, diff_data);
        std::cout << "current_node" << current_node->identifier << endl;
        current_node = current_node->parent;
        std::cout << "current_node" << current_node->identifier << endl;
        std::cout << "mrca " << mrca->identifier << endl;
        //if (current_node->identifier != mrca->identifier ) {
        //}
    } 
    std::cout <<   "OTHER SIDE" << std::endl;

    current_node = leaf->parent; 
    std::cout << "OTHER SIDE NODE" << current_node->identifier << std::endl;
    while (current_node->identifier != mrca->identifier ) {
        std::cout << "OTHER SIDE WHILE LOOP" << std::endl;
        nodeComp(current_node, leaf, diff_data);
        std::cout << "current_node" << current_node->identifier << endl;
        current_node = current_node->parent;
        std::cout << "current_node" << current_node->identifier << endl;
        std::cout << "mrca" << mrca->identifier << endl;
        //if (current_node->identifier != mrca->identifier ) {
        //}
    } 
    nodeComp(mrca, leaf, diff_data);

    

    //this code assumes that all mutations entered are snps, will not work if mutations are longer than 1bp (which is currently usher's capability)
    
}

      
*/




//different version of masking, not subtree, just max distance 
void localMask (uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file, std::string filename, uint32_t num_threads) {
    //collect all missing data for each leaf
    std::map<std::string, std::map<int, int>> diff_data = readDiff(diff_file);
    //idk if i need this
    std::set<std::string> leaves;
    //retrieve leaves 
    auto all_leaves = T.get_leaves();
    std::map<std::string, std::set<std::string>> comparisons;

    //i dont think parallel should happen here
    //#pragma omp parallel for num_threads(num_threads)
    for (auto l: all_leaves) {

        //turn this into a way to avoid recomparing a node to a set of missing data
        //std::set<std::string> visited;
        //string of node id maybe delete
        std::string samp = l->identifier;
        //retrieve leaf branch length, do i need this? 
        int bl = l->branch_length;
        //set current node == leaf , do i need this? 
        MAT::Node* current_node = l;
        //error check,delete eventually
        std::cout << "sample: " << samp << std::endl;

        //need? 
        //leaves.insert(l->identifier);
        //note that my leaf is going to be masking all nodes on the path between it and it's closest neighbors
        //i.e. i care about the missing data from my leaf, NOT the mutations        
        
        //determine if branchlen immediately disqualifies leaf from masking
        if (l->branch_length < max_snp_distance && diff_data.find(samp) != diff_data.end() ) {

        //get nearest neighbors of leaf
        std::pair<std::vector<std::string>, size_t> neighbors = get_closest_samples(&T, l->identifier, true, max_snp_distance);            
        
        //iterate through neibours, is there a better way to do this? 
        //can i parallelize this? 
        std::cout <<   "NEIGBORS of " << l->identifier << std::endl; 
        
        //const int num_threads = 4; // You can adjust this based on your system
        //#pragma omp parallel for num_threads(num_threads)
        // iterate through nearest neighbors, figure out LCA between them
        for (const auto& neigh : neighbors.first) {
            //identifier for neighbor
            std::cout << neigh << std::endl;
            //get mrca for neighbor and leaf so we can find path between
            auto mrca = MAT::LCA(T, l->identifier, neigh);

            //send leaf node and neighbor node to next function
            //i thnk i need to add visited to getDistance function
            //i need leaf and node to be treated as follows: i care about missing data for leaf,
            //i care about mutation data for node in get distance, leaf comes first and node comes second 
            //regarless of how they are subsequently named 
            getDistance(l, T.get_node(neigh), mrca, diff_data, comparisons);
            }
        }
        else {
            std::cout <<   "SAMPLE doesnt exist " << l->identifier << std::endl; 

        }
        //}
        }

        
        //find ancestor node
        //need to make a catch so that max_snp-_dist isnt 0 when snp distance is > 0
        // allows for checking all nearby neibors as long as comman anc is less than max branch distance away
        
        //this determines how many nodes deep we can go before maxing branch length 
        
        /*
        while (bl < max_snp_distance) {
            //add current branch to total branch length
            //int current_branch = current_node->branch_length; 
            
            std::cout << "current node" << current_node->identifier << std::endl;
            std::cout << "branch len" << bl << std::endl;
            //std::cout << "current branch" << current_branch << std::endl;
            

            //make this catch things with root ancestor node 
            if (current_node->parent == NULL) {
                std::cout <<   "ROOOOOOTTTTTT"  << std::endl;
                break;
            }
            
            //update current node to parent if branch length is not met 
            
            //else if ((bl + current_branch ) < snp_distance) 
            else {
                //bl += current_branch;
                //std::cout << "is this the right branch?" << current_node->parent->branch_length << std::endl;
                std::cout <<  "leaf node" << l->identifier << std::endl;
                std::cout <<  "old current node" << current_node->identifier << std::endl;
                std::cout << "old branch len" << bl << std::endl;
                //std::cout << "old current branch" << current_branch << std::endl;
                
                
                if (bl + current_node->parent->branch_length < max_snp_distance) {
                    std::cout <<  "old bl " << bl << std::endl;
                    current_node = current_node->parent;
                    
                    std::cout <<  "new current node " << current_node->identifier << std::endl;
                    std::cout <<  "current node branch len" << current_node->branch_length << std::endl;
                    
                    std::cout <<  "leaf node (this shouldn't change)" << l->identifier << std::endl;
                    std::cout << "new current node" << current_node->identifier << std::endl;
                    std::cout << "new branch len (should stay the same)" << bl << std::endl;
                    //current branch shouldn't change 
                    std::cout << "new current branch (this should probably change)" << current_node->branch_length << std::endl;
                    std::cout <<  "dfs call leaf node" << l->identifier << std::endl;
                    std::cout <<  "dfs call branch length" << bl << std::endl;
                    std::cout <<  "dfs call current node " << current_node -> identifier << std::endl;
                    
                    
                    //auto desc = T.get_leaves_ids(current_node->identifier);
                    //dfs(l, bl, current_node, diff_data, max_snp_distance, visited, leaves);
                    std::pair<std::vector<std::string>, size_t> neighbors = get_closest_samples(&T, l->identifier, false, max_snp_distance);
                    
                    std::cout <<   "NEIGBORS of " << l->identifier << std::endl;    
                    for (const auto& leaf : neighbors.first) {
                        std::cout << leaf << std::endl;
                    }
                    bl += current_node->branch_length;
                    //if more nodes can be checked, update bl, if not, exit while loop
                    //if (bl + current_node->branch_length < snp_distance) {
                }
                else {
                    current_node = current_node->parent;
                    break;
                }
            }   
        }
        std::cout <<  "dfs call leaf node" << l->identifier << std::endl;
        std::cout <<  "dfs call branch length" << bl << std::endl;
        std::cout <<  "dfs call current node " << current_node -> identifier << std::endl;

        //dfs(l, bl, current_node, diff_data, max_snp_distance, visited, leaves);
        
        }
*/
MAT::save_mutation_annotated_tree(T, filename);
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
        if (target->par_nuc != query->par_nuc) {
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
    size_t total_masked = 0;
    timer.Start();
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
            // The expected common case is to not match any mutations and have nothing to remove.
            std::vector<MAT::Mutation> muts_to_remove;
            for (auto& mut: n->mutations) {
                if (match_mutations(mutobj, &mut)) {
                    instances_masked++;
                    muts_to_remove.push_back(mut);
                }
            }
            for (auto mut: muts_to_remove) {
                auto iter = std::find(n->mutations.begin(), n->mutations.end(), mut);
                n->mutations.erase(iter);
            }
        }
        total_masked += instances_masked;
    }
    fprintf(stderr, "Completed in %ld msec \n", timer.Stop());
    infile.close();
    fprintf(stderr, "Masked a total of %lu mutations.  Collapsing tree...\n", total_masked);
    timer.Start();
    T->collapse_tree();
    fprintf(stderr, "Completed in %ld msec \n", timer.Stop());
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

std::unordered_set<std::string> mutation_set_from_node(MAT::Tree* T, MAT::Node* node, bool include_node, bool include_ancestors) {
    std::unordered_set<std::string> mutations;
    if (include_ancestors) {
        for (auto an: T->rsearch(node->identifier, include_node)) {
            for (auto mut: an->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                MAT::Mutation mut_opposite = mut.copy();
                //we check for whether this mutation is going to negate with something in the set
                //by identifying its opposite and checking whether the opposite is already present on the traversal.
                mut_opposite.par_nuc = mut.mut_nuc;
                mut_opposite.mut_nuc = mut.par_nuc;
                auto cml = mutations.find(mut_opposite.get_string());
                if (cml != mutations.end()) {
                    mutations.erase(cml);
                } else {
                    mutations.insert(mut.get_string());
                }
            }
        }
    } else if (include_node) {
        for (auto mut: node->mutations) {
            if (mut.is_masked()) {
                continue;
            }
            mutations.insert(mut.get_string());
        }
    } else {
        fprintf(stderr, "ERROR: mutation_set_from_node: at least one of include_node and include_ancestors should be true.\n");
        exit(1);
    }
    return mutations;
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
        std::unordered_set<std::string> curr_mutations = mutation_set_from_node(T, mn, false, true);
        std::unordered_set<std::string> new_mutations = mutation_set_from_node(T, np, true, true);
        if (curr_mutations == new_mutations) {
            //we can now proceed with the move.
            T->move_node(mn->identifier, np->identifier);
            fprintf(stderr, "Move of node %s to node %s successful.\n", words[0].c_str(), words[1].c_str());
        } else {
            // Not quite the same; figure out whether we need to add a node under new parent.
            fprintf(stderr, "The current (%s) and new (%s) node paths do not involve the same set of mutations.\n",
                    mn->identifier.c_str(), np->identifier.c_str());

            std::unordered_set<std::string> extra_mutations;
            size_t curr_in_new_count = 0;
            for (auto mut: curr_mutations) {
                if (new_mutations.find(mut) != new_mutations.end()) {
                    curr_in_new_count++;
                } else {
                    extra_mutations.insert(mut);
                }
            }
            if (extra_mutations.size() == 0 || curr_in_new_count != new_mutations.size()) {
              fprintf(stderr, "ERROR: the new parent (%s) has mutations not found in the current node (%s); %ld in common, %ld in new\n",
                      np->identifier.c_str(), mn->identifier.c_str(), curr_in_new_count, new_mutations.size());
              exit(1);
            }
            // Look for a child of np that already has extra_mutations.  If there is such a child
            // then move mn to that child.  Otherwise add those mutations to mn and move it to np.
            MAT::Node *child_with_muts = NULL;
            for (auto child: np->children) {
              std::unordered_set<std::string> mut_set = mutation_set_from_node(T, child, true, false);
                if (mut_set == extra_mutations) {
                    child_with_muts = child;
                    fprintf(stderr, "Found child with extra_mutations: %s\n", child->identifier.c_str());
                    break;
                }
            }
            if (child_with_muts != NULL) {
                T->move_node(mn->identifier, child_with_muts->identifier);
            } else {
                // Preserve chronological order expected by add_mutation by adding mn's mutations
                // after extra_mutations instead of vice versa.
                std::vector<MAT::Mutation>mn_mutations;
                for (auto mut: mn->mutations) {
                    mn_mutations.push_back(mut);
                }
                mn->mutations.clear();
                for (auto mut: extra_mutations) {
                    mn->add_mutation(*(MAT::mutation_from_string(mut)));
                }
                for (auto mut: mn_mutations) {
                    mn->add_mutation(mut);
                }
                T->move_node(mn->identifier, np->identifier);
            }
        }
    }
    fprintf(stderr, "All requested moves complete.\n");
}


