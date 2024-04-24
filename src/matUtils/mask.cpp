#include "mask.hpp"
#include <random>
#include <algorithm>

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
        fprintf(stderr, "made it to here");
        localMask(max_snp_distance, T, diff_file);
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
    fprintf(stderr, "made it to readDiff");
    std::map<std::string, std::map<int, int>> data;


    try {
        std::ifstream file(diff_file);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + diff_file);
        }

        std::string line;
        std::string current_sample;
        while (std::getline(file, line)) {
            // Parse the line and store data in the map (example)
            // For demonstration, assume each line contains key-value pairs separated by '='
            std::vector<std::string> substrings;
            size_t startPos = 0;
            size_t endPos;

            if (line[0] == '>') {
                current_sample = line.erase(0,1);
                auto it = data.find(current_sample);
                if (it == data.end()) {
                    // Inner map doesn't exist for the current sample
                    // Create a new inner map for the current sample
                    data[current_sample] = std::map<int, int>();
                    //std::cout << "Inner map created for sample: " << current_sample << std::endl;
                } else {
                    // Inner map already exists for the current sample
                    // throw an error 
                            throw std::runtime_error("Duplicate samples detected, inspect diff file for sample: " + current_sample);
                }
                //std::cout << "sample" << current_sample << endl;
                // add error checking to this
                //data[current_sample] = std::map<int, int>(); 
                
            }
            // Find the position of the first occurrence of the separator
            else if (line[0] == '-'){
                
                while ((endPos = line.find('\t', startPos)) != std::string::npos) {
                    // Extract the substring between startPos and endPos and add it to the substrings vector
                    //std::string cell = line.substr(startPos, endPos - startPos)
                    
                    substrings.push_back(line.substr(startPos, endPos - startPos));
                    //std::cout << line.substr(startPos, endPos - startPos) << endl;
                    //std::cout << "Vector contents: " << endl;
                    


                    // Update startPos to the position after the separator
                    startPos = endPos + 1;
                    //std::cout << "start " << startPos << " " << "end: " << endPos << " " << "line " << line << std::endl;
                    //std::cout << "line " << line << std::endl;
                    
                }
            

            // Extract the substring after the last occurrence of the separator
            substrings.push_back(line.substr(startPos));
            //std::cout << "Vector contents: " << endl;
            //for (const auto& word : substrings) {
            //        std::cout << "sample" << current_sample << endl;
            //    }

            //std::cout << "subs" << substrings[1] << endl;
            //values = std::vector<string> {substrings[1], substrings[2]};
            int position = std::stoi(substrings[1]);
            //if (position > 4500000) {
            //    std::cout << "current sample" << current_sample << position << endl;

                
            //}
            int length = std::stoi(substrings[2]);
            data[current_sample][position] = length;
                
            /*
            size_t pos = line.find('\t');
            if (pos != std::string::npos) {
                std::string key = line.substr(0, pos);
                std::string position = line.substr(pos + 1, pos +2);
                std::string length = line.substr(pos + 3);
                if (key == "-") {
                    std::cout << "Key: " << key << ", Value: " << position << 'len' << length << std::endl;
                }

                //fprintf(stderr, "oh shit this works %s, %s.\n", key, value);
                //data[key] = value;
            }*/
                
            }
        }
        //std::cout << "Vector contents: " << endl;
        //for (const auto& word : substrings) {
        //    std::cout << word << endl;
        //    }

        file.close();
    } catch (const std::exception& e) {
        std::cerr << "Exception caught while reading file: " << e.what() << std::endl;
        throw;  // Re-throw the exception to be handled where the function is called
    }
    
    /*
    for (const auto& word : data) {
            std::cout << "sample" << word.first << endl;
            for (const auto& pos : word.second) {
                std::cout << "pos" << pos.first << '\t' << pos.second << endl;
            }
            

    }
    */
    
    return data;
}

void dfs(MAT::Node* l, int bl, Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data, uint32_t snp_distance) {
    // Initialize a set to keep track of visited nodes
    //visited <- empty set
    //MAT::Node* anc_node = node;
    int distance_from_l = bl;
    std::set<std::string> visited;
    //std::vector<std::tuple<int, int>> missing;
    //std::vector<std::tuple<int, int>> new_missing;
    //fprintf(stderr, "dfs.\n");
    std::cout <<   "dfs current node " << node->identifier << std::endl;
    // Call the recursive DFS function
    /*
    for (auto i: diff_data[node->identifier]) {
        //subtreeMask(node);
        std::cout <<   "diff data" << i.first << std::endl;
        std::cout <<   "diff data" << i.second << std::endl;
    }
    */
    dfsUtil(l, distance_from_l, node, diff_data, visited, snp_distance);
}

// Define a utility function for DFS traversal
void dfsUtil(MAT::Node* l, int distance_from_l, Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data, std::set<std::string>& visited, uint32_t snp_distance) {
    // Mark the current node as visited
    
    visited.insert(node->identifier);

    //std::cout << "Length of the set: " << visited.size() << std::endl;
    //std::cout << "ancestor node: " << anc_node->identifier << std::endl;
    std::cout << "dfsutil current node " << node->identifier << std::endl;
    //fprintf(stderr, "dfsutil.\n");


    //process node here
    // Perform any operation on the current node
    // (e.g., print its value, process it, etc.)
   //processNode(node, anc_node, new_missing, missing, diff_data);
    
    // Explore all adjacent nodes of the current node
    for (auto c: node->children) {
        //subtreeMask(node);
        if (c-> identifier != l->identifier) {
            std::cout <<   "made it to here" << c->identifier << std::endl;

            std::cout <<   "child of node " << node->identifier << ": " << c->identifier << std::endl;
            std::cout <<   "child branch length: " << c->branch_length << std::endl;
            std::cout <<   "old branch length: " << distance_from_l << std::endl;
            int new_bl = distance_from_l + c->branch_length;
            //bl+=c->branch_length;
            std::cout <<   "new branch length: " << distance_from_l << std::endl;
            
            std::cout <<   c->identifier << "total distance from l " << l-> identifier << ": " << new_bl << std::endl;
            
            if (visited.find(c->identifier) == visited.end() && c->is_leaf() == false ) {
                std::cout <<   "unvisited: " << c->identifier << std::endl;    
                for (std::string element : visited) {
                    std::cout << element << ", ";
                }
                std::cout << std::endl;

                std::cout << "new node" << c->identifier << std::endl;
                dfsUtil(l, new_bl, c, diff_data, visited, snp_distance);
            }
        }
    /*
    for each neighbor in graph[currentNode]:
        
        if neighbor not in visited:
            // Recursive call to DFSUtil for unvisited neighbors
            dfsUtil(graph, neighbor, visited)
            */
    }
}

/*
    std::ifstream file(diff_file);
    
    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1; // Return an error code
    }

    // Read the file line by line
    std::string line;
    while (std::getline(file, line)) {
        // Process each line (in this example, simply print it)
        std::cout << line << std::endl;
    }
    
    // Close the file
    file.close();
    
}
*/

//old version, subtree focused 
/*
void dfs(Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data) {
    // Initialize a set to keep track of visited nodes
    //visited <- empty set
    MAT::Node* anc_node = node;
    std::set<std::string> visited;
    std::vector<std::tuple<int, int>> missing;
    std::vector<std::tuple<int, int>> new_missing;
    fprintf(stderr, "dfs.\n");
    std::cout <<   "node" << node->identifier << std::endl;
    // Call the recursive DFS function
    
    //for (auto i: diff_data[node->identifier]) {
        //subtreeMask(node);
    //    std::cout <<   "diff data" << i.first << std::endl;
    //    std::cout <<   "diff data" << i.second << std::endl;
    }
    
    dfsUtil(node, anc_node, diff_data, new_missing, missing, visited);
}

// Define a utility function for DFS traversal
void dfsUtil(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* anc_node, std::map<std::string, std::map<int, int>>& diff_data, std::vector<std::tuple<int, int>>& new_missing,  std::vector<std::tuple<int, int>>& missing, std::set<std::string>& visited) {
    // Mark the current node as visited
    
    visited.insert(node->identifier);

    std::cout << "Length of the set: " << visited.size() << std::endl;
    std::cout << "ancestor node: " << anc_node->identifier << std::endl;
    std::cout << "current node " << node->identifier << std::endl;
    fprintf(stderr, "dfsutil.\n");


    //process node here
    // Perform any operation on the current node
    // (e.g., print its value, process it, etc.)
    processNode(node, anc_node, new_missing, missing, diff_data);
    
    // Explore all adjacent nodes of the current node
    for (auto c: node->children) {
        //subtreeMask(node);
        std::cout <<   "child" << c->identifier << std::endl;
        if (visited.find(c->identifier) == visited.end()) {
            dfsUtil(c, anc_node, diff_data, new_missing, missing, visited);
        }
    
    //for each neighbor in graph[currentNode]:
        
        //if neighbor not in visited:
            // Recursive call to DFSUtil for unvisited neighbors
        //    dfsUtil(graph, neighbor, visited)
            
    }
}
*/
// Define a function to process a visited node
// collect all MISSING INFORMATION from a subtree
void processNode(Mutation_Annotated_Tree::Node* node, Mutation_Annotated_Tree::Node* anc_node, std::vector<std::tuple<int, int>>& new_missing, std::vector<std::tuple<int, int>>& missing, std::map<std::string, std::map<int, int>>& diff_data) {
    //int missing_len = diff_data[anc_node->identifier].size();
    //int muts_len = node->mutations.size();
    //auto missing_counter = diff_data[anc_node->identifier].begin();
    int temp_counter = 0;
    int missing_counter = 0;
    int new_missing_idx = 0;
    std::vector<std::tuple<int, int>> temp_list;

    //std::cout << "ancestor node: " << anc_node->identifier << std::endl;
    //std::cout << "current node " << node->identifier << std::endl;
    //std::cout << "missing size " << missing.size() << std::endl;

    if (missing.size() == 0){
        for (auto i: diff_data[node->identifier]) {
            missing.push_back(std::make_tuple(i.first, i.second));
        }
    }
    
    else {
        for (auto i: diff_data[node->identifier]) {
            temp_list.push_back(std::make_tuple(i.first, i.second));
        
        }
        std::cout << "temp size " << temp_list.size() << std::endl;

        while (temp_counter < temp_list.size() or missing_counter < missing.size()){
            //no new missing to add 
            if (temp_list.size() == 0) {
                break;
            }
            //iterate through both lists still 
            if (temp_counter < temp_list.size() and missing_counter < missing.size()){
                int temp_start = std::get<0>(temp_list[temp_counter]);
                int temp_end = temp_start + std::get<1>(temp_list[temp_counter]);
                int missing_start = std::get<0>(missing[missing_counter]);
                int missing_end = missing_start + std::get<1>(missing[missing_counter]);
                std::cout << "temp start " << temp_start << std::endl;
                std::cout << "temp end " << temp_end << std::endl;
                std::cout << "temp len " << std::get<1>(temp_list[temp_counter]) << std::endl;
                std::cout << "missing start " << missing_start << std::endl;
                std::cout << "missing end " << missing_end << std::endl;
                //temp and missing fully overlap with temp inside
                // missing doesnt need to be modified bc temp is included
                if (temp_start >= missing_start and temp_end <= missing_end) {
                    // may need to make this fancier later
                    std::cout << "temp is fully inside missing " << missing_end << std::endl;
                    new_missing.push_back(missing[missing_counter]);
                    new_missing_idx += 1;
                    temp_counter += 1;
                    missing_counter += 1;
                }
                //full overlap of temp and missing with missing inside
                // missing needs to be modified with temp start and temp end
                else if (temp_start <= missing_start and temp_end >= missing_end) {
                    std::cout << "temp start " << temp_start << std::endl;
                    std::cout << "temp end " << temp_end << std::endl;
                    // might have to make this fancier later 
                    std::cout << "missing is fully inside temp " << std::endl;
                    std::cout << "before missing beg " << missing_start << std::endl;
                    std::cout << "before missing end" << missing_end<< std::endl;
                    new_missing.push_back(std::make_tuple(temp_start, std::get<1>(temp_list[temp_counter])));
                    new_missing_idx += 1;
                    temp_counter += 1;
                    missing_counter += 1;
                    std::cout << "after missing beg " << std::get<0>(new_missing[-1]) << std::endl;
                    std::cout << "after  missing end" << std::get<1>(new_missing[-1]) << std::endl;
                }
                else if (temp_start < missing_start){
                    //logging.debug('line less than mask')
                    new_missing.push_back(std::make_tuple(temp_start, std::get<1>(temp_list[temp_counter])));
                    temp_counter += 1;

                }

                
                
            }

            // this part is temporary to prevent inf loop
            std::cout << "mut counter " << temp_counter << std::endl;
            //std::cout << "temp size " << temp_list.size() << std::endl;
            std::cout << "missing counter " << missing_counter << std::endl;
            //std::cout << "missing size " << missing.size() << std::endl;
            //std::cout << "missing: " << std::get<0>(missing[missing_counter]) << ", " << std::get<1>(missing[missing_counter]) << std::endl;
            //std::cout << "temp: " << std::get<0>(temp_list[temp_counter]) << ", " << std::get<1>(temp_list[temp_counter]) << std::endl;
            //temp_counter += 1;
            //missing_counter += 1;
            if (temp_counter == 5){
                break;
                //int temp_start = std::get<0>(temp_list[temp_counter]);
            }
            





        }
        /*
        for (auto i: diff_data[node->identifier]) {
            //subtreeMask(node);
            std::cout <<   "diff data" << i.first << std::endl;
            std::cout <<   "diff data" << i.second << std::endl;
        }
        */
    }
        
    


    /*
    //std::cout << "" << muts_counter << std::endl;
    //std::cout << "muts" << missing_counter->first << std::endl;
    
    if (anc_node == node) {
        std::cout << "this causes problems" << muts_counter << std::endl;
    } 
    
    //else {
        while (muts_counter <= muts_len or missing_counter != diff_data[anc_node->identifier].end()) {
            //std::cout << "type" << typeid(diff_data[anc_node->identifier]).name() << std::endl;
            muts_counter += 1;
            ++missing_counter;
            std::cout << "while: muts" << muts_counter << std::endl;
            std::cout << "while: missing" << missing_counter->first << std::endl;
            //break;
        //std::cout << "muts" << m << std::endl;
        //break;
    // Body of the while loop
    // Statements to be executed repeatedly as long as the condition is true
        }
    //std::cout << "i am outside the while loop"  << std::endl;
    //}
    /*
    for (const auto& mut: node->mutations) {
            std::cout << "pos" << mut.get_string() << endl;
            //std::cout << "len" << pos.second << endl;
    

    }
    
    
    for (const auto& pos : diff_data[node->identifier]) {
            std::cout << "pos" << pos.first << endl;
            std::cout << "len" << pos.second << endl;

    } 
    

    //print(node)  // Or perform any other operation on the node
*/
}

/*
void dfsUtil(Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data) {
    //visited.insert(node);
    cout << node << " "; // Or you can do any other operation you want with the current node

    for (int child : graph[node]) {
        if (visited.find(neighbor) == visited.end()) {
            dfsUtil(graph, neighbor, visited);
        }
    }
}


void dfs(Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data) {
    unordered_set<int> visited;
    fprintf(stderr, "oh shit this works.\n");
    dfsUtil(node, diff_data);
}

void subtreeMask (Mutation_Annotated_Tree::Node* node, std::map<std::string, std::map<int, int>>& diff_data) {
    fprintf(stderr, "oh shit this works.\n");
    
    //vec.push_back(node);
    for (auto c: node->children) {
        //subtreeMask(node);
        std::cout <<   "child" << c->identifier << std::endl;
    }
    //node->dfs_end_idx=vec.size();
}
*/

//different version of masking, not subtree, just max distance 
void localMask (uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file) {
    std::map<std::string, std::map<int, int>> diff_data = readDiff(diff_file);
    auto all_leaves = T.get_leaves();
    for (auto l: all_leaves) {
        //std::cout << "Data type of l: " << typeid(l).name() << std::endl;
        std::string samp = l->identifier;
        int bl = l->branch_length;
        MAT::Node* current_node = l->parent;
        //MAT::Node* parent = l
        std::cout << "sample: " << samp << std::endl;
        //std::cout <<   "made it here - lk - current node init " << current_node->identifier  << std::endl;
        //std::cout << " branch len " << l -> branch_length << std::endl;
        
        //skip samples coming from root
        /*
        if ((*current_node->parent).is_root()) {
            std::cout <<   "skipped" << current_node->identifier << std::endl;
            continue;
        }
        */

        //find ancestor node
        //need to make a catch so that max_snp-_dist isnt 0 when snp distance is > 0
        // allows for checking all nearby neibors as long as comman anc is less than max branch distance away
        
        //this determines how many nodes deep we can go before maxing branch length 
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
                
                
                if (bl + current_node->branch_length > max_snp_distance) {
                    std::cout <<  "old bl " << bl << std::endl;
                    bl += current_node->branch_length;
                    current_node = current_node->parent;
                    
                    std::cout <<  "new current node " << current_node->identifier << std::endl;
                    std::cout <<  "current node branch len" << current_node->branch_length << std::endl;
                    
                    std::cout <<  "leaf node (this shouldn't change)" << l->identifier << std::endl;
                    std::cout << "new current node" << current_node->identifier << std::endl;
                    std::cout << "new branch len (should stay the same)" << bl << std::endl;
                    //current branch shouldn't change 
                    std::cout << "new current branch (this should probably change)" << current_node->branch_length << std::endl;
                    
                    //if more nodes can be checked, update bl, if not, exit while loop
                    //if (bl + current_node->branch_length < snp_distance) {
                    
                    
                }

                else {
                    break;
                }
                //}
                //else {
                //    break;
                //}


                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
            }
            /*
            else {
                //std::cout << " done " << l->identifier << bl << std::endl;
                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
                break;
            }
            */
            
            //std::cout <<   "made it here - lk"  << std::endl;
            //break;

            
        }
        std::cout <<  "dfs call leaf node" << l->identifier << std::endl;
        std::cout <<  "dfs call branch length" << bl << std::endl;
        std::cout <<  "dfs call current node " << current_node -> identifier << std::endl;
        
        dfs(l, bl, current_node, diff_data, max_snp_distance);
        //std::cout << " done " << l->identifier << bl << std::endl;
        //std::cout << " branch len " << bl << std::endl;
        //std::cout << " new node" << current_node->identifier << std::endl;
        
        //add a warning for this?
        
        /*
        if ((current_node->children).size() >20 ) {
            std::cout <<  "LOTS" << (current_node->children).size() << std::endl;
        }
        */
        
        //need a function here to collect all missing data, 
        //dont pass diff data to dfs, pass a map of all missing data in the subtree

        //need error checking to confirm this is the ancestor node we want 
        //subtreeMask(current_node, diff_data);
        
        //old version: dfs done after anc node identified 
        //dfs(current_node, diff_data);

        //std::cout << " l: " << samp << std::endl;
        //for (const auto& word : diff_data[samp]) {
        //    std::cout << "pos " << word.first << " " << "len " << word.second << endl;
            //for (const auto& pos : word.second) {
            //    std::cout << "pos" << pos.first << '\t' << pos.second << endl;
        //    }
        //break;
        }
}


/*
void localMask (uint32_t snp_distance, uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file) {
    std::map<std::string, std::map<int, int>> diff_data = readDiff(diff_file);
    auto all_leaves = T.get_leaves();
    for (auto l: all_leaves) {
        //std::cout << "Data type of l: " << typeid(l).name() << std::endl;
        std::string samp = l->identifier;
        int bl = 0;
        MAT::Node* current_node = l;
        //MAT::Node* parent = l
        //std::cout << "Type of current_node: " << typeid(*current_node).name() << std::endl;
        //std::cout <<   "made it here - lk - current node init " << current_node->identifier  << std::endl;
        //std::cout << " branch len " << l -> branch_length << std::endl;
        
        //skip samples coming from root
        
        //if ((*current_node->parent).is_root()) {
        //    std::cout <<   "skipped" << current_node->identifier << std::endl;
        //    continue;
        //}
        

        //find ancestor node
        //need to make a catch so that max_snp-_dist isnt 0 when snp distance is > 0
        while (bl <= max_snp_distance) {
            //add current branch to total branch length
            int current_branch = current_node->branch_length; 
            
            //std::cout << " node" << current_node->identifier << std::endl;

            //make this catch things with root ancestor node 
            if (current_node->parent->parent == NULL) {
                std::cout <<   "ROOOOOOTTTTTT"  << std::endl;
            }

            //std::cout << " node" << current_node->parent->identifier << std::endl;
            //l = l -> parent;
            //std::cout << "Type of current_node: " << typeid(*current_node).name() << std::endl;
            
            //update current node to parent if branch length is not met 
            if ((bl + current_branch ) < snp_distance) {
                bl += current_branch;
                current_node = current_node->parent;
                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
            }

            else {
                //std::cout << " done " << l->identifier << bl << std::endl;
                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
                break;
            }
            
            //std::cout <<   "made it here - lk"  << std::endl;
            //break;

            
        }
        std::cout << " done " << l->identifier << bl << std::endl;
        std::cout << " branch len " << bl << std::endl;
        std::cout << " new node" << current_node->identifier << std::endl;
        //add a warning for this?
        if ((current_node->children).size() >20 ) {
            std::cout <<  "LOTS" << (current_node->children).size() << std::endl;
        }
        
        //need a function here to collect all missing data, 
        //dont pass diff data to dfs, pass a map of all missing data in the subtree

        //need error checking to confirm this is the ancestor node we want 
        //subtreeMask(current_node, diff_data);
        dfs(current_node, diff_data);
        //std::cout << " l: " << samp << std::endl;
        //for (const auto& word : diff_data[samp]) {
        //    std::cout << "pos " << word.first << " " << "len " << word.second << endl;
            //for (const auto& pos : word.second) {
            //    std::cout << "pos" << pos.first << '\t' << pos.second << endl;
        //    }
        //break;
        }
}
*/
//old version below 
/*
void localMask (uint32_t snp_distance, uint32_t max_snp_distance, MAT::Tree& T, std::string diff_file) {
    std::map<std::string, std::map<int, int>> diff_data = readDiff(diff_file);
    auto all_leaves = T.get_leaves();
    for (auto l: all_leaves) {
        //std::cout << "Data type of l: " << typeid(l).name() << std::endl;
        std::string samp = l->identifier;
        int bl = 0;
        MAT::Node* current_node = l;
        //MAT::Node* parent = l
        //std::cout << "Type of current_node: " << typeid(*current_node).name() << std::endl;
        //std::cout <<   "made it here - lk - current node init " << current_node->identifier  << std::endl;
        //std::cout << " branch len " << l -> branch_length << std::endl;
        
        //skip samples coming from root
        if ((*current_node->parent).is_root()) {
            std::cout <<   "skipped" << current_node->identifier << std::endl;
            continue;
        }

        //find ancestor node
        //need to make a catch so that max_snp-_dist isnt 0 when snp distance is > 0
        while (bl <= max_snp_distance) {
            //add current branch to total branch length
            int current_branch = current_node->branch_length; 
            
            //std::cout << " node" << current_node->identifier << std::endl;

            //make this catch things with root ancestor node 
            if (current_node->parent->parent == NULL) {
                std::cout <<   "ROOOOOOTTTTTT"  << std::endl;
            }

            //std::cout << " node" << current_node->parent->identifier << std::endl;
            //l = l -> parent;
            //std::cout << "Type of current_node: " << typeid(*current_node).name() << std::endl;
            
            //update current node to parent if branch length is not met 
            if ((bl + current_branch ) < snp_distance) {
                bl += current_branch;
                current_node = current_node->parent;
                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
            }

            else {
                //std::cout << " done " << l->identifier << bl << std::endl;
                //std::cout << " branch len " << bl << std::endl;
                //std::cout << " new node" << current_node->identifier << std::endl;
                break;
            }
            
            //std::cout <<   "made it here - lk"  << std::endl;
            //break;

            
        }
        std::cout << " done " << l->identifier << bl << std::endl;
        std::cout << " branch len " << bl << std::endl;
        std::cout << " new node" << current_node->identifier << std::endl;
        //add a warning for this?
        if ((current_node->children).size() >20 ) {
            std::cout <<  "LOTS" << (current_node->children).size() << std::endl;
        }
        
        //need a function here to collect all missing data, 
        //dont pass diff data to dfs, pass a map of all missing data in the subtree

        //need error checking to confirm this is the ancestor node we want 
        //subtreeMask(current_node, diff_data);
        dfs(current_node, diff_data);
        //std::cout << " l: " << samp << std::endl;
        //for (const auto& word : diff_data[samp]) {
        //    std::cout << "pos " << word.first << " " << "len " << word.second << endl;
            //for (const auto& pos : word.second) {
            //    std::cout << "pos" << pos.first << '\t' << pos.second << endl;
        //    }
        //break;
        }
}
*/

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


