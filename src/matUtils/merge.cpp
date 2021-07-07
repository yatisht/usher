  //Load mat files
 //uncondense leaves
 #include "merge.hpp"
 po::variables_map parse_merge_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description merge_desc("merge options");
    merge_desc.add_options()
        ("input-mat-1", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("input-mat-2", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("output-mat,o", po::value<std::string>()->required(),
        "Write output files to the target directory. Default is current directory.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str());
       
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(merge_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << merge_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}  

bool consistent(MAT::Tree A, MAT::Tree B){
    //vectors of all leaves in both input trees
    std::vector<std::string> A_leaves = A.get_leaves_ids();
    std::vector<std::string> B_leaves = B.get_leaves_ids();
    //creates a vector of common_leaves between two input trees
    std::vector<std::string> common_leaves(std::max(A_leaves.size(), B_leaves.size()));
    set_intersection(A_leaves.begin(),A_leaves.end(),B_leaves.begin(),B_leaves.end(), std::back_inserter(common_leaves));
    //creates two subtrees using the common_leaves
    auto Asub = MAT::get_subtree(A, common_leaves);
    auto Bsub = MAT::get_subtree(B, common_leaves);
    auto Adfs = Asub.depth_first_expansion();
    auto Bdfs = Bsub.depth_first_expansion();

    if (Adfs.size() != Bdfs.size()){
        return false;
    }
    bool verify = true;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, Adfs.size()),
            [&](tbb::blocked_range<size_t> r){
                    for (size_t k=r.begin(); k<r.end(); ++k){
                        verify = chelper(Adfs[k], Bdfs[k]);
                        if (verify == false){
                            return false;
                        }
                    }
            }, ap); 
    return true;
}

bool chelper(MAT::Node* a, MAT::Node* b){
    if (a->mutations.size() != b->mutations.size()){
        return false;
    }
    for (int x = 0; x< a->mutations.size(); x++){
            MAT::Mutation mut1 = a->mutations[x];
            MAT::Mutation mut2 = b->mutations[x];
            if (mut1.position != mut2.position){
                return false;
            }
            else if (mut1.ref_nuc != mut2.ref_nuc){
                return false;
            }
            else if (mut1.mut_nuc != mut2.mut_nuc){
                return false;
            }
            else{
                return true;
            }
    }
}

void merge_main(po::parsed_options parsed) {
    po::variables_map vm = parse_merge_command(parsed);
    std::string mat1_filename = vm["input-mat-1"].as<std::string>();
    std::string mat2_filename = vm["input-mat-2"].as<std::string>();
    std::string output_filename = vm["output-mat"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);
    MAT::Tree mat1 = MAT::load_mutation_annotated_tree(mat1_filename);
    MAT::Tree mat2 = MAT::load_mutation_annotated_tree(mat2_filename);
    MAT::Tree baseMat;
    MAT::Tree otherMat;

    if (mat1.condensed_nodes.size() > 0) {
      mat1.uncondense_leaves();
    }
    if (mat2.condensed_nodes.size() > 0) {
      mat2.uncondense_leaves();
    }
    if (consistent(mat1, mat2)==false){
        std::cout<<"false";
        return;
    }
    else{
        std::cout<<"true";
    }
    
   if (mat1.get_num_leaves() > mat2.get_num_leaves()){
        baseMat = mat1;
        otherMat = mat2;
    }
    else{
        baseMat = mat2;
        otherMat = mat1;
    }
    MAT::Tree finalMat = get_tree_copy(baseMat);
    std::vector<std::string> new_samples;
    auto otherLeaves = otherMat.get_leaves_ids();
    auto baseLeaves = baseMat.get_leaves_ids();
    sort( otherLeaves.begin(), otherLeaves.end() );
    sort( baseLeaves.begin(), baseLeaves.end() );

    //insert set
    std::set_difference(otherLeaves.begin(), otherLeaves.end(), baseLeaves.begin(), baseLeaves.end(), std::back_inserter(new_samples));
    auto expand = baseMat.breadth_first_expansion();
    //for every new sample, walk dopwn baseMat starting at the root
    for (auto x : new_samples){
        auto ancestors = otherMat.rsearch(x, true);
        auto curr = expand[0];
        std::vector<MAT::Mutation> diff_mutations;
        for (auto y : ancestors){
            //go through currents children and see if one of them matches the mutations
            for (auto z : curr->children){
               if (chelper(z,y)==true){
                   curr = z;
                   break;
               }  
               else{
                   std::vector<MAT::Mutation> tempdiff;
                   auto a = z->mutations;
                   auto b = y->mutations;
                   std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(tempdiff));
                   if (diff_mutations.size()>0){
                       if (tempdiff.size()<diff_mutations.size()){
                           diff_mutations = tempdiff;
                       }
                   }
                   if (diff_mutations.size() == 0){
                       diff_mutations = tempdiff;
                   }
               }  
            }
        }
        finalMat.create_node(x, curr->identifier, -1);
        MAT::Node* add = finalMat.get_node(x);
        add->mutations = diff_mutations;
    }
    finalMat.condense_leaves();
    MAT::save_mutation_annotated_tree(finalMat, output_filename);
    
}
