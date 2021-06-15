 //Load mat files
 //uncondense leaves
 #include "merge.hpp"
 po::variables_map parse_merge_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description merge_desc("merge options");
    merge_desc.add_options()
        ("input-mat-1,i1", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("input-mat-2,i2", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")
        ("output-mat,o", po::value<std::string>()->required(),
        "Write output files to the target directory. Default is current directory.")
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

bool consistent(MAT::Tree A, MAT::Tree B){
    //vectors of all leaves in both input trees
    auto A_leaves = A->get_leaves();
    auto B_leaves = B->get_leaves();
    //creates a vector of common_leaves between two input trees
    std::vector<auto> common_leaves;
    set_intersection(A_leaves.begin(),B_leaves.end(),B_leaves.begin(),B_leaves.end(), std::back_inserter(common_leaves));
    //creates two subtrees using the common_leaves
    auto Asub = get_subtree(A, common_leaves);
    auto Bsub = get_subtree(B, common_leaves);
    
    auto Adfs = Asub.depth_first_expansion();
    auto Bdfs = Bsub.depth_first_expansion();
    //rotate tree
    if (Adfs.size() != Bdfs.size()){
        return false;
    }
    for (int i = 0; i < Adfs.size(); i++){
        if ((Adfs[i].is_leave())){
            if (Bdfs[i].is_leave()){
                if (Adfs[i]->mutations == Bdfs[i]->mutations){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

}
void merge_main(po::parsed_options parsed) {
    po::variables_map vm = parse_mask_command(parsed);
    std::string mat1_filename = vm["input-mat-1"].as<std::string>();
    std::string mat2_filename = vm["input-mat-2"].as<std::string>();
    std::string output_filename = vm["output-mat"].as<std::string>();

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
        return;
    }
    if (mat1.get_num_leaves() > mat2.get_num_leaves()){
        baseMat = mat1;
        otherMat = mat2;
    }
    else{
        baseMat = mat2;
        otherMat = mat1;
    }

    auto new_samples = setdiff(otherMat->get_leaves(), baseMat->get_leaves());
    auto expand = baseMat.breadth_first_expansion();
    for (auto x : new_samples){
        auto ancestors = rsearch(x->identifier, true);
        auto curr = expand[0];
        auto diff_mutations;
        for (auto y : ancestors){
            for (z : curr->children){
               if (z->mutations == y->mutations){
                   curr = z;
               }  
               else{
                   diff_mutations = setdiff(z->mutations, y->mutations);
                   break;
               }  
            }
        }
        baseMat->create_node(x->identifier, curr->identifier, 0);
        x->mutations = diff_mutations;
    }

}
