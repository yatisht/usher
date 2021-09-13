//Load mat files
//uncondense leaves
#include "merge.hpp"
concurMap consistNodes;
//tbb::reader_writer_lock rd_wr_lock;
po::variables_map parse_merge_command(po::parsed_options parsed)
{
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description merge_desc("merge options");
    merge_desc.add_options()("input-mat-1", po::value<std::string>()->required(),
                             "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")("input-mat-2", po::value<std::string>()->required(),
                                                                                                                                                               "Input mutation-annotated tree file [REQUIRED]. If only this argument is set, print the count of samples and nodes in the tree.")("output-mat,o", po::value<std::string>()->required(),
                                                                                                                                                                                                                                                                                                 "Write output files to the target directory. Default is current directory.")("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str());

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try
    {
        po::store(po::command_line_parser(opts)
                      .options(merge_desc)
                      .run(),
                  vm);
        po::notify(vm);
    }
    catch (std::exception &e)
    {
        std::cerr << merge_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}
/**
 * Checks for consistency between both MAT files to ensure that they are able to merge
 **/
bool consistent(MAT::Tree A, MAT::Tree B)
{

    //vectors of all leaves in both input trees
    std::vector<std::string> A_leaves = A.get_leaves_ids();
    std::vector<std::string> B_leaves = B.get_leaves_ids();

    //creates a vector of common_leaves between two input trees
    std::vector<std::string> common_leaves;

    set_intersection(B_leaves.begin(), B_leaves.end(), A_leaves.begin(), A_leaves.end(), std::back_inserter(common_leaves));
    if (common_leaves.size() == 0)
    {
        return true;
    }
    //creates two subtrees using the common_leaves

    auto Asub = subtree(A, common_leaves);
    auto Bsub = subtree(B, common_leaves);
    Asub.rotate_for_display();
    Bsub.rotate_for_display();
    auto Adfs = Asub.depth_first_expansion();
    auto Bdfs = Bsub.depth_first_expansion();
    if (Adfs.size() != Bdfs.size())
    {
        return false;
    }
    concurMap::accessor ac;
    bool verify = true;
    static tbb::affinity_partitioner ap;
    /**
     * Parallel loop that parses through the depth first expansion of
     * both subtrees and compares each node's mutations
     **/
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, Adfs.size()),
        [&](tbb::blocked_range<size_t> r)
        {
            for (size_t k = r.begin(); k < r.end(); ++k)
            {

                verify = chelper(Adfs[k], Bdfs[k]);

                if (verify == false)
                {
                    return false;
                }
            }
        },
        ap);

    return true;
}
/**
 * Creates subtrees by removing all uncommon nodes from a 
 * copy of the original tree
 **/
MAT::Tree subtree(MAT::Tree tree, std::vector<std::string> common)
{

    tbb::mutex m;
    auto tree_copy = get_tree_copy(tree);
    auto all_leaves = tree_copy.get_leaves_ids();
    std::vector<std::string> to_remove;
    std::set_difference(all_leaves.begin(), all_leaves.end(), common.begin(), common.end(), std::back_inserter(to_remove));

    for (auto k : to_remove)
    {
        tree_copy.remove_node(k, true);
        // m.unlock();
    }
    return tree_copy;
}
/**
 * Helper function for consistency that individually compares each node's mutation
 * to make sure that the mutation paths are identical
 * Maps all consistent nodes using the consistNodes hashmap
 **/
bool chelper(MAT::Node *a, MAT::Node *b)
{
    tbb::mutex m;
    concurMap::accessor ac;
    if (a->is_root() && b->is_root())
    {
        return true;
    }
    if (a->mutations.size() != b->mutations.size())
    {
        return false;
    }
    for (int x = 0; x < a->mutations.size(); x++)
    {
        MAT::Mutation mut1 = a->mutations[x];
        MAT::Mutation mut2 = b->mutations[x];
        if (mut1.position != mut2.position)
        {
            return false;
        }
        else if (mut1.ref_nuc != mut2.ref_nuc)
        {
            return false;
        }
        else if (mut1.mut_nuc != mut2.mut_nuc)
        {
            return false;
        }
    }
    //Maps all nodes that are consistent
    consistNodes.insert(ac, b->identifier);
    ac->second = a->identifier;
    ac.release();
    return true;
}

void merge_main(po::parsed_options parsed)
{
    timer.Start();
    //Takes user input and loads the specified MAT files
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
    concurMap::accessor ac;

    if (mat1.condensed_nodes.size() > 0)
    {
        mat1.uncondense_leaves();
    }
    if (mat2.condensed_nodes.size() > 0)
    {
        mat2.uncondense_leaves();
    }
    //Assigns biggest MAT to baseMat and smaller one to otherMat
    if (mat1.get_num_leaves() > mat2.get_num_leaves())
    {
        baseMat = mat1;
        otherMat = mat2;
    }
    else
    {
        baseMat = mat2;
        otherMat = mat1;
    }
    //Checks for consistency

    if (consistent(baseMat, otherMat) == false)
    {
        fprintf(stderr, "ERROR: MAT files are not consistent!\n");
        exit(1);
    }
    
    //Makes a copy of the baseMat files to add the samples to later
    MAT::Tree finalMat = get_tree_copy(baseMat);
    std::vector<std::string> samples;
    auto otherLeaves = otherMat.get_leaves_ids();
    auto baseLeaves = baseMat.get_leaves_ids();
    
    sort(otherLeaves.begin(), otherLeaves.end());
    sort(baseLeaves.begin(), baseLeaves.end());

    //creates vector of new samples to be added
    std::set_difference(otherLeaves.begin(), otherLeaves.end(), baseLeaves.begin(), baseLeaves.end(), std::back_inserter(samples));

    //concurMap::accessor ac;
    static tbb::affinity_partitioner ap;
    fprintf(stderr, "Found %lu samples to merge\n", samples.size());
    tbb::mutex m;
    //parallel loop that concurrently parses through all the new samples
    int i = 0;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r)
        {
            for (size_t k = r.begin(); k < r.end(); k++)
            {

                std::unordered_map<int, std::string> anc;
                std::string x = samples[k];
                MAT::Node *s = otherMat.get_node(x);
                auto ancestors = otherMat.rsearch(x, true);
                std::string curr = finalMat.root->identifier;
                std::vector<int> indices;
                /** Parses through the the ancestors of each sample on otherMat and finds 
                 * Closest consistent node on baseMat 
                 **/
                for (int y = 0; y < ancestors.size(); y++)
                {
                    m.lock();
                    if (consistNodes.find(ac, ancestors[y]->identifier) == true)
                    {

                        anc[y] = ac->second;
                        indices.push_back(y);

                        ac.release();
                    }
                    m.unlock();
                }

                int min;

                if (indices.size() != 0)
                {
                    min = *std::min_element(indices.begin(), indices.end());
                    curr = anc[min];
                }

                //Roots bfs at closest consistent node identified in previous loop
                //Restricts tree search to a smaller subtree
                auto bfs = finalMat.breadth_first_expansion(curr);
                mapper2_input inp;
                //std::cout<<bfs[0]->mutations.size()<<std::endl;

                std::vector<std::vector<MAT::Mutation> > node_excess_mutations(bfs.size());
                std::vector<std::vector<MAT::Mutation> > node_imputed_mutations(bfs.size());
                size_t best_node_num_leaves = 0;
                int best_set_difference = s->mutations.size() + bfs[0]->mutations.size() + 1;
                size_t best_j = 0;
                size_t num_best = 1;
                bool best_node_has_unique = false;
                MAT::Node *best_node = bfs[0];
                std::vector<bool> node_has_unique(bfs.size(), false);
                std::vector<size_t> best_j_vec;
                best_j_vec.emplace_back(0);
                for (int k = 0; k < bfs.size(); k++)
                {
                    inp.T = &finalMat;
                    inp.node = bfs[k];
                    inp.missing_sample_mutations = &s->mutations;
                    inp.excess_mutations = &node_excess_mutations[k];
                    inp.imputed_mutations = &node_imputed_mutations[k];

                    inp.best_node_num_leaves = &best_node_num_leaves;
                    inp.best_set_difference = &best_set_difference;
                    inp.best_node = &best_node;
                    inp.best_j = &best_j;
                    inp.num_best = &num_best;
                    inp.j = k;
                    inp.has_unique = &best_node_has_unique;
                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);
                    mapper2_body(inp, false);
                }

                if (finalMat.get_node(x) == NULL)
                {
                    // Is placement as sibling
                    if (best_node->is_leaf() || best_node_has_unique)
                    {
                        std::string nid = finalMat.new_internal_node_id();
                        m.lock();
                        finalMat.create_node(nid, best_node->parent->identifier);
                        finalMat.create_node(x, nid);
                        finalMat.move_node(best_node->identifier, nid);
                        m.unlock();
                        // common_mut stores mutations common to the
                        // best node branch and the sample, l1_mut
                        // stores mutations unique to best node branch
                        // and l2_mut stores mutations unique to the
                        // sample not in best node branch
                        std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                        std::vector<MAT::Mutation> curr_l1_mut;

                        // Compute current best node branch mutations
                        for (auto m1 : best_node->mutations)
                        {
                            MAT::Mutation m = m1.copy();
                            curr_l1_mut.emplace_back(m);
                        }
                        // Clear mutations on the best node branch which
                        // will be later replaced by l1_mut
                        best_node->clear_mutations();

                        // Compute l1_mut
                        for (auto m1 : curr_l1_mut)
                        {
                            bool found = false;
                            for (auto m2 : node_excess_mutations[best_j])
                            {
                                if (m1.is_masked())
                                {
                                    break;
                                }
                                if (m1.position == m2.position)
                                {
                                    if (m1.mut_nuc == m2.mut_nuc)
                                    {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (!found)
                            {
                                MAT::Mutation m = m1.copy();
                                l1_mut.emplace_back(m);
                            }
                        }
                        // Compute l2_mut
                        for (auto m1 : node_excess_mutations[best_j])
                        {
                            bool found = false;
                            for (auto m2 : curr_l1_mut)
                            {
                                if (m1.is_masked())
                                {
                                    break;
                                }
                                if (m1.position == m2.position)
                                {
                                    if (m1.mut_nuc == m2.mut_nuc)
                                    {
                                        found = true;
                                        MAT::Mutation m = m1.copy();
                                        common_mut.emplace_back(m);
                                        break;
                                    }
                                }
                            }
                            if (!found)
                            {
                                MAT::Mutation m = m1.copy();
                                l2_mut.emplace_back(m);
                            }
                        }

                        // Add mutations to new node using common_mut
                        for (auto m : common_mut)
                        {
                            finalMat.get_node(nid)->add_mutation(m);
                        }
                        // Add mutations to best node using l1_mut
                        for (auto m : l1_mut)
                        {
                            finalMat.get_node(best_node->identifier)->add_mutation(m);
                        }
                        // Add new sample mutations using l2_mut

                        for (auto m : l2_mut)
                        {
                            finalMat.get_node(x)->add_mutation(m);
                        }
                    }
                    // Else placement as child
                    else
                    {
                        m.lock();
                        finalMat.create_node(x, best_node->identifier);
                        MAT::Node *node = finalMat.get_node(x);
                        m.unlock();
                        std::vector<MAT::Mutation> node_mut;

                        std::vector<MAT::Mutation> curr_l1_mut;

                        for (auto m1 : best_node->mutations)
                        {
                            MAT::Mutation m = m1.copy();
                            curr_l1_mut.emplace_back(m);
                        }

                        for (auto m1 : node_excess_mutations[best_j])
                        {
                            bool found = false;
                            for (auto m2 : curr_l1_mut)
                            {
                                if (m1.is_masked())
                                {
                                    break;
                                }
                                if (m1.position == m2.position)
                                {
                                    if (m1.mut_nuc == m2.mut_nuc)
                                    {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (!found)
                            {
                                MAT::Mutation m = m1.copy();
                                node_mut.emplace_back(m);
                            }
                        }

                        for (auto m : node_mut)
                        {

                            node->add_mutation(m);
                        }
                    }

                    if (node_imputed_mutations[best_j].size() > 0)
                    {
                        fprintf(stderr, "Imputed mutations:\t");
                        size_t tot = node_imputed_mutations[best_j].size();
                        for (size_t curr = 0; curr < tot; curr++)
                        {
                            if (curr < tot - 1)
                            {
                                fprintf(stderr, "%i:%c;", node_imputed_mutations[best_j][curr].position, MAT::get_nuc(node_imputed_mutations[best_j][curr].mut_nuc));
                            }
                            else
                            {
                                fprintf(stderr, "%i:%c", node_imputed_mutations[best_j][curr].position, MAT::get_nuc(node_imputed_mutations[best_j][curr].mut_nuc));
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                }
                m.lock();
                fprintf(stderr, "\rAdded %zu of %zu samples", ++i, samples.size());
                m.unlock();
            }
        },
        ap);
    //Save final MAT file
    finalMat.condense_leaves();

    MAT::save_mutation_annotated_tree(finalMat, output_filename);
    fprintf(stderr, "\nCompleted in %ld msec \n\n", timer.Stop());
}
