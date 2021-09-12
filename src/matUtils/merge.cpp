//Load mat files
//uncondense leaves
#include "merge.hpp"
concurMap consistNodes;
tbb::reader_writer_lock rd_wr_lock;
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
    std::cout<<A_leaves.size()<<std::endl;
    std::vector<std::string> B_leaves = B.get_leaves_ids();
        std::cout<<B_leaves.size()<<std::endl;

    //creates a vector of common_leaves between two input trees
    std::vector<std::string> common_leaves;

    set_intersection(B_leaves.begin(), B_leaves.end(), A_leaves.begin(), A_leaves.end(), std::back_inserter(common_leaves));
    std::cout << common_leaves.size() << std::endl;
    if (common_leaves.size() == 0)
    {
        std::cout << common_leaves.size() << std::endl;
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
    std::cout << samples.size() << std::endl;

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
                placement_input inp;
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
                    placement(inp, false, true);
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
                i++;
                std::cout << i << ": " << finalMat.get_node(x)->parent->identifier << std::endl;
                m.unlock();
                /**
                else
                {
                    min = anc.size() - 1;
                }
                m.lock();
                finalMat.create_node(x, curr, -1);
                MAT::Node *add = finalMat.get_node(x);
                m.unlock();
                for (int i = min; i >= 0; i--)
                {
                    for (auto m : ancestors[i]->mutations)
                    {
                        add->add_mutation(m);
                    }
                }**/
            }
        },
        ap);
    std::cout << baseMat.get_num_leaves() << std::endl;
    std::cout << otherMat.get_num_leaves() << std::endl;
    std::cout << finalMat.get_num_leaves() << std::endl;

    //Save final MAT file
    finalMat.condense_leaves();

    MAT::save_mutation_annotated_tree(finalMat, output_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
void placement(placement_input &input, bool compute_parsimony_scores, bool compute_vecs)
{
    //    TIMEIT();

    // Variable to store the number of parsimony-increasing mutations to
    // place sample at the current node
    int set_difference = 0;

    // Current smallest value of the number of parsimony-increasing mutations
    // during the parallel search to place the same at some node in the tree
    int best_set_difference = *input.best_set_difference;

    std::vector<int> anc_positions;
    std::vector<MAT::Mutation> ancestral_mutations;

    // if node has some unique mutations not in new sample, placement should be
    // done as a sibling
    bool has_unique = false;
    int node_num_mut = 0;
    int num_common_mut = 0;

    // For non-root nodes, add mutations common to current node (branch) to
    // excess mutations. Set has_unique to true if a mutation unique to current
    // node not in new sample is found.

    if (!input.node->is_root())
    {
        size_t start_index = 0;
        for (auto m1 : input.node->mutations)
        {
            node_num_mut++;
            auto anc_nuc = m1.mut_nuc;
            // if mutation is masked, treat it as a unique mutation (add as
            // sibling)
            if (m1.is_masked())
            {
                has_unique = true;
                break;
            }
            assert(((anc_nuc - 1) & anc_nuc) == 0);
            bool found = false;
            bool found_pos = false;
            for (size_t k = start_index; k < input.missing_sample_mutations->size(); k++)
            {
                auto m2 = (*input.missing_sample_mutations)[k];
                start_index = k;
                if (m1.position == m2.position)
                {
                    found_pos = true;
                    if (m2.is_missing)
                    {
                        found = true;
                        num_common_mut++;
                    }
                    else
                    {
                        auto nuc = m2.mut_nuc;
                        if ((nuc & anc_nuc) != 0)
                        {
                            MAT::Mutation m;
                            m.chrom = m1.chrom;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
                            m.par_nuc = m1.par_nuc;
                            m.mut_nuc = anc_nuc;

                            ancestral_mutations.emplace_back(m);
                            anc_positions.emplace_back(m.position);
                            assert((m.mut_nuc & (m.mut_nuc - 1)) == 0);
                            if (compute_vecs)
                            {
                                (*input.excess_mutations).emplace_back(m);
                            }

                            // Ambiguous base
                            //if ((nuc & (nuc-1)) != 0) {
                            //    (*input.imputed_mutations).emplace_back(m);
                            //}
                            found = true;
                            num_common_mut++;
                            break;
                        }
                    }
                }
                if (m1.position < m2.position)
                {
                    break;
                }
            }
            if (!found)
            {
                if (!found_pos && (anc_nuc == m1.ref_nuc))
                {
                    MAT::Mutation m;
                    m.position = m1.position;
                    m.chrom = m1.chrom;
                    m.ref_nuc = m1.ref_nuc;
                    m.par_nuc = m1.par_nuc;
                    m.mut_nuc = anc_nuc;

                    ancestral_mutations.emplace_back(m);
                    anc_positions.emplace_back(m.position);
                    assert((m.mut_nuc & (m.mut_nuc - 1)) == 0);
                    if (compute_vecs)
                    {
                        (*input.excess_mutations).emplace_back(m);
                    }

                    num_common_mut++;
                }
                else
                {
                    has_unique = true;
                }
            }
        }
    }
    else
    {
        for (auto m : input.node->mutations)
        {
            ancestral_mutations.emplace_back(m);
            anc_positions.emplace_back(m.position);
        }
    }

    // Add ancestral mutations to ancestral mutations. When multiple mutations
    // at same position are found in the path leading from the root to the
    // current node, add only the most recent mutation to the vector

    {
        auto n = input.node;
        while (n->parent != NULL)
        {
            n = n->parent;
            for (auto m : n->mutations)
            {
                if (!m.is_masked() && (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end()))
                {
                    ancestral_mutations.emplace_back(m);
                    anc_positions.emplace_back(m.position);
                }
            }
        }
    }

    // sort by position. This helps speed up the search
    std::sort(ancestral_mutations.begin(), ancestral_mutations.end());

    // Iterate over missing sample mutations
    for (auto m1 : (*input.missing_sample_mutations))
    {
        // Missing bases (Ns) are ignored
        if (m1.is_missing)
        {
            continue;
        }
        size_t start_index = 0;
        bool found_pos = false;
        bool found = false;
        bool has_ref = false;
        auto anc_nuc = m1.ref_nuc;
        if ((m1.mut_nuc & m1.ref_nuc) != 0)
        {
            has_ref = true;
        }
        // Check if mutation is found in ancestral_mutations
        for (size_t k = start_index; k < ancestral_mutations.size(); k++)
        {
            auto m2 = ancestral_mutations[k];
            // Masked mutations don't match anything
            if (m2.is_masked())
            {
                continue;
            }
            start_index = k;
            if (m1.position == m2.position)
            {
                found_pos = true;
                anc_nuc = m2.mut_nuc;
                if ((m1.mut_nuc & anc_nuc) != 0)
                {
                    found = true;
                }
                break;
            }
        }
        if (found)
        {
            // If mutation is found in ancestral_mutations
            // and if the missing sample base was ambiguous,
            // add it to imputed_mutations
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0))
            {
                MAT::Mutation m;
                m.chrom = m1.chrom;
                m.position = m1.position;
                m.ref_nuc = m1.ref_nuc;
                m.par_nuc = anc_nuc;
                m.mut_nuc = anc_nuc;
                input.imputed_mutations->emplace_back(m);
            }
        }
        // If neither the same mutation nor another mutation at the same
        // position is found in ancestor but if the missing sample can carry
        // the reference allele, add a mutation with reference allele to
        // imputed_mutations for the sample (it's not a parsimony-increasing
        // mutation)
        else if (!found_pos && has_ref)
        {
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0))
            {
                MAT::Mutation m;
                m.chrom = m1.chrom;
                m.position = m1.position;
                m.ref_nuc = m1.ref_nuc;
                m.par_nuc = anc_nuc;
                m.mut_nuc = m1.ref_nuc;
                input.imputed_mutations->emplace_back(m);
            }
        }
        // In all other cases, it is a parsimony-increasing mutation. Return
        // early if number of parsimony-increasing mutations exceeds the current
        // best. Otherwise add the mutation to excess_mutations and to
        // imputed_mutations, if base was originally ambiguous
        else
        {
            MAT::Mutation m;
            m.chrom = m1.chrom;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            m.par_nuc = anc_nuc;
            if (has_ref)
            {
                m.mut_nuc = m1.ref_nuc;
            }
            else
            {
                for (int j = 0; j < 4; j++)
                {
                    if (((1 << j) & m1.mut_nuc) != 0)
                    {
                        m.mut_nuc = (1 << j);
                        ;
                        break;
                    }
                }
            }
            assert((m.mut_nuc & (m.mut_nuc - 1)) == 0);
            // If the missing sample base is ambiguous, add it to
            // imputed_mutations
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0))
            {
                input.imputed_mutations->emplace_back(m);
            }
            if (m.mut_nuc != m.par_nuc)
            {
                if (compute_vecs)
                {
                    input.excess_mutations->emplace_back(m);
                }
                set_difference += 1;
                if (!compute_parsimony_scores && (set_difference > best_set_difference))
                {
                    return;
                }
            }
        }
    }

    // For loop to add back-mutations for cases in which a mutation from the
    // root to the current node consists of a non-reference allele but no such
    // variant is found in the missing sample
    for (auto m1 : ancestral_mutations)
    {
        size_t start_index = 0;
        bool found = false;
        bool found_pos = false;
        auto anc_nuc = m1.mut_nuc;
        for (size_t k = start_index; k < input.missing_sample_mutations->size(); k++)
        {
            // If ancestral mutation is masked, terminate the search for
            // identical mutation
            if (m1.is_masked())
            {
                break;
            }
            auto m2 = (*input.missing_sample_mutations)[k];
            start_index = k;
            if (m1.position == m2.position)
            {
                found_pos = true;
                // Missing bases (Ns) are ignored
                if (m2.is_missing)
                {
                    found = true;
                    break;
                }
                if ((m2.mut_nuc & anc_nuc) != 0)
                {
                    found = true;
                }
            }
        }
        // If ancestral mutation is found, do nothing. Else, if last ancestor
        // with mutation in the same position is found having reference allele, do
        // nothing. If same position is found but not with the same allele, do
        // nothing since the mutation must be added to excess mutations already
        // in the previous loop. In all remaining cases, add the mutation to
        // excess_mutations.
        if (found)
        {
        }
        else if (!found_pos && !m1.is_masked() && (anc_nuc == m1.ref_nuc))
        {
        }
        else if (found_pos && !found)
        {
        }
        else
        {
            MAT::Mutation m;
            m.chrom = m1.chrom;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            m.par_nuc = anc_nuc;
            m.mut_nuc = m1.ref_nuc;
            assert(m.is_masked() || ((m.mut_nuc & (m.mut_nuc - 1)) == 0));
            if (m.mut_nuc != m.par_nuc)
            {
                set_difference += 1;
                if (!compute_parsimony_scores && (set_difference > best_set_difference))
                {
                    return;
                }
                if (compute_vecs)
                {
                    (*input.excess_mutations).emplace_back(m);
                }
            }
        }
    }

    // Set the number of parsimony-increasing mutations
    if (compute_parsimony_scores)
    {
        *input.set_difference = set_difference;
    }

    // if sibling of internal node or leaf, ensure it is not equivalent to placing under parent
    // if child of internal node, ensure all internal node mutations are present in the sample
    if (input.node->is_root() || ((has_unique && !input.node->is_leaf() && (num_common_mut > 0) && (node_num_mut != num_common_mut)) ||
                                  (input.node->is_leaf() && (num_common_mut > 0)) || (!has_unique && !input.node->is_leaf() && (node_num_mut == num_common_mut))))
    {
        rd_wr_lock.lock_read();
        if (set_difference > *input.best_set_difference)
        {
            rd_wr_lock.unlock();
            return;
        }
        rd_wr_lock.unlock();

        rd_wr_lock.lock();
        size_t num_leaves = input.T->get_num_leaves(input.node);
        if (set_difference < *input.best_set_difference)
        {
            *input.best_set_difference = set_difference;
            *input.best_node = input.node;
            *input.best_node_num_leaves = num_leaves;
            *input.best_j = input.j;
            *input.num_best = 1;
            *input.has_unique = has_unique;
            *input.best_distance = input.distance;
            (*input.node_has_unique)[input.j] = has_unique;
            input.best_j_vec->clear();
            input.best_j_vec->emplace_back(input.j);
        }
        else if (set_difference == *input.best_set_difference)
        {
            // Tie breaking strategy when multiple parsimony-optimal placements
            // are found. it picks the node with a greater number of descendant
            // leaves for placement. However, if the choice is between a parent
            // and its child node, it picks the parent node if the number of
            // descendant leaves of the parent that are not shared with the child
            // node exceed the number of descendant leaves of the child.
            if (((input.distance == *input.best_distance) &&
                 ((num_leaves > *input.best_node_num_leaves) ||
                  ((num_leaves == *input.best_node_num_leaves) && (*input.best_j < input.j)))) ||
                (input.distance < *input.best_distance))
            {
                *input.best_set_difference = set_difference;
                *input.best_node = input.node;
                *input.best_node_num_leaves = num_leaves;
                *input.best_j = input.j;
                *input.has_unique = has_unique;
                *input.best_distance = input.distance;
            }
            *input.num_best += 1;
            (*input.node_has_unique)[input.j] = has_unique;
            input.best_j_vec->emplace_back(input.j);
        }
        rd_wr_lock.unlock();
    }
    else if (compute_parsimony_scores)
    {
        // Add 1 to the parsimony score for this node
        // as its current best placement is equivalent
        // to placing at parent/child
        *input.set_difference = set_difference + 1;
    }
}
