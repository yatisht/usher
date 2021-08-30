//Load mat files
//uncondense leaves
#include "merge.hpp"
concurMap consistNodes;

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
    std::vector<std::string> newleaves;

    set_intersection(A_leaves.begin(), A_leaves.end(), B_leaves.begin(), B_leaves.end(), std::back_inserter(common_leaves));
    std::cout << common_leaves.size() << std::endl;
    //creates two subtrees using the common_leaves
    if (common_leaves.size() == 0)
    {
        return true;
    }
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
    // m.lock();
    consistNodes.insert(ac, b->identifier);
    ac->second = a->identifier;
    ac.release();
    //  m.unlock();
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
        return;
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
     std::cout<<samples.size()<<std::endl;

    tbb::mutex m;
    //parallel loop that concurrently goes through all samples
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r)
        {
            for (size_t k = r.begin(); k < r.end(); k++)
            {

                std::unordered_map<int, std::string> anc;

                std::string x = samples[k];
                auto ancestors = otherMat.rsearch(x, true);
                std::string curr = finalMat.root->identifier;
                std::vector<int> indices;
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
                else
                {
                    min = anc.size() - 1;
                }
                m.lock();
                finalMat.create_node(x, curr, -1);
                MAT::Node *add = finalMat.get_node(x);
                m.unlock();
                for (int i = min; i > 0; i--)
                {
                    for (auto m : ancestors[i]->mutations)
                    {
                        add->add_mutation(m);
                    }
                }
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
