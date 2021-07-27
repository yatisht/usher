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
    std::cout << "called consistency" << std::endl;
    //vectors of all leaves in both input trees
    std::vector<std::string> A_leaves = A.get_leaves_ids();
    std::vector<std::string> B_leaves = B.get_leaves_ids();
    //creates a vector of common_leaves between two input trees
    std::vector<std::string> common_leaves(std::max(A_leaves.size(), B_leaves.size()));
    set_intersection(A_leaves.begin(), A_leaves.end(), B_leaves.begin(), B_leaves.end(), std::back_inserter(common_leaves));
    std::cout << "found common samples" << std::endl;
    //creates two subtrees using the common_leaves
    auto Asub = MAT::get_subtree(A, common_leaves);
    auto Bsub = MAT::get_subtree(B, common_leaves);
    auto Adfs = Asub.depth_first_expansion();
    auto Bdfs = Bsub.depth_first_expansion();

    if (Adfs.size() != Bdfs.size())
    {
        return false;
    }
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
    std::cout << "finished loops" << std::endl;
    return true;
}
/**
 * Checks each node to make sure their mutations are consistent with each other
 * Adds all consistent nodes to a concurrent hash map 
 **/
bool chelper(MAT::Node *a, MAT::Node *b)
{
    concurMap::accessor ac;
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
        else
        {
            consistNodes.insert(ac, b);
            ac->second = a;
            ac.release();
            return true;
        }
    }
}
void merge_main(po::parsed_options parsed)
{
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
    std::vector<std::string> new_samples;
    auto otherLeaves = otherMat.get_leaves_ids();
    auto baseLeaves = baseMat.get_leaves_ids();
    sort(otherLeaves.begin(), otherLeaves.end());
    sort(baseLeaves.begin(), baseLeaves.end());

    //creates vector of new samples to be added
    std::set_difference(otherLeaves.begin(), otherLeaves.end(), baseLeaves.begin(), baseLeaves.end(), std::back_inserter(new_samples));

    //creates a string of potential conflicting samples from the first round of merging
    std::vector<std::string> conflicts = mainhelper(new_samples, finalMat, baseMat, otherMat);

    //Continuously run mainhelper until the vector it returns is empty
    while (!(conflicts.empty()))
    {
        conflicts = mainhelper(conflicts, finalMat, baseMat, otherMat);
                std::cout<<"conflicts: ";

        std::cout<<conflicts.size()<<std::endl;
    }
    std::cout << baseMat.get_num_leaves() << std::endl;
    std::cout << otherMat.get_num_leaves() << std::endl;
    std::cout << finalMat.get_num_leaves() << std::endl;

    //Save final MAT file
    finalMat.condense_leaves();
    MAT::save_mutation_annotated_tree(finalMat, output_filename);
}

std::vector<std::string> mainhelper(std::vector<std::string> samples, MAT::Tree finalMat, MAT::Tree baseMat, MAT::Tree otherMat)
{
    concurMap::accessor ac;
    static tbb::affinity_partitioner ap;
    std::vector<MAT::Node *> parents;
    MAT::Node *curr;
    std::string x;
    int count = 0;
    std::mutex m;
    std::string child;
    std::vector<std::string> conflicts;
    auto expand = baseMat.breadth_first_expansion();
    std::cout << samples.size() << std::endl;
    //parallel loop that concurrently goes through all samples
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r)
        {
            for (size_t k = r.begin(); k < r.end(); ++k)
            {
                bool empty = true;
                std::vector<MAT::Mutation> diff_mutations;
                std::vector<MAT::Mutation> common_mutations;

                x = samples[k];
                std::cout << x << std::endl;
                //finds ancestors of new sample on otherMat
                auto ancestors = otherMat.rsearch(x, true);
                std::reverse(ancestors.begin(), ancestors.end());
                //sets curr to root of tree
                curr = expand[0];
                for (auto y : ancestors)
                {
                    //for each ancestor of the sample, checks concurrent hash map to see if there is a
                    //corresponding node in the baseMat
                    if (!(consistNodes.find(ac, y)))
                    {
                        //std::cout<<"reached ancestors"<<std::endl;
                        //after reaching a point of non-concordance, goes through children of current node
                        for (auto z : curr->children)
                        {
                            std::vector<MAT::Mutation> tempdiff;
                            std::vector<MAT::Mutation> tempcommon;
                            std::vector<MAT::Mutation> a = otherMat.get_node(x)->mutations;
                            std::vector<MAT::Mutation> b = z->mutations;
                            //std::cout<<a.size()<<std::endl;
                            // std::cout<<b.size()<<std::endl;

                            //finds child with least amount of differing mutations from sample
                            std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(tempdiff));
                            std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(tempcommon));
                            //std::cout<<tempdiff.size()<<std::endl;

                            if (tempdiff.size() < diff_mutations.size())
                            {
                                diff_mutations = tempdiff;
                                common_mutations = tempcommon;
                                child = z->identifier;
                            }

                            else if (diff_mutations.size() == 0)
                            {
                                diff_mutations = tempdiff;
                                common_mutations = tempcommon;
                                child = z->identifier;
                            }
                        }
                    }
                    else
                    {
                        curr = y;
                    }
                }
                //std::cout<<curr->identifier<<std::endl;
                //If samples intended parent is already in parents vector, adds to conflicts vector
                m.lock();
                bool c = false;
                if (common_mutations.size() > 0)
                {
                    std::string nid = std::to_string(++finalMat.curr_internal_node);
                    finalMat.create_node(nid, curr->identifier);
                    finalMat.move_node(child, nid);
                    curr = finalMat.get_node(nid);
                    for (auto m : common_mutations)
                    {
                        curr->add_mutation(m);
                    }
                    //std::cout << common_mutations.size() << std::endl;

                    common_mutations.clear();
                }
                for (auto p : parents)
                {
                    if (chelper(p, curr) == true)
                    {
                        conflicts.push_back(x);
                        c = true;
                        break;
                    }
                }
                //std::cout<<"sample is non conflicting"<<std::endl;
                //If there are any common mutations between sample and its intended neighbors,
                //creates an internal node with the common mutations and assigns it as the parent to both nodes

                if (finalMat.get_node(x) == NULL)
                {

                    if (c == false)
                    {

                        finalMat.create_node(x, curr->identifier, -1);
                        MAT::Node *add = finalMat.get_node(x);
                        //std::cout << diff_mutations.size() << std::endl;
                        for (auto m : diff_mutations)
                        {
                            add->add_mutation(m);
                        }
                        diff_mutations.clear();
                        parents.push_back(curr);
                    }
                }
                m.unlock();
            }
        },
        ap);
    //std::cout << count << std::endl;
    return conflicts;
}
