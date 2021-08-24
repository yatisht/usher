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

    //creates two subtrees using the common_leaves
    auto Asub = MAT::get_subtree(A, common_leaves);
    auto Bsub = MAT::get_subtree(B, common_leaves);
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
   m.lock();
    consistNodes.insert(ac, b->identifier);
    ac->second = a->identifier;
    ac.release();
    m.unlock();
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
        std::cout << "conflicts: ";

        std::cout << conflicts.size() << std::endl;
        //conflicts.clear();
        conflicts = mainhelper(conflicts, finalMat, baseMat, otherMat);
    }
    std::cout << baseMat.get_num_leaves() << std::endl;
    std::cout << otherMat.get_num_leaves() << std::endl;
    std::cout << finalMat.get_num_leaves() << std::endl;

    //Save final MAT file
    finalMat.condense_leaves();
    MAT::save_mutation_annotated_tree(finalMat, output_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

std::vector<std::string> mainhelper(std::vector<std::string> samples, MAT::Tree finalMat, MAT::Tree baseMat, MAT::Tree otherMat)
{
    //tbb::reader_writer_lock rw;
    concurMap::accessor ac;
    static tbb::affinity_partitioner ap;
    std::vector<std::string> parents;
    std::string curr;
    std::string x;
    tbb::mutex m;
    std::string child;
    std::vector<std::string> conflicts;
    auto expand = finalMat.breadth_first_expansion();

    //parallel loop that concurrently goes through all samples
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r)
        {
            for (size_t k = r.begin(); k < r.end(); k++)
            {

                curr = (expand[0])->identifier;
                m.lock();
                x = samples[k];
                auto ancestors = otherMat.rsearch(x, true);
                std::vector<MAT::Mutation> diff_mutations;
                for (auto y : ancestors)
                {
                    if (consistNodes.find(ac, y->identifier) == true)
                    {
                        curr = ac->second;
                        ac.release();
                        break;
                    }
                    for (auto z : y->mutations)
                    {
                        diff_mutations.push_back(z);
                    }
                }
                bool c = false;
                if (std::find(parents.begin(), parents.end(), curr) != parents.end())
                {
                    conflicts.push_back(x);
                    c = true;
                }
                m.unlock();
                if (c == false)
                {
                    addhelper(x, curr, finalMat, diff_mutations);
                    parents.push_back(curr);
                    diff_mutations.clear();
                }
            }
        },
        ap);
    return conflicts;
}
void addhelper(std::string x, std::string curr, MAT::Tree finalMat, std::vector<MAT::Mutation> diff)
{
    tbb::mutex m;
    finalMat.create_node(x, curr, -1);
    MAT::Node *add = finalMat.get_node(x);

    for (auto m : diff)
    {
        add->add_mutation(m);
    }
    diff.clear();
}
