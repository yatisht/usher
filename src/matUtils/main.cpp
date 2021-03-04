#include "annotate.hpp"
#include "mask.hpp"
#include "convert.hpp"
#include "filter.hpp"
#include "describe.hpp"
#include "uncertainty.hpp"
#include "select.hpp"
#include "summary.hpp"

Timer timer; 

po::variables_map parse_extract_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    //it will attempt to do something dumb with casting if the default prefix is ./ for the output directory
    //so I need to predefine that default as an actual string.

    po::variables_map vm;
    po::options_description conv_desc("extract options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("samples,s", po::value<std::string>()->default_value(""),
        "Select samples by explicitly naming them. One per line")
        ("clade,c", po::value<std::string>()->default_value(""),
        "Select samples by membership in at least one of the indicated clade(s), comma delimited.")
        ("mutation,m", po::value<std::string>()->default_value(""),
        "Select samples by whether they contain any of the indicated mutation(s), comma delimited.")
        ("max-epps,e", po::value<size_t>()->default_value(0),
        "Select samples by whether they have less than the maximum indicated number of equally parsimonious placements. Note: calculation adds significantly to runtime.")
        ("max-parsimony,a", po::value<int>()->default_value(-1),
        "Select samples by whether they have less than the maximum indicated parsimony score (terminal branch length)")
        ("get-representative,r", po::bool_switch(),
        "Automatically select two representative samples per clade in the tree after other selection steps and prune all other samples.")
        ("prune,p", po::bool_switch(),
        "Remove the selected samples instead of keeping them in the output files.")
        ("output-directory,d", po::value<std::string>()->default_value("./"),
        "Write output files to the target directory. Default is current directory. NOTE: The target directory must already exist.")
        ("sample-paths,S", po::value<std::string>()->default_value(""),
        "Write the path of mutations defining each selected sample to the target file (all samples if no selection arguments)")
        ("clade-paths,C", po::value<std::string>()->default_value(""),
        "Write the path of mutations defining each clade in the tree after sample selection to the target file.")
        ("all-paths,A", po::value<std::string>()->default_value(""),
        "Write mutations assigned to each node in the selected tree to the target file.")
        ("write-vcf,v", po::value<std::string>()->default_value(""),
         "Output VCF file representing selected subtree. Default is full tree")
        ("no-genotypes,n", po::bool_switch(),
        "Do not include sample genotype columns in VCF output. Used only with the write-vcf option")
        ("collapse-tree,O", po::bool_switch(),
        "Collapse the MAT before writing it to output. Used only with the write-mat option")
        ("write-mat,o", po::value<std::string>()->default_value(""),
        "Write the selected tree as a new protobuf to the target file.")
        ("write-tree,t", po::value<std::string>()->default_value(""),
         "Use to write a newick tree to the indicated file.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
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

void extract_main (po::parsed_options parsed) {
    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parse_extract_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string input_samples_file = vm["samples"].as<std::string>();
    std::string clade_choice = vm["clade"].as<std::string>();
    std::string mutation_choice = vm["mutation"].as<std::string>();
    int max_parsimony = vm["max-parsimony"].as<int>();
    size_t max_epps = vm["max-epps"].as<size_t>();
    bool prune_samples = vm["prune"].as<bool>();
    bool get_representative = vm["get-representative"].as<bool>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    //if the input on dir_prefix doesn't end with a /, go ahead and add that manually
    if (dir_prefix[dir_prefix.size()-1] != '/') {
        dir_prefix.append("/");
    }
    std::string sample_path_filename = dir_prefix + vm["sample-paths"].as<std::string>();
    std::string clade_path_filename = dir_prefix + vm["clade-paths"].as<std::string>();
    std::string all_path_filename = dir_prefix + vm["all-paths"].as<std::string>();
    std::string tree_filename = dir_prefix + vm["write-tree"].as<std::string>();
    std::string vcf_filename = dir_prefix + vm["write-vcf"].as<std::string>();
    std::string output_mat_filename = dir_prefix + vm["write-mat"].as<std::string>();
    bool collapse_tree = vm["collapse-tree"].as<bool>();
    bool no_genotypes = vm["no-genotypes"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    //check that at least one of the output filenames (things which take dir_prefix)
    //are set before proceeding. 
    std::vector<std::string> outs = {sample_path_filename, clade_path_filename, all_path_filename, tree_filename, vcf_filename, output_mat_filename};
    if (!std::any_of(outs.begin(), outs.end(), [=](std::string f){return f != dir_prefix;})) {
        fprintf(stderr, "ERROR: No output files requested!\n");
        exit(1);
    }

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str()); 
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    //first step is to select a set of samples to work with
    //if any sample selection arguments are set.
    //each of the sample selection arguments is run sequentially
    //so that you can choose, for example, "from among this clade, take samples with EPPs = 1..."
    //explicit setting goes first, then clade, then mutation, then epps- in order of increasing runtime for checking sample membership
    //to maximize efficiency (no need to calculate epps for a sample you're going to throw out in the next statement)
    timer.Start();
    fprintf(stderr, "Checking for and applying sample selection arguments\n");
    //TODO: sample select code could take a previous list of samples to check in some cases like what EPPs does
    //which could save on runtime compared to getting instances across the whole tree and intersecting
    //though the current setup would enable more complex things, like saying this OR this in arguments, it's less efficient
    //and arguably efficiency should be prioritized over flexibility for features that are not implemented
    std::vector<std::string> samples;
    if (input_samples_file != "") {
        samples = read_sample_names(input_samples_file);
    } 
    if (clade_choice != "") {
        //process the clade choice string
        std::vector<std::string> clades;
        std::stringstream cns(clade_choice);
        std::string n;
        while (std::getline(cns,n,',')) {
            clades.push_back(n);
        }
        assert (clades.size() > 0);
        //this has a bit different setup from the rest of selection because these are OR logic
        //e.g. if I pass 20B,19D I want samples that belong to one OR the other
        //so I create a vector which represents all samples in both, then intersect that vector with the samples selected otherwise
        std::vector<std::string> samples_in_clade;
        for (auto cname: clades) {
            fprintf(stderr, "Getting member samples of clade %s\n", cname.c_str());
            auto csamples = get_clade_samples(T, cname);
            if (csamples.size() == 0) {
                //warning because they may want the other clade they indicated and it can proceed with that
                //itll error down the line if this was the only one they passed in and it leaves them with no samples
                fprintf(stderr, "WARNING: Clade %s is not detected in the input tree!\n", cname.c_str());
            }
            samples_in_clade.insert(samples_in_clade.end(), csamples.begin(), csamples.end());
        }
        //remove duplicate samples
        std::set<std::string> tset(samples_in_clade.begin(), samples_in_clade.end());
        samples_in_clade.assign(tset.begin(),tset.end());

        //proceed to the normal intersection code
        if (samples.size() == 0) {
            samples = samples_in_clade;
        } else {
            samples = sample_intersect(samples, samples_in_clade);
        }
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (mutation_choice != "") {
        //same situation as clades, above (comma-sep mutations are OR)
        std::vector<std::string> mutations;
        std::stringstream mns(mutation_choice);
        std::string m;
        while (std::getline(mns,m,',')) {
            mutations.push_back(m);
        }
        assert (mutations.size() > 0);
        //std::cerr << mutations[0];

        std::vector<std::string> samples_with_mutation;
        for (auto mname: mutations) {
            fprintf(stderr, "Getting samples with mutation %s\n", mname.c_str());
            auto msamples = get_mutation_samples(T, mname);
            if (msamples.size() == 0) {
                fprintf(stderr, "WARNING: No samples with mutation %s found in tree!\n", mname.c_str());
            }
            samples_with_mutation.insert(samples_with_mutation.end(), msamples.begin(), msamples.end());
        }
        //remove duplicate samples
        std::set<std::string> tset(samples_with_mutation.begin(), samples_with_mutation.end());
        samples_with_mutation.assign(tset.begin(),tset.end());

        if (samples.size() == 0) {
            samples = samples_with_mutation;
        } else {
            samples = sample_intersect(samples, samples_with_mutation);
        }
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (max_parsimony >= 0) {
        if (samples.size() == 0) {
            samples = get_parsimony_samples(T, max_parsimony);
        } else {
            auto psamples =  get_parsimony_samples(T, max_parsimony);
            samples = sample_intersect(samples, psamples);
            //check to make sure we haven't emptied our sample set; if we have, throw an error
            if (samples.size() == 0) {
                fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
                exit(1);
            }
        }
    } 
    if (max_epps > 0) {
        //this specific sample parser only calculates for values present in samples argument
        //so it doesn't need any intersection code
        //if no samples are indicated, you get it for the whole tree. This can take several hours.
        //check to make sure we haven't emptied our sample set; if we have, throw an error
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    //the final step of selection is to invert the set if prune is set
    //this is done by getting all sample names which are not in the samples vector.
    if (prune_samples) {
        fprintf(stderr, "Sample pruning requested...\n");
        std::vector<std::string> nsamples;
        for (auto s: T.get_leaves_ids()) {
            //for every sample in the tree, if that sample is NOT in the selected set, save it
            if (std::find(samples.begin(), samples.end(), s) == samples.end()) {
                nsamples.push_back(s);
            }
        }
        //and overwrite our choice of samples
        samples = nsamples;
        //check to make sure we haven't emptied our sample set; if we have, throw an error
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria (tree is left empty). Change arguments and try again\n");
            exit(1);
        }
    }
    //if we have no samples at the end without throwing an error, 
    //probably because no selection arguments were set,
    //we want the samples to be all samples in the tree
    //and don't bother pruning.

    MAT::Tree subtree; 
    if (samples.size() == 0) {
        samples = T.get_leaves_ids();
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        //if no selection was set, then there's no need to filter.
        //subtree is identical to tree.
        subtree = T;
    } else {
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        //timer.Start();
        //fprintf(stderr, "Extracting subtree...\n");
        //second step is to filter the input based on the samples

        //TODO: filter_master can support pruning directly rather than generating an inverse set and pruning all but,
        //but downstream stuff wants the inverse set of sample names to work with
        //so there's a better way to do this in at least some cases.
        subtree = filter_master(T, samples, false);
        //fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    //getting clade representative samples is a special additional step after constructing a subtree
    //so that it selects the first two nodes which fulfilled all other conditions and were retained in the initial filtering
    //TODO: there should be a better way to integrate this information that doesn't involve multiple steps of filtering
    //its basically just that the selection of two samples is very dependent on the other samples which were valid
    //add a valid samples vector argument to the clade representatives function and check membership before selection, maybe
    if (get_representative) {
        //fprintf(stderr, "Filtering again to a clade representative tree...\n");
        auto rep_samples = get_clade_representatives(subtree);
        //run filter master again
        subtree = filter_master(subtree, rep_samples, false);
    }
    //if additional information was requested, save it to the target files
    if (sample_path_filename != dir_prefix || clade_path_filename != dir_prefix || all_path_filename != dir_prefix) {
        timer.Start();
        fprintf(stderr,"Retriving path information...\n");
        if (sample_path_filename != dir_prefix) {
            std::ofstream outfile (sample_path_filename);
            auto mpaths = mutation_paths(subtree, samples);
            for (auto mstr: mpaths) {
                outfile << mstr << "\n";
            }
            outfile.close();
        }
        if (all_path_filename != dir_prefix) {
            std::ofstream outfile (all_path_filename);
            auto apaths = all_nodes_paths(subtree);
            for (auto astr: apaths) {
                outfile << astr << "\n";
            }
            outfile.close();
        }
        if (clade_path_filename != dir_prefix) {
            //need to get the set of all clade annotations currently in the tree for the clade_paths function
            //TODO: refactor summary so I can just import the clade counter function from there (disentangle from file printing)
            //also this block of code just generally sucks.
            std::vector<std::string> cladenames;
            auto dfs = subtree.depth_first_expansion();
            for (auto s: dfs) {
                std::vector<std::string> canns = s->clade_annotations;
                if (canns.size() > 0) {
                    if (canns.size() > 1 || canns[0] != "") {
                        for (auto c: canns) {
                            if (c != "" && std::find(cladenames.begin(), cladenames.end(), c) == cladenames.end()) {
                                cladenames.push_back(c);
                            }
                        }
                    }
                }
            }
            std::ofstream outfile (clade_path_filename);
            //TODO: maybe better to update clade_paths to take an "all current clades" option.
            //should be more computationally efficient at least.
            auto cpaths = clade_paths(subtree, cladenames);
            for (auto cstr: cpaths) {
                outfile << cstr;
            }
            outfile.close(); 
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //last step is to convert the subtree to other file formats
    if (vcf_filename != dir_prefix) {
        fprintf(stderr, "Generating VCF of final tree\n");
        make_vcf(subtree, vcf_filename, no_genotypes);
    }
    if (tree_filename != dir_prefix) {
        fprintf(stderr, "Generating Newick file of final tree\n");
        FILE *tree_file = fopen(tree_filename.c_str(), "w");
        fprintf(tree_file, "%s\n",
            MAT::get_newick_string(subtree, true, true, true).c_str());
        fclose(tree_file);        
    }
    //and save a MAT if that was set
    if (output_mat_filename != dir_prefix) {
        fprintf(stderr, "Saving output MAT file %s.\n", output_mat_filename.c_str()); 
        subtree.condense_leaves();
        if (collapse_tree) {
            subtree.collapse_tree();
        }
        MAT::save_mutation_annotated_tree(subtree, output_mat_filename);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}


int main (int argc, char** argv) {
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, mask, convert, prune and describe.")
        ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    std::string cmd;
    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();
    //this help string shows up over and over, lets just define it once
    std::string helpstr = "matUtils has several valid subcommands: \n\n"
        "extract: subsets the input MAT on various conditions and/or converts to other formats (MAT, newick, VCF, etc)\n\n"
        "summary: calculates basic statistics and counts members in the input MAT\n\n"
        "annotate: assigns clade identities to nodes, directly or by inference\n\n"
        "uncertainty: calculates sample placement uncertainty metrics and writes the results to tsv\n\n"
        "mask: masks the input samples\n\n"
        "Individual command options can be accessed with matUtils command --help, e.g. matUtils annotate --help will show annotation-specific help messages.\n\n";
    
    try {
        po::store(parsed, vm);
        cmd = vm["command"].as<std::string>();
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Help follows:\n\n");
        fprintf(stderr, "%s", helpstr.c_str());
        //0 when no command is selected because that's what passes tests.
        exit(0);
    }
    if (cmd == "extract") {
        extract_main(parsed);
    } else if (cmd == "annotate") {
        annotate_main(parsed);
    } else if (cmd == "mask"){
        mask_main(parsed); 
    } else if (cmd == "uncertainty") {
        uncertainty_main(parsed);
    } else if (cmd == "summary") {
        summary_main(parsed);
    } else if (cmd == "help") { 
        fprintf(stderr, "%s", helpstr.c_str());
        exit(0);
    } else {
        fprintf(stderr, "Invalid command. Help follows:\n\n");
        fprintf(stderr, "%s", helpstr.c_str());
        exit(1);
    }

    return 0;
}

