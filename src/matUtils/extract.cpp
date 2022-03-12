#include "extract.hpp"

po::variables_map parse_extract_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("extract options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("input-gtf,g", po::value<std::string>()->default_value(""),
     "Input GTF annotations file (only used with --taxodium)")
    ("input-fasta,f", po::value<std::string>()->default_value(""),
     "Input FASTA reference file (only used with --taxodium)")
    ("samples,s", po::value<std::string>()->default_value(""),
     "Select samples by explicitly naming them. One per line")
    ("metadata,M", po::value<std::string>()->default_value(""),
     "Comma-delineated paths to metadata tsv/csvs containing categorical metadata values for a json or taxodium output. Used with -j and -l only")
    ("clade,c", po::value<std::string>()->default_value(""),
     "Select samples by membership in at least one of the indicated clade(s), comma delimited.")
    ("mutation,m", po::value<std::string>()->default_value(""),
     "Select samples by whether they contain any of the indicated mutation(s), comma delimited.")
    ("match,H", po::value<std::string>()->default_value(""),
     "Select samples by whether their identifier matches the indicated regex pattern.")
    ("max-epps,e", po::value<size_t>()->default_value(0),
     "Select samples by whether they have less than the maximum indicated number of equally parsimonious placements. Note: calculation adds significantly to runtime.")
    ("max-parsimony,a", po::value<int>()->default_value(-1),
     "Select samples by whether they have less than the maximum indicated parsimony score (terminal branch length)")
    ("max-branch-length,b", po::value<int>()->default_value(-1),
     "Select samples which don't have any branches with than the indicated length in their ancestry.")
    ("max-path-length,P", po::value<int>()->default_value(-1),
     "Select samples which have a total path length (number of mutations different from reference) less than or equal to P.")
    ("nearest-k,k", po::value<std::string>()->default_value(""),
     "Select a sample ID and the nearest k samples to it, formatted as sample:k. E.g. -k sample_1:50 gets sample 1 and the nearest 50 samples to it as a subtree.")
    ("nearest-k-batch,K", po::value<std::string>()->default_value(""),
     "Pass a text file of sample IDs and a number of the number of context samples, formatted as sample_file.txt:k.")
    ("set-size,z", po::value<size_t>()->default_value(0),
     "Automatically add or remove samples at random from the selected sample set until it is the indicated size.")
    ("limit-to-lca,Z", po::bool_switch(),
     "Use to limit random samples chosen with -z or -w to below the most recent common ancestor of all other samples.")
    ("get-internal-descendents,I", po::value<std::string>()->default_value(""),
     "Select the set of samples descended from the indicated internal node.")
    ("from-mrca,U", po::bool_switch(),
     "Select all samples which are descended from the most recent common ancestor of the indicated set of samples. Applied before filling background with random samples.")
    ("get-representative,r", po::value<size_t>()->default_value(0),
     "Automatically select the indicated number of representative samples per clade in the tree after other selection steps and prune all other samples (minimum 2).")
    ("prune,p", po::bool_switch(),
     "Remove the selected samples instead of keeping them in the output files.")
    ("resolve-polytomies,R", po::bool_switch(),
     "Resolve all polytomies by assigning branch length 0 relationships arbitrarily. Applied after selection; prevents recondensing of the MAT.")
    ("output-directory,d", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("used-samples,u", po::value<std::string>()->default_value(""),
     "Write a simple text file of selected sample ids.")
    ("sample-paths,S", po::value<std::string>()->default_value(""),
     "Write the path of mutations defining each sample in the subtree.")
    ("clade-paths,C", po::value<std::string>()->default_value(""),
     "Write the path of mutations defining each clade in the subtree to the target file.")
    ("all-paths,A", po::value<std::string>()->default_value(""),
     "Write mutations assigned to each node in the subtree in depth-first traversal order to the target file.")
    ("write-vcf,v", po::value<std::string>()->default_value(""),
     "Output VCF file representing selected subtree. Default is full tree")
    ("no-genotypes,n", po::bool_switch(),
     "Do not include sample genotype columns in VCF output. Used only with the write-vcf option")
    ("collapse-tree,O", po::bool_switch(),
     "Collapse the MAT before writing it to output. Used only with the write-mat option")
    ("write-mat,o", po::value<std::string>()->default_value(""),
     "Write the selected tree as a new protobuf to the target file.")
    ("write-json,j", po::value<std::string>()->default_value(""),
     "Write the tree as a JSON to the indicated file.")
    ("write-tree,t", po::value<std::string>()->default_value(""),
     "Use to write a newick tree to the indicated file.")
    ("write-taxodium,l", po::value<std::string>()->default_value(""),
     "Write protobuf in alternate format consumed by Taxodium.")
    ("x-scale,G", po::value<float>()->default_value(0.2),
     "Specifies custom X-axis scaling value for Taxodium output. Not necessary for UShER SARS-CoV-2 trees.")
    ("title,B", po::value<std::string>()->default_value("mutation_annotated_tree"),
     "Title of MAT to display in Taxodium or Auspice (used with --write-taxodium or -j).")
    ("description,D", po::value<std::string>()->default_value(""),
     "Description of MAT to display in Taxodium (used with --write-taxodium).")
    ("include-nt,J", po::bool_switch(),
     "Include nucleotide mutations in addition to protein mutations (used with --write-taxodium)")
    ("extra-fields,F", po::value<std::string>()->default_value(""),
     "Comma-separated list of additional metadata column names beyond the defaults of: strain, genbank_accession, country, date, pangolin_lineage (used with --write-taxodium).")
    ("retain-branch-length,E", po::bool_switch(),
     "Use to not recalculate branch lengths when saving newick output. Used only with -t")
    ("minimum-subtrees-size,N", po::value<size_t>()->default_value(0),
     "Use to generate a series of JSON or Newick format files representing subtrees of the indicated size which cover all queried samples. Uses and overrides -j and -t output arguments.")
    ("reroot,y", po::value<std::string>()->default_value(""),
     "Indicate an internal node ID to reroot the output tree to. Applied before all other manipulation steps.")
    ("usher-single-subtree-size,X", po::value<size_t>()->default_value(0),
     "Use to produce an usher-style single sample subtree of the indicated size with all selected samples plus random samples to fill. Produces .nh and .txt files.")
    ("usher-minimum-subtrees-size,x", po::value<size_t>()->default_value(0),
     "Use to produce an usher-style minimum set of subtrees of the indicated size which include all of the selected samples. Produces .nh and .txt files.")
    ("usher-clades-txt", po::bool_switch(),
     "When producing usher-style subtree(s), also write an usher-style clades.txt file with clade annotations for selected samples, if the tree has clade annotations.")
    ("add-random,W", po::value<size_t>()->default_value(0),
     "Add exactly W samples at random to your selection. Affected by -Z and overridden by -z.")
    ("select-nearest,Y", po::value<size_t>()->default_value(0),
     "Set to add to the sample selection the y nearest samples to each of your samples, without duplicates.")
    ("closest-relatives,V", po::value<std::string>()->default_value(""),
     "Write a tsv file of the closest relative(s) (in mutations) of each selected sample to the indicated file. All equidistant closest samples are included unless --break-ties is set.")
    ("break-ties,q", po::bool_switch(),
     "Only output one closest relative per sample (used only with --closest-relatives). If multiple closest relatives are equidistant, the lexicographically smallest sample ID is chosen.")
    ("dump-metadata,Q", po::value<std::string>()->default_value(""),
     "Set to write all final stored metadata to a tsv.")
    ("whitelist,L", po::value<std::string>()->default_value(""),
     "Pass a list of samples, one per line, to always retain regardless of any other parameters.")
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
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
    std::string whitelist_samples_file = vm["whitelist"].as<std::string>();
    std::string nearest_k = vm["nearest-k"].as<std::string>();
    std::string nearest_k_batch_file = vm["nearest-k-batch"].as<std::string>();
    std::string clade_choice = vm["clade"].as<std::string>();
    std::string mutation_choice = vm["mutation"].as<std::string>();
    std::string match_choice = vm["match"].as<std::string>();
    std::string internal_choice = vm["get-internal-descendents"].as<std::string>();
    std::string reroot_node = vm["reroot"].as<std::string>();
    int max_parsimony = vm["max-parsimony"].as<int>();
    int max_branch = vm["max-branch-length"].as<int>();
    int max_path = vm["max-path-length"].as<int>();
    size_t max_epps = vm["max-epps"].as<size_t>();
    bool prune_samples = vm["prune"].as<bool>();
    bool mrca = vm["from-mrca"].as<bool>();
    size_t get_representative = vm["get-representative"].as<size_t>();
    bool resolve_polytomies = vm["resolve-polytomies"].as<bool>();
    bool retain_branch = vm["retain-branch-length"].as<bool>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    size_t usher_single_subtree_size = vm["usher-single-subtree-size"].as<size_t>();
    size_t usher_minimum_subtrees_size = vm["usher-minimum-subtrees-size"].as<size_t>();
    bool usher_clades_txt = vm["usher-clades-txt"].as<bool>();
    size_t setsize = vm["set-size"].as<size_t>();
    size_t minimum_subtrees_size = vm["minimum-subtrees-size"].as<size_t>();
    bool limit_lca = vm["limit-to-lca"].as<bool>();
    size_t add_random = vm["add-random"].as<size_t>();
    size_t select_nearest = vm["select-nearest"].as<size_t>();
    float x_scale = vm["x-scale"].as<float>();
    bool break_ties = vm["break-ties"].as<bool>();
    bool include_nt = vm["include-nt"].as<bool>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string used_sample_filename = dir_prefix + vm["used-samples"].as<std::string>();
    std::string sample_path_filename = dir_prefix + vm["sample-paths"].as<std::string>();
    std::string clade_path_filename = dir_prefix + vm["clade-paths"].as<std::string>();
    std::string all_path_filename = dir_prefix + vm["all-paths"].as<std::string>();
    std::string closest_relatives_filename = dir_prefix + vm["closest-relatives"].as<std::string>();
    std::string tree_filename = dir_prefix + vm["write-tree"].as<std::string>();
    std::string vcf_filename = dir_prefix + vm["write-vcf"].as<std::string>();
    std::string output_mat_filename = dir_prefix + vm["write-mat"].as<std::string>();
    std::string output_tax_filename = dir_prefix + vm["write-taxodium"].as<std::string>();
    std::string tax_title = vm["title"].as<std::string>();
    std::string tax_description = vm["description"].as<std::string>();
    std::string raw_meta_fields = vm["extra-fields"].as<std::string>();
    std::string json_filename = dir_prefix + vm["write-json"].as<std::string>();
    std::string meta_filename = vm["metadata"].as<std::string>();
    std::string gtf_filename = dir_prefix + vm["input-gtf"].as<std::string>();
    std::string fasta_filename = dir_prefix + vm["input-fasta"].as<std::string>();
    std::string dump_metadata = dir_prefix + vm["dump-metadata"].as<std::string>();


    std::vector<std::string> additional_meta_fields;
    MAT::string_split(raw_meta_fields, ',', additional_meta_fields);

    bool collapse_tree = vm["collapse-tree"].as<bool>();
    bool no_genotypes = vm["no-genotypes"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    //check that at least one of the output filenames (things which take dir_prefix)
    //are set before proceeding.
    std::vector<std::string> outs = {sample_path_filename, clade_path_filename, all_path_filename, tree_filename, vcf_filename, output_mat_filename, output_tax_filename, json_filename, used_sample_filename};
    if (!std::any_of(outs.begin(), outs.end(), [=](std::string f) {
    return f != dir_prefix;
}) &&
usher_single_subtree_size == 0 && usher_minimum_subtrees_size == 0) {
        if (nearest_k_batch_file == "" && closest_relatives_filename == dir_prefix) {
            fprintf(stderr, "ERROR: No output files requested!\n");
            exit(1);
        }
    }

    tbb::task_scheduler_init init(num_threads);

    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    // Load input MAT and uncondense tree
    MAT::Tree T;
    if (input_mat_filename.find(".pb\0") != std::string::npos) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
    } else if (input_mat_filename.find(".json\0") != std::string::npos) {
        T = load_mat_from_json(input_mat_filename);
    } else {
        fprintf(stderr, "ERROR: Input file ending not recognized. Must be .json or .pb\n");
        exit(1);
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    //first step is to select a set of samples to work with
    //if any sample selection arguments are set.
    //each of the sample selection arguments is run sequentially
    //so that you can choose, for example, "from among this clade, take samples with EPPs = 1..."
    //explicit setting goes first, then they are calculated in rough order of time to compute, with EPPs at the end, so that
    //expensive computation isn't performed on samples that are going to be removed anyways.
    timer.Start();
    fprintf(stderr, "Checking for and applying sample selection arguments\n");
    std::vector<std::string> samples;
    if (input_samples_file != "") {
        samples = read_sample_names(input_samples_file);
    }
    std::string sample_id;
    if (nearest_k != "") {
        auto split_point = nearest_k.find(":");
        if (split_point == std::string::npos) {
            fprintf(stderr, "ERROR: Invalid formatting of -k argument. Requires input in the form of 'sample_id:k' to get k nearest samples to sample_id\n");
            exit(1);
        }
        sample_id = nearest_k.substr(0, split_point);
        std::string nkstr = nearest_k.substr(split_point+1, nearest_k.size() - split_point);
        int nk = std::stoi(nkstr);
        if (nk <= 0) {
            fprintf(stderr, "ERROR: Invalid neighborhood size. Please choose a positive nonzero integer.\n");
            exit(1);
        }
        auto nk_samples = get_nearby(&T, sample_id, nk);
        assert ( nk_samples.size() > 0 ) ;
        if (samples.size() == 0) {
            samples = nk_samples;
        } else {
            samples = sample_intersect(samples, nk_samples);
        }
    }
    if (internal_choice != "") {
        auto ic_samples = T.get_leaves_ids(internal_choice);
        if (samples.size() == 0) {
            samples = ic_samples;
        } else {
            samples = sample_intersect(samples, ic_samples);
        }
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
            auto csamples = get_clade_samples(&T, cname);
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

        std::vector<std::string> samples_with_mutation;
        for (auto mname: mutations) {
            fprintf(stderr, "Getting samples with mutation %s\n", mname.c_str());
            auto msamples = get_mutation_samples(&T, mname);
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
    if (match_choice != "") {
        samples = get_sample_match(&T, samples, match_choice);
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (max_parsimony >= 0) {
        samples = get_parsimony_samples(&T, samples, max_parsimony);
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (max_branch >= 0) {
        samples = get_short_steppers(&T, samples, max_branch);
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (max_path >= 0) {
        samples = get_short_paths(&T, samples, max_path);
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (max_epps > 0) {
        //this specific sample parser only calculates for values present in samples argument
        //so it doesn't need any intersection code
        //if no samples are indicated, you get it for the whole tree. This can take several hours.
        //check to make sure we haven't emptied our sample set; if we have, throw an error
        samples = get_samples_epps(&T, max_epps, samples);
        if (samples.size() == 0) {
            fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
            exit(1);
        }
    }
    if (mrca) {
        //get all descendents from the most recent common ancestor of the current samples.
        //there's no reason to check if the set is empty because this set is always a minimum what's passed in (its additive only)
        //perform this step before filling in background with random samples
        samples = get_mrca_samples(&T, samples);
    }
    std::set<std::string> not_nearest;
    if (select_nearest > 0) {
        fprintf(stderr, "Selecting context samples...\n");
        //store the reason they were included as a map to add to metadata later.
        //in this case, just store the original list.
        not_nearest.insert(samples.begin(), samples.end());
        std::set<std::string> nsamples;
        for (auto s: samples) {
            auto nearest = get_nearby(&T, s, select_nearest);
            nsamples.insert(nearest.begin(), nearest.end());
        }
        samples.assign(nsamples.begin(), nsamples.end());
    }
    std::set<std::string> not_random;
    if ((setsize > 0) || (add_random > 0)) {
        not_random.insert(samples.begin(), samples.end());
        size_t tv;
        if (setsize > 0) {
            tv = setsize;
        } else {
            fprintf(stderr, "Adding background samples...\n");
            tv = add_random + samples.size();
        }
        samples = fill_random_samples(&T, samples, tv, limit_lca);
    }
    //the final step of selection is to invert the set if prune is set
    //this is done by getting all sample names which are not in the samples vector.
    if (prune_samples) {
        fprintf(stderr, "Sample pruning requested...\n");
        std::vector<std::string> nsamples;
        std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
        for (auto s: T.get_leaves_ids()) {
            //for every sample in the tree, if that sample is NOT in the selected set, save it
            if (sample_set.find(s) == sample_set.end()) {
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
    //make sure all whitelisted samples are included in the output after all other selection is performed.
    if (whitelist_samples_file != "") {
        fprintf(stderr, "Whitelisting samples...\n");
        auto wsamples = read_sample_names(whitelist_samples_file);
        for (auto w: wsamples) {
            if (std::find(samples.begin(),samples.end(),w) == samples.end()) {
                samples.push_back(w);
            }
        }
    }
    //before performing any other action, reroot the tree if requested.
    if (reroot_node != "") {
        reroot_tree(&T, reroot_node);
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    //retrive path information for samples, clades, everything before pruning occurs. Behavioral change
    //to get the paths post-pruning, will need to save a new tree .pb and then repeat the extract command on that
    if (sample_path_filename != dir_prefix) {
        timer.Start();
        if (sample_path_filename != dir_prefix) {
            std::cerr << "Sample mutation path information requested; writing pre-pruning paths to " << sample_path_filename << " with usher-style output.\n";
            if (samples.size() > 0) {
                MAT::get_sample_mutation_paths(&T, samples, sample_path_filename);
            } else {
                fprintf(stderr, "No samples selected; writing paths for all samples...\n");
                MAT::get_sample_mutation_paths(&T, T.get_leaves_ids(), sample_path_filename);
            }
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //before proceeding to actually applying sample selection,
    //if usher-style subtree output is requested,
    //produce that.

    if (usher_minimum_subtrees_size > 0) {
        timer.Start();
        fprintf(stderr, "Random minimum sample subtrees of size %ld requested.\n", usher_minimum_subtrees_size);
        if (samples.size() > 0) {
            MAT::get_random_sample_subtrees(&T, samples, dir_prefix, usher_minimum_subtrees_size, 0, false, retain_branch);
        } else {
            fprintf(stderr, "ERROR: Minimum sample subtree output requested with no valid samples! Check selection parameters\n");
            exit(1);
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    if (usher_single_subtree_size > 0) {
        timer.Start();
        fprintf(stderr, "Random single encompassing subtree of size %ld requested.\n", usher_single_subtree_size);
        if (samples.size() > 0) {
            MAT::get_random_single_subtree(&T, samples, dir_prefix, usher_single_subtree_size, 0, false, retain_branch);
        } else {
            fprintf(stderr, "ERROR: Encompassing subtree output requested with no valid samples! Check selection parameters\n");
            exit(1);
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    size_t num_annotations = T.get_num_annotations();
    if (usher_clades_txt && (usher_minimum_subtrees_size > 0 || usher_single_subtree_size > 0) && num_annotations > 0) {
        // Also make a clades.txt file for the samples like usher does.
        timer.Start();
        auto annotations_filename = dir_prefix + "/clades.txt";
        fprintf(stderr, "Writing clade annotations to file %s \n", annotations_filename.c_str());
        FILE* annotations_file = fopen(annotations_filename.c_str(), "w");
        if (! annotations_file) {
            fprintf(stderr, "ERROR: could not open file %s\n", annotations_filename.c_str());
            exit(1);
        }
        for (std::string sample: samples) {
            fprintf(annotations_file, "%s", sample.c_str());
            for (size_t k=0; k < num_annotations; k++) {
                MAT::Node* node = T.get_node(sample);
                if (node == NULL) {
                    fprintf(stderr, "ERROR: no node for sample %s\n", sample.c_str());
                    exit(1);
                }
                std::string clade = T.get_clade_assignment(node, k, false);
                fprintf(annotations_file, "\t%s", clade.c_str());
            }
            fprintf(annotations_file, "\n");
        }
        fclose(annotations_file);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    //if we have no samples at the end without throwing an error,
    //probably because no selection arguments were set,
    //we want the samples to be all samples in the tree
    //and don't bother pruning.
    MAT::Tree subtree;
    if (samples.size() == 0) {
        timer.Start();
        fprintf(stderr, "No sample selection arguments passed; using full input tree for further output.\n");
        samples = T.get_leaves_ids();
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
        //if no selection was set, then there's no need to filter.
        //subtree is identical to tree.
        subtree = T;
    } else {
        //second step is to filter the input based on the samples

        //TODO: filter_master can support pruning directly rather than generating an inverse set and pruning all but,
        //but downstream stuff wants the inverse set of sample names to work with
        //so there's a better way to do this in at least some cases.
        subtree = filter_master(T, samples, false, true);
    }
    //if polytomy resolution was requested, apply it
    if (resolve_polytomies) {
        timer.Start();
        fprintf(stderr, "Resolving polytomies...\n");
        subtree = resolve_all_polytomies(subtree);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //getting clade representative samples is a special additional step after constructing a subtree
    //so that it selects the first two nodes which fulfilled all other conditions and were retained in the initial filtering
    //TODO: there should be a better way to integrate this information that doesn't involve multiple steps of filtering
    //its basically just that the selection of two samples is very dependent on the other samples which were valid
    //add a valid samples vector argument to the clade representatives function and check membership before selection, maybe
    if (get_representative > 0) {
        if (get_representative < 2) {
            fprintf(stderr, "ERROR: value of --get-representative must be at least 2.\n");
            exit(1);
        }
        //fprintf(stderr, "Filtering again to a clade representative tree...\n");
        auto rep_samples = get_clade_representatives(&subtree, get_representative);
        //run filter master again
        subtree = filter_master(subtree, rep_samples, false, true);
        //overwrite samples with new subset
        samples = rep_samples;
    }
    //if additional information was requested, save it to the target files
    //starting with dumping the set of samples used to a plain text file
    //this is intended for use with uncertainty and perhaps repeated extract commands
    if (used_sample_filename != dir_prefix) {
        //not bothering timing it because this will always be very fast since its an extremely simple non-MAT operation
        fprintf(stderr, "Dumping selected samples to file...\n");
        std::ofstream outfile (used_sample_filename);
        for (auto s: samples) {
            outfile << s << "\n";
        }
    }
    if (all_path_filename != dir_prefix) {
        //print all node mutations *for this subtree*
        timer.Start();
        fprintf(stderr, "Writing full node mutation information to %s\n", all_path_filename.c_str());
        std::ofstream outfile (all_path_filename);
        auto apaths = all_nodes_paths(&T);
        for (auto astr: apaths) {
            outfile << astr << "\n";
        }
        outfile.close();
        fprintf(stderr, "Completed in %ld msec\n\n", timer.Stop());
    }
    if (clade_path_filename != dir_prefix) {
        timer.Start();
        fprintf(stderr, "Writing clade paths to %s\n", clade_path_filename.c_str());
        std::ofstream outfile (clade_path_filename);
        auto cpaths = clade_paths(&T);
        for (auto cstr: cpaths) {
            outfile << cstr;
        }
        outfile.close();
        fprintf(stderr, "Completed in %ld msec\n\n", timer.Stop());
    }

    std::vector<std::unordered_map<std::string,std::unordered_map<std::string,std::string>>> catmeta;
    std::vector<std::string> metav;
    if (meta_filename != "") {
        std::stringstream mns(meta_filename);
        std::string m;
        while (std::getline(mns,m,',')) {
            metav.push_back(m);
        }
        assert (metav.size() > 0);
        std::set<std::string> samples_included(samples.begin(), samples.end());
        if (output_tax_filename == dir_prefix) {
            // Don't do this in the case of taxodium pb output
            for (auto mv: metav) {
                auto scm = read_metafile(mv, samples_included);
                catmeta.emplace_back(scm);
            }
        }
    }
    //if json output AND mutation context is requested, add an additional metadata column indicating whether each branch contains
    //the mutation of interest. the metadata map is not limited to leaf nodes.
    if ((json_filename != "") && (mutation_choice != "")) {
        std::unordered_map<std::string,std::string> mutmap;
        std::vector<std::string> mutations;
        std::stringstream mns(mutation_choice);
        std::string m;
        while (std::getline(mns,m,',')) {
            mutations.push_back(m);
        }
        assert (mutations.size() > 0);
        for (auto n: subtree.depth_first_expansion()) {
            std::string metastr = "";
            for (auto m: n->mutations) {
                for (auto mstr: mutations ) {
                    if (m.get_string() == mstr) {
                        if (metastr == "") {
                            metastr = mstr;
                        } else {
                            metastr = metastr + ',' + mstr;
                        }
                    }
                }
            }
            mutmap[n->identifier] = metastr;
        }
        std::unordered_map<std::string,std::unordered_map<std::string,std::string>> submet;
        submet["mutation_of_interest"]=mutmap;
        catmeta.push_back(submet);
    }
    if (minimum_subtrees_size > 0) {
        fprintf(stderr, "Finding minimum covering sample subtrees.\n");
        timer.Start();
        //references the original tree for getting nearest background.
        get_minimum_subtrees(&T, samples, minimum_subtrees_size, dir_prefix, &catmeta, json_filename, tree_filename, retain_branch);
        fprintf(stderr, "Minimum subtree files written in %ld msec; exiting\n", timer.Stop());
        exit(0);
    }
    if (nearest_k_batch_file != "") {
        fprintf(stderr, "Batch sample context writing requested.\n");
        auto split_point = nearest_k_batch_file.find(":");
        if (split_point == std::string::npos) {
            fprintf(stderr, "ERROR: Invalid formatting of -K argument. Requires input in the form of 'sample_file.txt:k' to generate json context files\n");
            exit(1);
        }
        std::string sample_file = nearest_k_batch_file.substr(0, split_point);
        std::string nkstr = nearest_k_batch_file.substr(split_point+1, nearest_k_batch_file.size() - split_point);
        int nk = std::stoi(nkstr);
        if (nk <= 0) {
            fprintf(stderr, "ERROR: Invalid neighborhood size. Please choose a positive nonzero integer.\n");
            exit(1);
        }
        auto batch_samples = read_sample_names(sample_file);
        timer.Start();
        static tbb::affinity_partitioner ap;

        tbb::parallel_for(tbb::blocked_range<size_t>(0, batch_samples.size() ),
        [&](const tbb::blocked_range<size_t> r) {

            for (auto s = r.begin() ; s < r.end() ; ++s ) {
                auto cs = get_nearby(&T, batch_samples[s], nk);
                if ( cs.size() == 0 ) {
                    continue ;
                }
                MAT::Tree subt = filter_master(T, cs, false);
                //remove forward slashes from the string, replacing them with underscores, so that it works for file names
                size_t pos = 0;
                while ((pos = batch_samples[s].find("/")) != std::string::npos) {
                    batch_samples[s].replace(pos, 1, "_");
                }
                write_json_from_mat(&subt, batch_samples[s] + "_context.json", &catmeta, tax_title);
            }
        }, ap) ;
        fprintf(stderr, "%ld batch sample jsons written in %ld msec.\n\n", batch_samples.size(), timer.Stop());
    }
    if (closest_relatives_filename != dir_prefix) {
        fprintf(stderr, "Per-sample closest relative(s) requested. Computing...\n");
        if (break_ties) {
            fprintf(stderr, "Storing one closest relative per sample.\n");
        }
        std::ofstream out(closest_relatives_filename);

        timer.Start();
        std::vector<std::string> closest_relatives;

        for (std::string sample : samples) {
            std::pair<std::vector<std::string>, size_t> closest_relatives_pair = get_closest_samples(&T, sample);
            std::vector<std::string> closest_relatives = closest_relatives_pair.first;
            size_t dist = closest_relatives_pair.second;
            if (closest_relatives.size() > 0) {
                closest_relatives = closest_relatives_pair.first;
                std::string s = "";
                s += sample + '\t';
                std::string lex_smallest_sample = closest_relatives[0];
                for (std::string relative : closest_relatives) {
                    if (break_ties) {
                        if (relative < lex_smallest_sample) {
                            lex_smallest_sample = relative;
                        }
                    } else {
                        s += relative + ',';
                    }
                }
                if (break_ties) {
                    s += lex_smallest_sample;
                } else {
                    s = s.substr(0, s.size() - 1);
                }
                s += '\t' + std::to_string(dist);
                out << s << "\n";
            }
        }
        fprintf(stderr, "TSV of closest relative written to %s in %ld msec.\n\n", closest_relatives_filename.c_str(), timer.Stop());
    }
    //if json output AND sample context is requested, add an additional metadata column which simply indicates the focal sample versus context
    if ((json_filename != "") && (nearest_k != "")) {
        std::unordered_map<std::string,std::string> conmap;
        for (auto s: samples) {
            if (s == sample_id) {
                conmap[s] = "focal";
            }
        }
        std::unordered_map<std::string,std::unordered_map<std::string,std::string>> submet;
        submet["focal_view"] = conmap;
        catmeta.emplace_back(submet);
    }
    //same for random samples and nearest samples.
    if ((json_filename != "") && ((not_random.size() > 0) || (not_nearest.size() > 0))) {
        std::unordered_map<std::string,std::string> ranmap;
        for (auto s: samples) {
            //samples that are random will not appear in either set
            //samples that are nearest will appear in not_random but not not_nearest
            if ((not_random.size() > 0) && (not_random.find(s) == not_random.end())) {
                ranmap[s] = "random";
            } else if ((not_nearest.size() > 0) && (not_nearest.find(s) == not_nearest.end())) {
                ranmap[s] = "nearest";
            } else {
                ranmap[s] = "primary";
            }
        }
        std::unordered_map<std::string,std::unordered_map<std::string,std::string>> submet;
        submet["selected_for"] = ranmap;
        catmeta.emplace_back(submet);
    }

    //last step is to convert the subtree to other file formats
    if (vcf_filename != dir_prefix) {
        fprintf(stderr, "Generating VCF of final tree\n");
        make_vcf(subtree, vcf_filename, no_genotypes, samples);
    }
    if (json_filename != dir_prefix) {
        fprintf(stderr, "Generating JSON of final tree\n");
        write_json_from_mat(&subtree, json_filename, &catmeta, tax_title);
    }
    if (tree_filename != dir_prefix) {
        fprintf(stderr, "Generating Newick file of final tree\n");
        FILE *tree_file = fopen(tree_filename.c_str(), "w");
        fprintf(tree_file, "%s\n",
                MAT::get_newick_string(subtree, true, true, retain_branch).c_str());
        fclose(tree_file);
    }
    //and save a MAT if that was set
    if (output_mat_filename != dir_prefix) {
        timer.Start();
        fprintf(stderr, "Saving output MAT file %s.\n", output_mat_filename.c_str());
        //only recondense the tree if polytomies weren't resolved.
        if (collapse_tree) {
            subtree.collapse_tree();
        }
        if (!resolve_polytomies) {
            subtree.condense_leaves();
        }

        MAT::save_mutation_annotated_tree(subtree, output_mat_filename);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    if (output_tax_filename != dir_prefix) {
        timer.Start();
        bool quit = false;
        if (gtf_filename == dir_prefix) {
            fprintf(stderr, "ERROR: You must specify a GTF file with -g\n");
            quit = true;
        }
        if (fasta_filename == dir_prefix) {
            fprintf(stderr, "ERROR: You must specify a FASTA reference file with -f\n");
            quit = true;
        }
        if (quit) {
            exit(1);
        }
        fprintf(stderr, "Saving output MAT file in Taxodium format: %s.\n",  output_tax_filename.c_str());
        if (collapse_tree) {
            subtree.collapse_tree();
        }
        if (!resolve_polytomies) {
            subtree.condense_leaves();
        }
        save_taxodium_tree(subtree, output_tax_filename, metav, gtf_filename, fasta_filename, tax_title, tax_description, additional_meta_fields, x_scale, include_nt);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    if (dump_metadata != dir_prefix) {
        timer.Start();
        fprintf(stderr, "Dumping final metadata.\n");
        std::ofstream mout (dump_metadata);
        mout << "strain";
        //first, all column names (except for sample, which is always first).
        //ordered map so iteration is consistent with column names.
        std::map<std::string,size_t> column_order;
        std::vector<std::string> row_order;
        size_t i = 0;
        for (auto& omap: catmeta) {
            for (auto& cv: omap) {
                column_order[cv.first] = i;
                mout << "\t" << cv.first;
            }
            i++;
        }
        for (auto s: samples) {
            mout << "\n" << s;
            for (auto cv: column_order) {
                auto sv = catmeta[cv.second][cv.first].find(s);
                if (sv != catmeta[cv.second][cv.first].end()) {
                    mout << "\t" << sv->second;
                } else {
                    mout << "\tmissing";
                }
            }
        }
        mout << "\n";
        mout.close();
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
}
