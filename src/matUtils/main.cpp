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

    po::variables_map vm;
    po::options_description conv_desc("convert options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("samples,s", po::value<std::string>()->default_value(""),
        "Select samples by explicitly naming them. One per line")
        ("clade,c", po::value<std::string>()->default_value(""),
        "Select samples by membership in the indicated clade.")
        ("mutation,m", po::value<std::string>()->default_value(""),
        "Select samples by whether they contain the indicated mutation.")
        ("max-epps,e", po::value<size_t>()->default_value(0),
        "Select samples by whether they have less than the maximum indicated number of equally parsimonious placements. Note: calculation adds significantly to runtime.")
        ("prune,p", po::bool_switch(),
        "Remove the selected samples instead of keeping them in the output files.")
        ("sample-paths", po::value<std::string>()->default_value(""),
        "Write the path of mutations defining each selected sample to the target file (all samples if no selection arguments)")
        ("clade-paths", po::value<std::string>()->default_value(""),
        "Write the path of mutations defining each clade in the tree after sample selection to the target file.")
        ("write-vcf,v", po::value<std::string>()->default_value(""),
         "Output VCF file representing selected subtree. Default is full tree")
        ("no-genotypes,n", po::bool_switch(),
        "Do not include sample genotype columns in VCF output. Used only with the vcf option")
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
    size_t max_epps = vm["max-epps"].as<size_t>();
    bool prune_samples = vm["prune"].as<bool>();
    std::string sample_path_filename = vm["sample-paths"].as<std::string>();
    std::string clade_path_filename = vm["clade-paths"].as<std::string>();
    std::string tree_filename = vm["write-tree"].as<std::string>();
    std::string vcf_filename = vm["write-vcf"].as<std::string>();
    std::string output_mat_filename = vm["write-mat"].as<std::string>();
    bool no_genotypes = vm["no-genotypes"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

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
    fprintf(stderr, "Parsing and applying sample selection arguments\n");
    //TODO: the intersection code is repetitive. There's probably some way to restructure it to reuse
    //the code blocks for checking intersection and whether the sample selection is empty after intersection
    std::vector<std::string> samples;
    if (input_samples_file != "") {
        samples = read_sample_names(input_samples_file);
    } 
    if (clade_choice != "") {
        if (samples.size() == 0) {
            samples = get_clade_samples(T, clade_choice);
        } else {
            auto cladesamples = get_clade_samples(T, clade_choice);
            //only retain samples that are in both the current and new sets.
            std::vector<std::string> nsamples;
            for (auto s: samples) {
                if (std::find(cladesamples.begin(), cladesamples.end(), s) != cladesamples.end()) {
                    nsamples.push_back(s);
                }
            }
            samples = nsamples;
            //check to make sure we haven't emptied our sample file; if we have, throw an error
            if (samples.size() == 0) {
                fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
                exit(1);
            }
        }
    }
    if (mutation_choice != "") {
        if (samples.size() == 0) {
            samples = get_mutation_samples(T, mutation_choice);
        } else {
            auto mutsamples = get_mutation_samples(T, mutation_choice);
            //only retain samples that are in both the current and new sets.
            std::vector<std::string> nsamples;
            for (auto s: samples) {
                if (std::find(mutsamples.begin(), mutsamples.end(), s) != mutsamples.end()) {
                    nsamples.push_back(s);
                }
            }
            samples = nsamples;
            //check to make sure we haven't emptied our sample file; if we have, throw an error
            if (samples.size() == 0) {
                fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
                exit(1);
            }
        }
    }
    if (max_epps > 0) {
        if (samples.size() == 0) {
            samples = get_samples_epps(T, max_epps);
        } else {
            auto esamples = get_samples_epps(T, max_epps);
            //only retain samples that are in both the current and new sets.
            std::vector<std::string> nsamples;
            for (auto s: samples) {
                if (std::find(esamples.begin(), esamples.end(), s) != esamples.end()) {
                    nsamples.push_back(s);
                }
            }
            samples = nsamples;
            //check to make sure we haven't emptied our sample file; if we have, throw an error
            if (samples.size() == 0) {
                fprintf(stderr, "ERROR: No samples fulfill selected criteria. Change arguments and try again\n");
                exit(1);
            }
        }
    }
    //the final step of selection is to invert the set if prune is set
    //this is done by getting all sample names which are not in the samples vector.
    if (prune_samples) {
        std::vector<std::string> nsamples;
        for (auto s: T.get_leaves_ids()) {
            //for every sample in the tree, if that sample is NOT in the selected set, save it
            if (std::find(samples.begin(), samples.end(), s) == samples.end()) {
                nsamples.push_back(s);
            }
        }
        //and overwrite our choice of samples
        samples = nsamples;
        //check to make sure we haven't emptied our sample file; if we have, throw an error
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
        timer.Start();
        fprintf(stderr, "Extracting subtree...\n");
        //second step is to filter the input based on the samples
        //TODO: filter_master supports pruning directly rather than generating an inverse set,
        //but downstream stuff wants the inverse set of sample names to work with
        //so there's a better way to do this probably.
        subtree = filter_master(T, samples, false);
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }

    //if additional information was requested, save it to the target files
    if (sample_path_filename != "" || clade_path_filename != "") {
        timer.Start();
        fprintf(stderr,"Retriving path information...\n");
        if (sample_path_filename != "") {
            std::ofstream outfile (sample_path_filename);
            auto mpaths = mutation_paths(subtree, samples);
            for (auto mstr: mpaths) {
                outfile << mstr << "\n";
            }
            outfile.close();
        }
        if (clade_path_filename != "") {
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
                            if (c != "" && std::find(cladenames.begin(), cladenames.end(), c) != cladenames.end()) {
                                cladenames.push_back(c);
                            }
                        }
                    }
                }
            }
            std::ofstream outfile (clade_path_filename);
            //maybe better to update clade_paths to take an "all current clades" option.
            //should be more computationally efficient at least.
            auto cpaths = clade_paths(subtree, cladenames);
            for (auto cstr: cpaths) {
                outfile << cstr << "\n";
            }
            outfile.close(); 
        }
        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    }
    //TODO: add argument to save all output in a specified directory.

    //last step is to convert the subtree to other file formats
    if (vcf_filename != "") {
        fprintf(stderr, "Generating VCF of final tree\n");
        make_vcf(subtree, vcf_filename, no_genotypes);
    }
    if (tree_filename != "") {
        fprintf(stderr, "Generating Newick file of final tree\n");
        auto tree_filepath = "./" + tree_filename; //for simplicity, write it to this directory
        FILE *tree_file = fopen(tree_filepath.c_str(), "w");
        fprintf(tree_file, "%s\n",
            MAT::get_newick_string(subtree, true, true, true).c_str());
        fclose(tree_file);        
    }
    //and save a MAT if that was set
    if (output_mat_filename != "") {
        fprintf(stderr, "Saving output MAT file %s.\n", output_mat_filename.c_str()); 
        subtree.condense_leaves();
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
    
    try {
        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();

        po::store(parsed, vm);
        std::string cmd = vm["command"].as<std::string>();
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
        } else if (cmd == "help" || cmd == "--help" || cmd == "-h") { 
            // TODO: improve this message
            fprintf(stderr, "matUtils has several major subcommands: annotate, mask, extract, summary, and uncertainty.\nIndividual command options can be accessed with matUtils command --help, e.g. matUtils annotate --help will show annotation-specific help messages.");
            exit(0);
        } else {
            fprintf(stderr, "Invalid command. Please choose from annotate, mask, extract, summary, uncertainty, or help and try again.\n");
            exit(1);
        }
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Please choose from annotate, mask, extract, summary, uncertainty, or help and try again.\n");
        exit(0);
    }

    return 0;
}

