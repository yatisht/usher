#include "translate.hpp"

po::variables_map parse_translate_command(po::parsed_options parsed) {

    po::variables_map vm;
    po::options_description filt_desc("introduce options");
    filt_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("population-samples,s", po::value<std::string>()->required(), 
         "Names of samples from the population of interest [REQUIRED].") 
        ("additional-info,a", po::bool_switch(),
        "Set to calculate additional phylogenetic trait association statistics for whole regions and individual introductions. WARNING: Adds significantly to runtime.")
        ("clade-regions,c", po::value<std::string>()->default_value(""),
        "Set to optionally record, for each clade root in the tree, the support for that clade root being IN each region in the input, as a tsv with the indicated name.")
        ("date-metadata,M", po::value<std::string>()->default_value(""),
        "Pass a TSV or CSV containing a 'date' column to use for date information. If not used, date will be inferred from the sample name where possible.")
        ("output,o", po::value<std::string>()->required(),
        "Name of the file to save the introduction information to.")
        ("origin-confidence,C", po::value<float>()->default_value(0.5),
        "Set the threshold for recording of putative origins of introductions. Default is 0.5")
        ("evaluate-metadata,E", po::bool_switch(),
        "Set to assign each leaf a confidence value based on ancestor distance and confidence.")
        ("dump-assignments,D", po::value<std::string>()->default_value(""),
        "Indicate a directory to which two-column text files containing node assignment values should be dumped for downstream processing.")
        // ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << filt_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void translate_main(po::parsed_options parsed) {
 //   po::variables_map vm = parse_introduce_command(parsed);
 //   std::string input_mat_filename = vm["input-mat"].as<std::string>();
//    std::string samples_filename = vm["population-samples"].as<std::string>();
//    std::string clade_regions = vm["clade-regions"].as<std::string>();
//    std::string metafile = vm["date-metadata"].as<std::string>();
//    bool add_info = vm["additional-info"].as<bool>();
  //  std::string output_file = vm["output"].as<std::string>();
//    std::string dump_assignments = vm["dump-assignments"].as<std::string>();
//    float moconf = vm["origin-confidence"].as<float>();
//    bool leafconf = vm["evaluate-metadata"].as<bool>();
    // int32_t num_threads = vm["threads"].as<uint32_t>();

    // Load input MAT and uncondense tree

    MAT::Tree T = MAT::load_mutation_annotated_tree("test.pb.gz");
    //T here is the actual object.

    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }


    std::ifstream infile("test.gff");
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the GFF file: %s!\n", "filename");
        exit(1);
    }
    std::string gff_line;

    while (std::getline(infile, gff_line)) {
        std::cout << gff_line;
    };


    std::vector<Codon *> codonPtrs = {
        new Codon("Orf1", 0, 'A', 'T', 'G'),
        new Codon("Orf1", 3, 'T', 'T', 'C'),
        new Codon("Orf1", 6, 'C', 'C', 'C'),
        new Codon("Orf1", 5, 'C', 'C', 'C')
    };

    std::map<int, std::vector<Codon *>> codons = {
        {0, {codonPtrs[0]}},
        {1, {codonPtrs[0]}},
        {2, {codonPtrs[0]}},
        {3, {codonPtrs[1]}},
        {4, {codonPtrs[1]}},
        {5, {codonPtrs[1], codonPtrs[3]}},
        {6, {codonPtrs[2], codonPtrs[3]}},
        {7, {codonPtrs[2], codonPtrs[3]}},
        {8, {codonPtrs[2]}}  
    };

    // for (auto c : codons) {
    //     std::cout << c.get_string() << '\n';
    // }
    auto dfs = T.depth_first_expansion();
    int count = 0;
    std::string parentNuc;
    std::string mutatedNuc;
    std::string pos;
    for (auto s: dfs) {
        std::cout << '\n' << s->identifier << '\n';
        for (auto m: s->mutations) {
            parentNuc = MAT::get_nuc(m.par_nuc);
            mutatedNuc = MAT::get_nuc(m.mut_nuc);
            pos = m.position;
            std::cout << m.get_string() << '\n';
        }
        count += 1;
        if (count > 10) {
            break;
        }
    }

    for (auto c: codonPtrs) {
        delete c;
    }
}
