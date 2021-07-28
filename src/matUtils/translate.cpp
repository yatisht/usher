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

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> result;
	std::stringstream ss(s);
	std::string item;
	while (getline(ss, item, delim)) {
		result.push_back(item);
	}
	return result;
}

char translate_codon(std::string nt) {
    if (nt == "GCT" || nt == "GCC" || nt == "GCA" || nt == "GCG" || nt == "GCN") {
        return 'A';
    } else if (nt == "TGT" || nt == "TGC" || nt == "TGY") {
            return 'C';
    } else if (nt == "GAT" || nt == "GAC" || nt == "GAY") {
            return 'D';
    } else if (nt == "GAA" || nt == "GAG" || nt == "GAR") {
            return 'E';
    } else if (nt == "TTT" || nt == "TTC" || nt == "TTY") {
            return 'F';
    } else if (nt == "GGT" || nt == "GGC" || nt == "GGA" || nt == "GGG" || nt == "GGN") {
            return 'G';
    } else if (nt == "CAT" || nt == "CAC" || nt == "CAY") {
            return 'H';
    } else if (nt == "ATT" || nt == "ATC" || nt == "ATA" || nt == "ATH") {
            return 'I';
    } else if (nt == "AAA" || nt == "AAG" || nt == "AAR") {
            return 'K';
    } else if (nt == "TTA" || nt == "TTG" || nt == "CTT" || nt == "CTC" || nt == "CTA" || nt == "CTG" || nt == "YTR" || nt == "CTN") {
            return 'L';
    } else if (nt == "ATG") {
            return 'M';
    } else if (nt == "AAT" || nt == "AAC" || nt == "AAY") {
            return 'N';
    } else if (nt == "CCT" || nt == "CCC" || nt == "CCA" || nt == "CCG" || nt == "CCN") {
            return 'P';
    } else if (nt == "CAA" || nt == "CAG" || nt == "CAR") {
            return 'Q';
    } else if (nt == "CGT" || nt == "CGC" || nt == "CGA" || nt == "CGG" || nt == "AGA" || nt == "AGG" || nt == "CGN" || nt == "MGR") {
            return 'R';
    } else if (nt == "TCT" || nt == "TCC" || nt == "TCA" || nt == "TCG" || nt == "AGT" || nt == "AGC" || nt == "TCN" || nt == "AGY") {
            return 'S';
    } else if (nt == "ACT" || nt == "ACC" || nt == "ACA" || nt == "ACG" || nt == "ACN") {
            return 'T';
    } else if (nt == "GTT" || nt == "GTC" || nt == "GTA" || nt == "GTG" || nt == "GTN") {
            return 'V';
    } else if (nt == "TGG") {
            return 'W';
    } else if (nt == "TAT" || nt == "TAC" || nt == "TAY") {
            return 'Y';
    } else if (nt == "TAG" || nt == "TAA" || nt == "TGA") {
        return '*';
    } else { //ambiguous
        return 'X';
    }
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

    std::map<int, std::vector<Codon *>> codonMap;
    std::vector<Codon> codons;

    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }


    std::ifstream fasta_file("ref.fasta");
    if (!fasta_file) {
        fprintf(stderr, "ERROR: Could not open the fasta file: %s!\n", "filename");
        exit(1);
    }
    std::ifstream gff_file("test.gff");
    if (!gff_file) {
        fprintf(stderr, "ERROR: Could not open the GFF file: %s!\n", "filename");
        exit(1);
    }

    std::string fasta_line;
    std::string reference = "";
    size_t line_length;
    while(std::getline(fasta_file, fasta_line)) {
        if (fasta_line[0] == '>' or fasta_line[0] == '\n') {
            continue;
        } else {
            line_length = fasta_line.length();
            if (fasta_line[line_length-1] == '\r') {
                fasta_line.erase(line_length-1);
            }
            reference += fasta_line;
        }
    }

    
    std::string gff_line;
    std::string feature;
    std::string attribute;
    int start;
    int stop;
    while (std::getline(gff_file, gff_line)) {
        if (gff_line[0] == '#' || gff_line[0] == '\n') {
            continue;
        }

        std::vector<std::string> split_line = split(gff_line, '\t');

        if(split_line.size() <= 1) {
            continue;
        }
        
        feature = split_line[2];
        if (feature == "gene") {
            attribute = split_line[8];
            start = std::stoi(split_line[3]);
            stop = std::stoi(split_line[4]);
            std::cout << feature << ',' << attribute << ',' << start << ',' << stop << '\n';
            for (int pos = start - 1; pos < stop; pos += 3) {
                
                std::string nt = "";
                nt +=  reference[pos];
                nt += reference[pos+1];
                nt += reference[pos+2];
                Codon *c = new Codon(
                    attribute,
                    pos,
                    nt,
                    translate_codon(nt)
                );
                
                // pos not in map yet
                if (codonMap.find(pos) == codonMap.end()) {
                    codonMap.insert({pos, {c}});
                } else { // add to existing list
                    codonMap[pos].push_back(c);
                }
                if (codonMap.find(pos+1) == codonMap.end()) {
                    codonMap.insert({pos+1, {c}});
                } else {
                    codonMap[pos+1].push_back(c);
                }
                if (codonMap.find(pos+2) == codonMap.end()) {
                    codonMap.insert({pos+2, {c}});
                } else {
                    codonMap[pos+2].push_back(c);
                }
             }
        }
                
    }
    for (auto const& [key, val] : codonMap) {
        std::cout << '\n' << key << "=>";
        for (auto v : val) {
            std::cout << '\n' << "    " <<v->get_string() << '\n';
        }
    }
}

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

