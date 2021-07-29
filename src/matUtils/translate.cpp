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

    MAT::Tree T = MAT::load_mutation_annotated_tree("public-latest.all.masked.pb");
    //T here is the actual object.

    std::map<int, std::vector<Codon *>> codon_map;
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
        
            for (auto & c: fasta_line) c = (char)toupper(c);
        
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
//            std::cout << feature << ',' << attribute << ',' << start << ',' << stop << '\n';
            for (int pos = start - 1; pos < stop; pos += 3) {
                
                char nt[3] = {
                    reference[pos],
                    reference[pos+1],
                    reference[pos+2]
                };

                Codon *c = new Codon(attribute, (pos - start + 1) / 3,  pos, nt);
                
                auto it = codon_map.find(pos);
                if (it == codon_map.end()) {
                    codon_map.insert({pos, {c}});
                } else { 
                    (it->second).push_back(c);
                }

                it = codon_map.find(pos+1);
                if (it == codon_map.end()) {
                    codon_map.insert({pos+1, {c}});
                } else {
                    (it->second).push_back(c);
                }

                it = codon_map.find(pos+2);
                if (it == codon_map.end()) {
                    codon_map.insert({pos+2, {c}});
                } else {
                    (it->second).push_back(c);
                }
             }
        }       
    }


    // std::string orf1a = "";
    // std::string orf1b = "";
    // std::string orf7a = "";
    // std::string orf7b = "";
    // std::string S = "";
    

//    validate proteins constructed correctly
//     for (auto const& [key, val] : codon_map) {
// //        std::cout << '\n' << key << "=>";
//         for (auto v : val) {
//             if (v->orf_name == " gene_name \"ORF1a\"") {
//                 orf1a += v->protein;
//             } else if (v->orf_name == " gene_name \"ORF1b\"") {
//                 orf1b += v->protein;
//             } else if (v->orf_name == " gene_name \"ORF7a\"") {
//                 orf7a += v->protein;
//             } else if (v->orf_name == " gene_name \"ORF7b\"") {
//                 orf7b += v->protein;
//             } else if (v->orf_name == " gene_name \"S\"") {
//                 S += v->protein;
//             }
//         }
//     }
    


    auto dfs = T.depth_first_expansion();
    int count = 0;
    MAT::Node *last_visited = nullptr;
    bool first = true;
    for (auto node: dfs) {
        std::string mutation_result = "";
        if (first) {
//            std::cout << "\nBeginning DFS for translate\n";
 //           std::cout << "NODE: root\n";
            first = false;
            std::cout << node->identifier << "\t.\t.\n";

        } else {

            
            //std::cout << "\n--------\nNODE: " << node->identifier << '\n';
            if (last_visited == node->parent) {
                ;
            //    std::cout << "\nvisiting a child\n";
            } else {
                // Jumping across a branch, so we need to revert codon mutations up to
                // the LCA of this node and the last visited
                MAT::Node *last_common_ancestor = MAT::LCA(T, node->identifier, last_visited->identifier);
                MAT::Node *trace_to_lca = last_visited;
            //    std::cout << "retracing..." << '\n';
                while (trace_to_lca != last_common_ancestor) {
                    undo_mutations(trace_to_lca->mutations, codon_map);        
                    trace_to_lca = trace_to_lca->parent;
                }
           //     std::cout << "reached LCA"  << '\n';
            }
          //  std::cout << "\ndoing mutations:" << '\n';
            mutation_result = do_mutations(node->mutations, codon_map);

            count += 1;
         //   if (count > 20) {
         //      break;
         //   }

         std::cout << node->identifier << '\t' << mutation_result;

        }
        

        last_visited = node;

    }
}



std::string do_mutations(std::vector<MAT::Mutation> &mutations, std::map<int, std::vector<Codon *>> &codon_map) {

    std::string prot_string = "";
    std::string nuc_string = "";
    for (auto m: mutations) {
        nuc_string += m.get_string();
        nuc_string += ',';
        char mutated_nuc = MAT::get_nuc(m.mut_nuc);
        int pos = m.position;
        auto it = codon_map.find(pos);
        if (it == codon_map.end()) {
            continue; // Not a coding mutation
        } else {
            // Mutate each codon associated with this position
            for (auto codon_ptr : it->second) {
                prot_string += codon_ptr->orf_name + ":";
                prot_string += codon_ptr->protein;
                codon_ptr->mutate(pos, mutated_nuc);
                prot_string += std::to_string(codon_ptr->codon_number);
                prot_string += codon_ptr->protein;
                prot_string += ',';                
             }

        }
    }
    if (!prot_string.empty() && prot_string.back() == ',') {
        prot_string.resize(prot_string.length() - 1); //remove trailing ',' 
    }
    if (!nuc_string.empty() && nuc_string.back() == ',') {
        nuc_string.resize(nuc_string.length() - 1);
    } else if (nuc_string.empty()) {
        nuc_string = ".";
    }
    if (prot_string.empty()) {
        prot_string = ".";
    }

    return prot_string + '\t' + nuc_string + '\n';
}            

void undo_mutations(std::vector<MAT::Mutation> &mutations, std::map<int, std::vector<Codon *>> &codon_map) {
    for (auto m: mutations) {
        char parent_nuc = MAT::get_nuc(m.par_nuc);
        int pos = m.position;
        auto it = codon_map.find(pos);
        if (it == codon_map.end()) {
            continue;
            // Not a coding mutation
        } else {
            // Revert the mutation by mutating to the parent nucleotide
            for (auto codon_ptr : it->second) {
                codon_ptr->mutate(pos, parent_nuc);
            }
        }
    //    std::cout << "undoing " << m.get_string() << '\n';
    }
}
