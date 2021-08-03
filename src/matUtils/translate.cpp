#include "translate.hpp"


po::variables_map parse_translate_command(po::parsed_options parsed) {
  
    po::variables_map vm;
    po::options_description filt_desc("translate options");
    filt_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("output,o", po::value<std::string>()->required(),
        "Name of the file to save the tsv output to [REQUIRED]")
        ("input-fasta,f", po::value<std::string>()->required(),
         "Input reference sequence fasta [REQUIRED]")
        ("input-gff,g", po::value<std::string>()->required(),
        "Input GFF file [REQUIRED]")
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

std::string build_reference(std::ifstream &fasta_file) {
    std::string reference_output = "";
    std::string fasta_line;
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
            reference_output += fasta_line;
        }
    }
    return reference_output;
}

// Maps a genomic coordinate to a list of codons it is part of
std::map<int, std::vector<std::shared_ptr<Codon>>> build_codon_map(std::ifstream &gff_file, std::string reference) {
    std::map<int, std::vector<std::shared_ptr<Codon>>> codon_map;
    std::string gff_line;
    while (std::getline(gff_file, gff_line)) {
        if (gff_line[0] == '#' || gff_line[0] == '\n') {
            continue;
        }
        std::vector<std::string> split_line = split(gff_line, '\t');
        if(split_line.size() <= 1) {
            continue;
        }
        std::string feature = split_line[2];
        if (feature == "gene") {
            std::string attribute = split_line[8];
            int start = std::stoi(split_line[3]);
            int stop = std::stoi(split_line[4]);
            for (int pos = start - 1; pos < stop; pos += 3) {

                char nt[3] = {
                    reference[pos],
                    reference[pos+1],
                    reference[pos+2]
                };

                // Coordinates are 0-based at this point
                std::shared_ptr<Codon> c(new Codon(attribute, ((pos - start + 1) / 3),  pos, nt));
                
                // The current pos and the next positions
                // are associated with this codon
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
    return codon_map;
}


void translate_main(po::parsed_options parsed) {
    po::variables_map vm = parse_translate_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string input_fasta_filename =  vm["input-fasta"].as<std::string>();
    std::string input_gff_filename = vm["input-gff"].as<std::string>();
    std::string output_filename = vm["output"].as<std::string>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);

    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }

    std::ifstream fasta_file(input_fasta_filename);
    if (!fasta_file) {
        fprintf(stderr, "ERROR: Could not open the fasta file: %s!\n", input_fasta_filename.c_str());
        exit(1);
    }
    std::ifstream gff_file(input_gff_filename);
    if (!gff_file) {
        fprintf(stderr, "ERROR: Could not open the GFF file: %s!\n", input_gff_filename.c_str());
        exit(1);
    }
    std::ofstream output_file(output_filename);
    if (!output_file) {
        fprintf(stderr, "ERROR: Could not open file for writing: %s!\n", output_filename.c_str());
        exit(1);
    }
    
    std::string reference = build_reference(fasta_file);

    output_file << "Node ID\tAmino acid mutations\tNucleotide mutations\tLeaves sharing mutations" << '\n';

    // This maps each position in the reference to a vector of codons.
    // Some positions may be associated with multiple codons (frame shifts).
    // The Codons in the map are updated as the tree is traversed
    std::map<int, std::vector<std::shared_ptr<Codon>>> codon_map = build_codon_map(gff_file, reference);

    // Traverse the tree in depth-first order. As we descend the tree, mutations at
    // each node are applied to the respective codon(s) in codon_map.
    auto dfs = T.depth_first_expansion();
    MAT::Node *last_visited = nullptr;
    for (auto node: dfs) {
        std::string mutation_result = "";
        if(last_visited != node->parent) {
            // Jumping across a branch, so we need to revert codon mutations up to
            // the LCA of this node and the last visited node
            MAT::Node *last_common_ancestor = MAT::LCA(T, node->identifier, last_visited->identifier);
            MAT::Node *trace_to_lca = last_visited;
            while (trace_to_lca != last_common_ancestor) {
                undo_mutations(trace_to_lca->mutations, codon_map);        
                trace_to_lca = trace_to_lca->parent;
            }
        } // If we are visiting a child, we can continue mutating
        
        mutation_result = do_mutations(node->mutations, codon_map);
        if (mutation_result != ""){
            output_file << node->identifier << '\t' << mutation_result << '\t' << T.get_leaves(node->identifier).size() << '\n';
        }
        last_visited = node;
    }

}


std::string do_mutations(std::vector<MAT::Mutation> &mutations, std::map<int, std::vector<std::shared_ptr<Codon>>> &codon_map) {

    std::string prot_string = "";
    std::string nuc_string = "";
    for (auto m: mutations) {
        nuc_string += m.get_string();
        nuc_string += ';';
        char mutated_nuc = MAT::get_nuc(m.mut_nuc);
        int pos = m.position - 1;
        auto it = codon_map.find(pos);
        if (it == codon_map.end()) {
            continue; // Not a coding mutation
        } else {
            // Mutate each codon associated with this position
            for (auto codon_ptr : it->second) {
                prot_string += codon_ptr->orf_name + ":";
                prot_string += codon_ptr->protein;
                
                codon_ptr->mutate(pos, mutated_nuc);

                prot_string += std::to_string(codon_ptr->codon_number+1);
                prot_string += codon_ptr->protein;
                prot_string += ',';
             }
            prot_string.resize(prot_string.length() - 1);
            prot_string += ';';
        }
    }
    if (!prot_string.empty() && prot_string.back() == ';') {
        prot_string.resize(prot_string.length() - 1); //remove trailing ',' 
    }
    if (!nuc_string.empty() && nuc_string.back() == ';') {
        nuc_string.resize(nuc_string.length() - 1);
    }
    
    if (nuc_string.empty()) {
        return "";
    } else {
        return prot_string + '\t' + nuc_string;
    } 
}            

void undo_mutations(std::vector<MAT::Mutation> &mutations, std::map<int, std::vector<std::shared_ptr<Codon>>> &codon_map) {
    for (auto m: mutations) {
        char parent_nuc = MAT::get_nuc(m.par_nuc);
        int pos = m.position - 1;
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
    }
}