#include "translate.hpp" 
 
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
std::map<int, std::vector<std::shared_ptr<Codon>>> build_codon_map(std::ifstream &gtf_file, std::string reference) { 
    std::map<int, std::vector<std::shared_ptr<Codon>>> codon_map; 
    std::string gtf_line; 
    std::vector<std::string> gtf_lines; 
    std::vector<std::string> done; 
    while (std::getline(gtf_file, gtf_line)) { 
        gtf_lines.push_back(gtf_line); 
    } 
    for (std::string line_outer : gtf_lines) { 
        if (line_outer[0] == '#' || line_outer[0] == '\n') { 
            continue; 
        } 
        std::vector<std::string> split_line_outer = split(line_outer, '\t'); 
        if(split_line_outer.size() <= 1) { 
            continue; 
        } 
        if (split_line_outer[8].substr(0, 7) != "gene_id") { 
            fprintf(stderr, "ERROR: GTF file formatted incorrectly. Please see the wiki for details.\n"); 
            exit(1); 
        } 
        std::string feature_outer = split_line_outer[2]; 
        std::string gene_outer = split(split(split_line_outer[8], '\"')[1], '\"')[0]; 
 
        if (feature_outer == "CDS") { 
            bool found = (std::find(done.begin(), done.end(), gene_outer) != done.end()); 
            if (found) { 
                continue; 
            } else { 
                done.push_back(gene_outer); 
            } 
            // There may be multiple CDS features per gene. First build codons for the first CDS 
            int first_cds_start = std::stoi(split_line_outer[3]); // expect the GTF is ordered by start position 
            int first_cds_stop = std::stoi(split_line_outer[4]); 
            int codon_counter = 0; // the number of codons we have added so far 
            for (int pos = first_cds_start - 1; pos < first_cds_stop; pos += 3) { 
 
                char nt[3] = { 
                    reference[pos], 
                    reference[pos+1], 
                    reference[pos+2] 
                }; 
 
                // Coordinates are 0-based at this point 
                std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt)); 
                codon_counter += 1; 
 
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
            for (std::string line_inner : gtf_lines) { // find the rest of the CDS features, assuming they are in position order 
                std::vector<std::string> split_line_inner = split(line_inner, '\t'); 
                std::string feature_inner = split_line_inner[2]; 
                std::string gene_inner = split(split(split_line_inner[8], '\"')[1], '\"')[0]; 
                if (feature_inner == "CDS" && gene_outer == gene_inner) { 
                    int inner_cds_start = std::stoi(split_line_inner[3]); 
                    int inner_cds_stop = std::stoi(split_line_inner[4]); 
                    if (inner_cds_start != first_cds_start) { 
                        for (int pos = inner_cds_start - 1; pos < inner_cds_stop; pos += 3) { 
                            char nt[3] = { 
                                reference[pos], 
                                reference[pos+1], 
                                reference[pos+2] 
                            }; 
                            std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt)); 
                            codon_counter += 1; 
 
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
            } 
        } 
    } 
    return codon_map; 
} 
 
 
void translate_main(MAT::Tree *T, std::string output_filename, std::string gtf_filename, std::string fasta_filename ) { 
    std::ifstream fasta_file(fasta_filename); 
    if (!fasta_file) { 
        fprintf(stderr, "ERROR: Could not open the fasta file: %s!\n", fasta_filename.c_str()); 
        exit(1); 
    } 
    std::ifstream gtf_file(gtf_filename); 
    if (!gtf_file) { 
        fprintf(stderr, "ERROR: Could not open the gtf file: %s!\n", gtf_filename.c_str()); 
        exit(1); 
    } 
    std::ofstream output_file(output_filename); 
    if (!output_file) { 
        fprintf(stderr, "ERROR: Could not open file for writing: %s!\n", output_filename.c_str()); 
        exit(1); 
    } 
 
    if (T->condensed_nodes.size() > 0) { 
        T->uncondense_leaves(); 
    } 
 
    std::string reference = build_reference(fasta_file); 
 
    output_file << "node_id\taa_mutations\tnt_mutations\tleaves_sharing_mutations" << '\n'; 
 
    // This maps each position in the reference to a vector of codons. 
    // Some positions may be associated with multiple codons (frame shifts). 
    // The Codons in the map are updated as the tree is traversed 
    std::map<int, std::vector<std::shared_ptr<Codon>>> codon_map = build_codon_map(gtf_file, reference); 
 
    // Traverse the tree in depth-first order. As we descend the tree, mutations at 
    // each node are applied to the respective codon(s) in codon_map. 
    auto dfs = T->depth_first_expansion(); 
    MAT::Node *last_visited = nullptr; 
    for (auto node: dfs) { 
        std::string mutation_result = ""; 
        if(last_visited != node->parent) { 
            // Jumping across a branch, so we need to revert codon mutations up to 
            // the LCA of this node and the last visited node 
            MAT::Node *last_common_ancestor = MAT::LCA(*T, node->identifier, last_visited->identifier); 
            MAT::Node *trace_to_lca = last_visited; 
            while (trace_to_lca != last_common_ancestor) { 
                undo_mutations(trace_to_lca->mutations, codon_map); 
                trace_to_lca = trace_to_lca->parent; 
            } 
        } // If we are visiting a child, we can continue mutating 
        mutation_result = do_mutations(node->mutations, codon_map); 
        if (mutation_result != "") { 
            output_file << node->identifier << '\t' << mutation_result << '\t' << T->get_leaves(node->identifier).size() << '\n'; 
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
        prot_string.resize(prot_string.length() - 1); //remove trailing ';' 
    } 
    if (!nuc_string.empty() && nuc_string.back() == ';') { 
        nuc_string.resize(nuc_string.length() - 1); 
    } 
    if (nuc_string.empty() || prot_string.empty()) { 
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
