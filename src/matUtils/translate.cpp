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

char complement(char nt) {
    auto it = complement_map.find(nt);
    if (it == complement_map.end()) {
        return 'N'; // ambiguous, couldn't resolve nt
    } else {
        return it->second;
    }
}

// Maps a genomic coordinate to a list of codons it is part of
std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> build_codon_map(std::ifstream &gtf_file, std::string reference) {
    std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> codon_map;
    std::string gtf_line;
    std::vector<std::string> gtf_lines;
    std::vector<std::string> done;

    while (std::getline(gtf_file, gtf_line)) {
        gtf_lines.push_back(gtf_line);
    }
    int curr_line = -1;
    for (std::string line_outer : gtf_lines) {
        curr_line += 1;

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
        char strand_outer = split_line_outer[6][0];

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
            if (strand_outer == '+') {
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
            } else {
                for (int pos = first_cds_stop - 1; pos > first_cds_start; pos -= 3) {

                    char nt[3] = {
                        complement(reference[pos]),
                        complement(reference[pos-1]),
                        complement(reference[pos-2])
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

                    it = codon_map.find(pos-1);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos-1, {c}});
                    } else {
                        (it->second).push_back(c);
                    }

                    it = codon_map.find(pos-2);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos-2, {c}});
                    } else {
                        (it->second).push_back(c);
                    }
                }
            }
            for (std::string line_inner : gtf_lines) { // find the rest of the CDS features, assuming they are in position order

                if (line_inner[0] == '#' || line_inner[0] == '\n') {
                    continue;
                }
                std::vector<std::string> split_line_inner = split(line_inner, '\t');
                std::string feature_inner = split_line_inner[2];
                std::string gene_inner = split(split(split_line_inner[8], '\"')[1], '\"')[0];
                if (feature_inner == "CDS" && gene_outer == gene_inner) {
                    int inner_cds_start = std::stoi(split_line_inner[3]);
                    int inner_cds_stop = std::stoi(split_line_inner[4]);
                    char strand_inner = split_line_inner[6][0];
                    if (strand_inner == '+') {
                        if (inner_cds_start != first_cds_start || strand_outer != strand_inner) {
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
                    } else {
                        if (inner_cds_start != first_cds_start || strand_outer != strand_inner) {
                            for (int pos = inner_cds_stop - 1; pos > inner_cds_start; pos -= 3) {
                                char nt[3] = {
                                    complement(reference[pos]),
                                    complement(reference[pos-1]),
                                    complement(reference[pos-2])
                                };
                                std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt));
                                codon_counter += 1;

                                auto it = codon_map.find(pos);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos-1);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos-1, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos-2);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos-2, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return codon_map;
}

void translate_main(MAT::Tree *T, std::string output_filename, std::string gtf_filename, std::string fasta_filename) {
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

    output_file << "node_id\taa_mutations\tnt_mutations\tcodon_changes\tleaves_sharing_mutations" << '\n';

    // This maps each position in the reference to a vector of codons.
    // Some positions may be associated with multiple codons (frame shifts).
    // The Codons in the map are updated as the tree is traversed
    std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> codon_map = build_codon_map(gtf_file, reference);

    // Traverse the tree in depth-first order. As we descend the tree, mutations at
    // each node are applied to the respective codon(s) in codon_map.
    auto dfs = T->depth_first_expansion();
    MAT::Node *last_visited = nullptr;
    for (auto &node: dfs) {
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
        mutation_result = do_mutations(node->mutations, codon_map, false);
        if (mutation_result != "") {
            output_file << node->identifier << '\t' << mutation_result << '\t' << T->get_leaves(node->identifier).size() << '\n';
        }
        last_visited = node;
    }
}

// This is used for taxodium output. It translates each node and saves metadata to node_data along the way
void translate_and_populate_node_data(MAT::Tree *T, std::string gtf_filename, std::string fasta_filename, Taxodium::AllNodeData *node_data, Taxodium::AllData *all_data, std::unordered_map<std::string, std::vector<std::string>> &metadata, MetaColumns fixed_columns, std::vector<GenericMetadata> &generic_metadata, float x_scale, bool include_nt) {
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

    if (T->condensed_nodes.size() > 0) {
        T->uncondense_leaves();
    }

    T->rotate_for_display();
    std::string reference = build_reference(fasta_file);
    std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> codon_map = build_codon_map(gtf_file, reference);
    auto dfs = T->depth_first_expansion();

    MAT::Node *last_visited = nullptr;

    std::unordered_map<std::string, int32_t> seen_mutations_map;
    std::unordered_map<std::string, int32_t> index_map; // map node id to index in protobuf arrays
    std::unordered_map<std::string, float> y_map;
    std::unordered_map<std::string, float> branch_length_map;
    std::map<size_t, std::vector<MAT::Node *>, std::greater<int>> by_level;
    int32_t count = 0;
    float curr_x_value = 0;

    all_data->add_mutation_mapping(""); // no mutations
    std::vector<MAT::Node *> leaves;
    int32_t mutation_counter = 0;

    // DFS to translate aa mutations, adding to Taxodium pb objects along the way
    for (auto &node: dfs) {
        if (node->is_leaf()) {
            leaves.push_back(node);
        }
        by_level[node->level].push_back(node); // store nodes by level for later step
        index_map[node->identifier] = count;

        // If we are jumping across a branch relative to the last visited node, reset mutations
        // and x-value to LCA of last node and this node
        if(last_visited != node->parent) {
            MAT::Node *last_common_ancestor = MAT::LCA(*T, node->identifier, last_visited->identifier);
            MAT::Node *trace_to_lca = last_visited;
            while (trace_to_lca != last_common_ancestor) {
                undo_mutations(trace_to_lca->mutations, codon_map);
                trace_to_lca = trace_to_lca->parent;
            }
            curr_x_value = branch_length_map[trace_to_lca->identifier] + node->mutations.size();
        } else {
            curr_x_value += node->mutations.size();
        }
        branch_length_map[node->identifier] = curr_x_value;

        // First collect (syn and nonsyn) nucleotide mutations with a fake gene called nt
        std::string mutation_result = "";
        if (include_nt) {
            for (auto m : node->mutations) {
                mutation_result += "nt:";
                mutation_result += MAT::get_nuc(m.par_nuc);
                mutation_result +=  "_" + std::to_string(m.position) + "_";
                mutation_result += MAT::get_nuc(m.mut_nuc);
                mutation_result +=  ";";
            }
        }
        // Do mutations
        Taxodium::MutationList *mutation_list = node_data->add_mutations();

        // This string is a semicolon separated list of mutations in format
        // [orf]:[orig aa]_[orf num]_[new aa]
        // e.g. S:K_200_V;ORF1a:G_240_N
        mutation_result += do_mutations(node->mutations, codon_map, true);


        if (node->is_root()) {
            // For the root node, modify mutation_result with "fake" mutations,
            // to enable correct coloring by amino acid in Taxodium
            std::unordered_map<std::string, bool> done_codons = {}; // some codons are duplicated in codon_map, track them
            std::string root_mutations = ""; // add "mutations" at the root
            for (int32_t pos = 0; pos < (int32_t) reference.length(); pos++) {
                if (codon_map.find(pos) == codon_map.end()) {
                    continue;
                }
                for (auto codon_ptr : codon_map[pos]) {
                    std::string codon_id = codon_ptr->orf_name + ":" + std::to_string(codon_ptr->codon_number+1);
                    if (done_codons.find(codon_id) != done_codons.end()) {
                        continue;
                    }
                    done_codons[codon_id] = true;
                    root_mutations += codon_ptr->orf_name + ":X_" + std::to_string(codon_ptr->codon_number+1) + "_" + codon_ptr->protein + ";";
                }
            }
            mutation_result = root_mutations;
        }
        if (mutation_result != "") {
            // Add mutations to protobuf object
            for (auto m : split(mutation_result, ';')) {
                if (seen_mutations_map.find(m) == seen_mutations_map.end()) {
                    mutation_counter++;
                    seen_mutations_map[m] = mutation_counter;
                    all_data->add_mutation_mapping(m);
                    mutation_list->add_mutation(mutation_counter);
                } else {
                    mutation_list->add_mutation(seen_mutations_map[m]);
                }
            }
        }

        node_data->add_x(branch_length_map[node->identifier] * x_scale);
        node_data->add_y(0); // temp value, set later
        node_data->add_epi_isl_numbers(0); // not currently set
        node_data->add_num_tips(T->get_leaves(node->identifier).size());

        if (node->identifier.substr(0,5) == "node_") {
            //internal nodes don't have metadata, so populate with empty data
            node_data->add_names("");
            if (fixed_columns.date_column > -1) {
                node_data->add_dates(0);
            }
            if (fixed_columns.genbank_column > -1) {
                node_data->add_genbanks("");
            }
            for (const auto &m : generic_metadata) {
                m.protobuf_data_ptr->add_node_values(0); // no metadata for this node
            }
        } else if (metadata.find(node->identifier) == metadata.end()) {
            node_data->add_names(split(node->identifier, '|')[0]);

            if (fixed_columns.date_column > -1) {
                node_data->add_dates(0);
            }
            if (fixed_columns.genbank_column > -1) {
                node_data->add_genbanks("");
            }
            for (const auto &m : generic_metadata) {
                m.protobuf_data_ptr->add_node_values(0); // no metadata for this node
            }

        } else {

            // All of the metadata values (integer-encoded for those with mappings) for the current node
            std::vector<std::string> meta_fields = metadata[node->identifier];

            if (fixed_columns.date_column > -1) {
                int32_t date = std::stoi(meta_fields[fixed_columns.date_column]);
                node_data->add_dates(date);
            }
            if (fixed_columns.genbank_column > -1) {
                node_data->add_genbanks(meta_fields[fixed_columns.genbank_column]);
            }

            node_data->add_names(split(node->identifier, '|')[0]);

            for (const auto &m : generic_metadata) {
                m.protobuf_data_ptr->add_node_values(std::stoi(meta_fields[m.column])); // lookup the encoding for this value
            }

        }

        if (node->parent == nullptr) {
            node_data->add_parents(0); // root node
        } else {
            node_data->add_parents(index_map[node->parent->identifier]);
        }

        last_visited = node;
        count++;
    }

    // Set a y-value for each leaf
    int32_t i = 1;
    std::reverse(leaves.begin(), leaves.end());
    for (auto &leaf : leaves) {
        if (leaf->identifier == "CHN/YN-0306-466/2020|MT396241.1|2020-03-06") {
            node_data->set_y(index_map[leaf->identifier], 0.0);
            continue;
        }
        node_data->set_y(index_map[leaf->identifier], (float) i / 40000);
        i++;
    }

    // Travel by level up the tree assigning y-values to each internal node
    for (auto &node_list : by_level) { // in descending order by level
        for (auto &bylevel_node : node_list.second) {
            float children_mean_y = 0;
            int num_children = (bylevel_node->children).size();
            for (auto &child : bylevel_node->children) {
                children_mean_y += node_data->y(index_map[child->identifier]);
            }
            if (num_children > 0) {
                children_mean_y /= num_children;
            } else { // leaf
                continue;
            }
            // y-value of a node is the average y-value of its children
            node_data->set_y(index_map[bylevel_node->identifier], children_mean_y);
        }
    }
}
std::string do_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map, bool taxodium_format) {
    std::string prot_string = "";
    std::string nuc_string = "";
    std::string cchange_string = "";
    std::sort(mutations.begin(), mutations.end());
    std::unordered_map<std::string, std::set<MAT::Mutation>> codon_to_nt_map;
    std::unordered_map<std::string, std::string> codon_to_changestring_map;
    std::unordered_map<std::string, char> orig_proteins;
    std::vector<std::shared_ptr<Codon>> affected_codons;

    for (auto &m : mutations) {
        char mutated_nuc = MAT::get_nuc(m.mut_nuc);
        char par_nuc = MAT::get_nuc(m.par_nuc);
        int pos = m.position - 1;
        auto codon_map_it = codon_map.find(pos);
        if (codon_map_it == codon_map.end()) {
            continue; // Not a coding mutation
        } else {
            // Mutate each codon associated with this position
            for (auto codon_ptr : codon_map_it->second) {
                std::string codon_id = codon_ptr->orf_name + ':' + std::to_string(codon_ptr->codon_number+1);
                //first, update the codon to match the parent state instead of the reference state as part of the codon output
                codon_ptr->mutate(pos, par_nuc);
                auto orig_it = orig_proteins.find(codon_id);
                if (orig_it == orig_proteins.end()) {
                    orig_proteins.insert({codon_id, codon_ptr->protein});
                }
                if (std::find(affected_codons.begin(), affected_codons.end(), codon_ptr) == affected_codons.end()) {
                    affected_codons.push_back(codon_ptr);
                }
                std::string original_codon = codon_ptr->nucleotides;
                //then update it again to match the mutated state
                codon_ptr->mutate(pos, mutated_nuc);
                // store a string representing the original and new codons in nucleotides
                // this may incorporate multiple nucleotide mutations, just accounting for the original and end states.
                std::string changestring = original_codon + ">" + codon_ptr->nucleotides;
                codon_to_changestring_map.insert({codon_id, changestring});
                // Build a map of codons and their associated nt mutations
                auto to_nt_it = codon_to_nt_map.find(codon_id);
                if (to_nt_it == codon_to_nt_map.end()) {
                    codon_to_nt_map.insert({codon_id, {m}});
                } else {
                    to_nt_it->second.insert(m);
                }
            }
        }
    }

    for (auto codon_ptr : affected_codons) {
        std::string codon_id = codon_ptr->orf_name + ':' + std::to_string(codon_ptr->codon_number+1);
        char orig_protein = orig_proteins.find(codon_id)->second;
        if (taxodium_format) {
            if (orig_protein == codon_ptr->protein) { // exclude synonymous mutations
                continue;
            }
            prot_string += split(codon_id, ':')[0] + ':' + orig_protein + '_' + split(codon_id, ':')[1] + '_' + codon_ptr->protein + ';';
        } else {
            prot_string += split(codon_id, ':')[0] + ':' + orig_protein + split(codon_id, ':')[1] + codon_ptr->protein + ';';
        }
        auto codon_it = codon_to_nt_map.find(codon_id);
        for (auto &m : codon_it->second) {
            nuc_string += m.get_string() + ",";
        }

        if (!nuc_string.empty() && nuc_string.back() == ',') {
            nuc_string.resize(nuc_string.length() - 1); // remove trailing ','
            nuc_string += ';';
        }
        std::string changestring = codon_to_changestring_map.find(codon_id)->second;
        cchange_string += changestring + ";";
    }

    if (!nuc_string.empty() && nuc_string.back() == ';') {
        nuc_string.resize(nuc_string.length() - 1); // remove trailing ';'
    }
    if (!prot_string.empty() && prot_string.back() == ';') {
        prot_string.resize(prot_string.length() - 1); //remove trailing ';'
    }
    if (!cchange_string.empty() && cchange_string.back() == ';') {
        cchange_string.resize(cchange_string.length() - 1); //remove trailing ';'
    }
    if (nuc_string.empty() || prot_string.empty() || cchange_string.empty()) {
        return "";
    } else if(taxodium_format) { // format string for taxodium pb
        return prot_string;
    } else { // format string for TSV output
        return prot_string + '\t' + nuc_string + '\t' + cchange_string;
    }
}

void undo_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map) {
    for (auto &m: mutations) {
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



/// Taxodium protobuf output functions below


// Helper function to format one attribute into taxodium encoding for a SingleValuePerNode metadata type
void populate_generic_metadata(int attribute_column, std::vector<std::string> &attributes, std::unordered_map<std::string, std::string> &seen_map, int &encoding_counter, Taxodium::MetadataSingleValuePerNode *single) {
    if (attributes[attribute_column] != "") {
        std::string attr_val = attributes[attribute_column];
        if (seen_map.find(attr_val) == seen_map.end()) {
            encoding_counter++;
            std::string encoding_str = std::to_string(encoding_counter);
            seen_map[attr_val] = encoding_str;
            single->add_mapping(attr_val);
            attributes[attribute_column] = encoding_str;
        } else {
            attributes[attribute_column] = seen_map[attr_val];
        }
    } else {
        attributes[attribute_column] = "0";
    }
}

// Helper function to populate non-generic metadata types that have mapping encodings.
void populate_fixed_metadata(std::string name, int attribute_column, std::vector<std::string> &attributes, std::unordered_map<std::string, std::string> &seen_map, int &encoding_counter, Taxodium::AllData &all_data) {
    if (attributes[attribute_column] != "") {
        if (seen_map.find(attributes[attribute_column]) == seen_map.end()) {
            encoding_counter++;
            std::string encoding_str = std::to_string(encoding_counter);
            seen_map[attributes[attribute_column]] = encoding_str;
            if (name == "date") { // only date for now
                all_data.add_date_mapping(attributes[attribute_column]);
            }
            attributes[attribute_column] = encoding_str;
        } else {
            attributes[attribute_column] = seen_map[attributes[attribute_column]];
        }
    } else {
        attributes[attribute_column] = "0";
    }
}

// Store the metadata differently for taxodium pb format
std::unordered_map<std::string, std::vector<std::string>> read_metafiles_tax(std::vector<std::string> filenames, Taxodium::AllData &all_data, Taxodium::AllNodeData *node_data, MetaColumns &columns, std::vector<GenericMetadata> &generic_metadata, std::vector<std::string> additional_meta_fields) {

    /* We look for the following metadata fields:
     * strain, genbank_accession, country, date, pangolin_lineage
     * strain is the sample ID and is required. The rest may be missing.
     * Additional fields to look for are specified with -F
     */


    int32_t date_ct = 0;
    all_data.add_date_mapping("");
    std::unordered_map<std::string, std::string> seen_dates_map;

    int num_generic = 0;
    bool country_found = false;
    bool lineage_found = false;

    std::unordered_map<std::string, std::vector<std::string>> metadata;

    // First parse all files into metadata map
    std::vector<std::string> header;
    int additional_fields = 0; // Number of new fields in each metadata file
    for (std::string f : filenames) {
        std::ifstream infile(f);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the file: %s!\n", f.c_str());
            exit(1);
        }
        std::string line;
        char delim = '\t';
        if (f.find(".csv\0") != std::string::npos) {
            delim = ',';
        }
        bool first = true;
        int strain_column = -1;
        std::unordered_map<std::string, bool> seen_in_this_file; // ignore duplicates per file
        while (std::getline(infile, line)) {
            if (line.length() < 3) {
                continue;
            }
            std::vector<std::string> words;
            if (line[line.size()-1] == '\r') {
                line = line.substr(0, line.size()-1);
            }
            MAT::string_split(line, delim, words);
            if (first) { // header line
                int field_count = 0;
                for (int i = 0; i < (int) words.size(); i++) { // for each column name
                    header.push_back(words[i]);
                    field_count++;
                    if (words[i] == "strain") {
                        strain_column = i;
                    }
                }
                additional_fields = field_count;
                first = false;
                if (strain_column == -1) {
                    fprintf(stderr, "The column \"strain\" (sample ID) is missing from at least one metadata file.\n");
                    exit(1);
                }
                continue;
            }
            std::string key = words[strain_column];
            if ( seen_in_this_file.find(key) != seen_in_this_file.end() ) {
                continue; // ignore duplicates in each metadata file
            }
            seen_in_this_file[key] = true;

            int prev_header_size = header.size() - additional_fields;

            if (metadata.find(key) == metadata.end()) {
                // if we haven't seen this sample yet
                metadata[key] = std::vector<std::string>();
                while (metadata[key].size() < prev_header_size) {
                    metadata[key].push_back("");
                }
            }

            for (int i = 0; i < words.size(); i++) {
                // Check all metadata fields up to this one
                // If the same field exists earlier, copy non-empty
                // values into the first column of the field
                std::string word = words[i];
                metadata[key].push_back(word);
                for (int j = 0; j < prev_header_size; j++) {
                    if (header[j] == header[prev_header_size + i]) {
                        if (words[i] != "") {
                            metadata[key][j] = words[i];
                        }
                        break;
                    }
                }
            }

            // fills out empty columns
            while(metadata[key].size() < header.size()) {
                metadata[key].push_back("");
            }

        }
        infile.close();
    }
    for(auto &v : metadata) {
        // fill out empty columns
        while(metadata[v.first].size() < header.size()) {
            metadata[v.first].push_back(""); // handles empty metadata in the last columns
        }
    }
    // Then use map to make taxodium encodings and check for defined/generic fields
    // If the same column is present in multiple metadata files (or multiple times in a file),
    // the non-empty values are condensed into the first column of that name.
    std::unordered_map<std::string, bool> done_fields;
    for (int i = 0; i < (int) header.size(); i++) {
        std::string field = header[i];
        if (done_fields.find(field) != done_fields.end()) {
            continue; // already included this field
        }
        done_fields[field] = true;
        if (field == "strain") {
            columns.strain_column = i;
        } else if (field == "genbank_accession") {
            columns.genbank_column = i;
        } else if (field == "date") {
            columns.date_column = i;
        }

        // Encode everything else as generic
        if (field != "strain" && field != "genbank_accession" && field != "date") {
            if (std::find(additional_meta_fields.begin(), additional_meta_fields.end(), field) != additional_meta_fields.end() || field == "pangolin_lineage" || field == "country") {
                Taxodium::MetadataSingleValuePerNode *metadata_single = node_data->add_metadata_singles();

                // Lineage and country are expected to be named "Lineage" and "Country" in Taxodium, so
                // rename the standard column names. Other custom metadata fields will use their column
                // name as the pb field name
                if (field == "pangolin_lineage") {
                    lineage_found = true;
                    metadata_single->set_metadata_name("Lineage");
                } else if (field == "country") {
                    country_found = true;
                    metadata_single->set_metadata_name("Country");
                } else {
                    metadata_single->set_metadata_name(field);
                }
                metadata_single->add_mapping("");
                std::unordered_map<std::string, std::string> empty_map;
                GenericMetadata new_meta = {
                    .name = field,
                    .column = i,
                    .index = num_generic,
                    .count = 0,
                    .seen = empty_map,
                    .protobuf_data_ptr = metadata_single
                };
                generic_metadata.push_back(new_meta);
                num_generic ++;
            }
        }
    }
    char found[]  = "(FOUND)";
    char missing[] = "(MISSING)";
    fprintf(stderr, "\nLooking for the following metadata fields:\n");
    fprintf(stderr, "-- %s %s\n", found, "strain (this is the sample ID)");
    fprintf(stderr, "-- %s %s\n", columns.genbank_column > -1 ? found : missing, "genbank_accession");
    fprintf(stderr, "-- %s %s\n", columns.date_column > -1 ? found : missing, "date");
    fprintf(stderr, "-- %s %s\n", country_found ? found : missing, "country");
    fprintf(stderr, "-- %s %s\n", lineage_found ? found : missing, "pangolin_lineage");
    fprintf(stderr, "\nIf any of the above fields are missing, importing into Taxodium may not work properly.\n");

    if (additional_meta_fields.size() > 0) {
        additional_meta_fields.erase(std::remove(additional_meta_fields.begin(), additional_meta_fields.end(), "strain"), additional_meta_fields.end());
        fprintf(stderr, "\nThe following additional metadata fields were specified:\n");
        for (std::string name : additional_meta_fields) {
            fprintf(stderr, "-- %s\n", name.c_str());
        }
    }


    for (auto &v : metadata) { // for each metadata value
        // Build mappings for generic metadata types
        for (int i = 0; i < (int) generic_metadata.size(); i++) {
            populate_generic_metadata(generic_metadata[i].column, v.second, generic_metadata[i].seen, generic_metadata[i].count, node_data->mutable_metadata_singles(i));
        }

        // Build mappings for fixed metadata types
        if (columns.date_column > -1) {
            populate_fixed_metadata("date", columns.date_column, v.second, seen_dates_map, date_ct, all_data);
        }

    }


    fprintf(stderr, "\nPerforming conversion.\n");

    return metadata;
}
void save_taxodium_tree (MAT::Tree &tree, std::string out_filename, std::vector<std::string> meta_filenames, std::string gtf_filename, std::string fasta_filename, std::string title, std::string description, std::vector<std::string> additional_meta_fields, float x_scale, bool include_nt) {

    // These are the taxodium pb objects
    Taxodium::AllNodeData *node_data = new Taxodium::AllNodeData();
    Taxodium::AllData all_data;

    all_data.set_tree_title(title);
    all_data.set_tree_description(description);

    // For fixed fields
    MetaColumns columns = {
        .strain_column = -1,
        .date_column = -1,
        .genbank_column = -1,
    };

    std::vector<GenericMetadata> generic_metadata;

    std::unordered_map<std::string, std::vector<std::string>> metadata = read_metafiles_tax(meta_filenames, all_data, node_data, columns, generic_metadata, additional_meta_fields);
    TIMEIT();

    // Fill in the taxodium data while doing aa translations
    translate_and_populate_node_data(&tree, gtf_filename, fasta_filename, node_data, &all_data, metadata, columns, generic_metadata, x_scale, include_nt);
    all_data.set_allocated_node_data(node_data);

    // Boost library used to stream the contents to the output protobuf file in
    // uncompressed or compressed .gz format
    std::ofstream outfile(out_filename, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outbuf;

    if (out_filename.find(".gz\0") != std::string::npos) {
        try {
            outbuf.push(boost::iostreams::gzip_compressor());
            outbuf.push(outfile);
            std::ostream outstream(&outbuf);
            all_data.SerializeToOstream(&outstream);
            boost::iostreams::close(outbuf);
            outfile.close();
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        all_data.SerializeToOstream(&outfile);
        outfile.close();
    }
}
// end taxodium output functions
