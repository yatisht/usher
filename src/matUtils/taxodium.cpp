#include "taxodium.hpp"

#include <google/protobuf/text_format.h>


//#include taxodium pb


// std::vector<std::unordered_map<std::string,std::unordered_map<std::string,std::string>>> convert_meta(std::vector<std::map<std::string,std::map<std::string,std::string>>> &catmeta) {
//     std::vector<std::unordered_map<std::string,std::unordered_map<std::string,std::string>>> out;
//     for (auto &m : catmeta) {
        
//         out.push_back()
//     }
// }

void save_taxodium_tree (MAT::Tree &tree, std::string filename, std::vector<std::map<std::string,std::map<std::string,std::string>>> &catmeta) {

    auto metadata = convert_meta(catmeta); // convert to use unordered_maps for faster lookup


	Taxodium::AllNodeData node_data;
	Taxodium::AllData all_data;

    std::unordered_map<std::string, int32_t> country_map;
    std::unordered_map<std::string, int32_t> lineage_map;
    
    std::vector<std::string> country_list;
    std::vector<std::string> lineage_list; 


    for (auto m : metadata) {
        int ct = 0;
        for (auto &c : m["country"]) {
            if (std::find(country_list.begin(), country_list.end(), c.second) == country_list.end()) {
                // not yet in list
                country_list.push_back(c.second);
                country_map.insert({c.second, ct}); // map this country to a number
                all_data.add_country_mapping(c.second);
                ct += 1;
            }
        }
        ct = 0;
        for (auto &c : m["pangolin_lineage"]) {
            if (std::find(lineage_list.begin(), lineage_list.end(), c.second) == lineage_list.end()) {
                lineage_list.push_back(c.second);
                lineage_map.insert({c.second, ct}); 
                all_data.add_lineage_mapping(c.second);
                ct += 1;
            }
        }
    }
    

	int count = 0;
	TIMEIT();

    auto dfs = tree.depth_first_expansion();

    for (size_t idx = 0; idx < dfs.size(); idx++) {
		if (count > 20) {
			break;
		}
		count++;
		MAT::Node *node = dfs[idx];
		
		node_data.add_names(node->identifier);


        bool found_genbank = false;
        for (auto m : metadata) {
            if (m["genbank_accession"].find(node->identifier) != m["genbank_accession"].end()) {
                node_data.add_genbanks(m["genbank_accession"][node->identifier]);
                found_genbank = true;    
            }
            if (m["country"].find(node->identifier) != m["country"].end()) {
                node_data.add_countries(country_map[m["country"][node->identifier]]);
            }
        }
        if (!found_genbank) {
            fprintf(stderr, ("ERROR: No genbank field in metadata for " + node->identifier + "\n").c_str());
        }

//        node_data.add_countries;
		
        // auto meta = data.add_metadata();
        // for (size_t k = 0; k < node->clade_annotations.size(); k++) {
        //     meta->add_clade_annotations(node->clade_annotations[k]);
        // }
        // auto mutation_list = data.add_node_mutations();
        // for (auto m: node->mutations) {
        //     auto mut = mutation_list->add_mutation();
        //     mut->set_chromosome(m.chrom);
        //     mut->set_position(m.position);

        //     if (m.is_masked()) {
        //         mut->set_ref_nuc(-1);
        //         mut->set_par_nuc(-1);
        //     } else {
        //         int8_t j = get_nt(m.ref_nuc);
        //         assert (j >= 0);
        //         mut->set_ref_nuc(j);

        //         j = get_nt(m.par_nuc);
        //         assert(j >= 0);
        //         mut->set_par_nuc(j);

        //         mut->clear_mut_nuc();
        //         for (auto nuc: get_nuc_vec_from_id(m.mut_nuc)) {
        //             mut->add_mut_nuc(nuc);
        //         }
        //     }
        // }
    }

    // // Add condensed nodes
    // for (auto cn: tree.condensed_nodes) {
    //     auto cn_ptr = data.add_condensed_nodes();
    //     cn_ptr->set_node_name(cn.first);
    //     for (auto lid: cn.second) {
    //         cn_ptr->add_condensed_leaves(lid);
    //     }
    // }

    // Boost library used to stream the contents to the output protobuf file in
    // uncompressed or compressed .gz format
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outbuf;

    if (filename.find(".gz\0") != std::string::npos) {
        try {
            outbuf.push(boost::iostreams::gzip_compressor());
            outbuf.push(outfile);
            std::ostream outstream(&outbuf);
            node_data.SerializeToOstream(&outstream);
            std::string s;
            google::protobuf::TextFormat::PrintToString(node_data, &s);
            std::cout << s << '\n';
            boost::iostreams::close(outbuf);
            outfile.close();
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        node_data.SerializeToOstream(&outfile);
        std::string s;
        google::protobuf::TextFormat::PrintToString(node_data, &s);
        std::cout << s << '\n';
        outfile.close();
    }
}