#include "import_vcf.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <string>
#include <unordered_map>
namespace MAT = Mutation_Annotated_Tree;
VCF_Reader::VCF_Reader(boost::iostreams::filtering_istream& instream):instream(instream){
    std::string s;
    while (s.substr(0, 3) != "#CH") {
        std::getline(instream, s);
    }
    MAT::string_split(s, header);
}
Parsed_VCF_Line VCF_Reader::operator()(tbb::flow_control& fc) const {

    // check if reached end-of-file
    int curr_char = instream.peek();
    Parsed_VCF_Line output;
    if (curr_char == EOF){
        fc.stop();
        return output;
    }

    std::string s;
    std::getline(instream, s);
    std::vector<std::string> words;
    MAT::string_split(s, words);

    if (words.size() != header.size()) {
        fprintf(stderr, "ERROR! Incorrect VCF format.\n");
        exit(1);
    }

    std::vector<std::string> alleles;
    alleles.clear();
    output.mutation =
        MAT::Mutation(words[CHROM_IDX], std::stoi(words[POS_IDX]), 0, 0, 0, 1,
                      0, MAT::get_nuc_id(words[REF_IDX][0]));

    std::unordered_map<std::string, nuc_one_hot> *alt_map =
        new std::unordered_map<std::string, nuc_one_hot>;
    //mutated_positions.emplace(output.mutation,alt_map);
    output.mutated = alt_map;
    MAT::string_split(words[ALT_IDX], ',', alleles);
    std::vector<nuc_one_hot> allele_translated;
    allele_translated.reserve(alleles.size());
    for (const auto &allele : alleles) {
        allele_translated.push_back(MAT::get_nuc_id(allele[0]));
    }
    // Ref nuc id uses one-hot encoding (A:0b1, C:0b10, G:0b100,
    // T:0b1000)
    for (size_t j = SAMPLE_START_IDX; j < words.size(); j++) {
        if (isdigit(words[j][0])) {
            int allele_id = std::stoi(words[j]);
            if (allele_id > 0) {
                alt_map->emplace(header[j], allele_translated[allele_id-1]);
            }
        } else {
            alt_map->emplace(header[j], 0xf);
        }
    }
    return output;
}
