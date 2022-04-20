#include "common.hpp"
#include "taxodium.pb.h"

// column numbers of fixed metadata
typedef struct {
    int strain_column;
    int date_column;
    int genbank_column;
} MetaColumns;

typedef struct {
    std::string name;
    int column; // column in metadata file
    int index; // index in proto list
    int32_t count; // used for building mapping
    std::unordered_map<std::string, std::string> seen; // used for building mapping
    Taxodium::MetadataSingleValuePerNode *protobuf_data_ptr;
} GenericMetadata;

static std::unordered_map<std::string, char> translation_map = {
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'}, {"GCN", 'A'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGY", 'C'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAY", 'D'},
    {"GAA", 'E'}, {"GAG", 'E'}, {"GAR", 'E'},
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTY", 'F'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}, {"GGN", 'G'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAY", 'H'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATH", 'I'},
    {"AAA", 'K'}, {"AAG", 'K'}, {"AAR", 'K'},
    {"TTA", 'L'}, {"TTG", 'L'}, {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'}, {"YTR", 'L'}, {"CTN", 'L'},
    {"ATG", 'M'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAY", 'N'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'}, {"CCN", 'P'},
    {"CAA", 'Q'}, {"CAG", 'Q'}, {"CAR", 'Q'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'}, {"AGA", 'R'}, {"AGG", 'R'}, {"CGN", 'R'}, {"MGR", 'R'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'}, {"AGT", 'S'}, {"AGC", 'S'}, {"TCN", 'S'}, {"AGY", 'S'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'}, {"ACN", 'T'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'}, {"GTN", 'V'},
    {"TGG", 'W'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAY", 'Y'},
    {"TAG", '*'}, {"TAA", '*'}, {"TGA", '*'}
};

static std::unordered_map<char, char> complement_map = {
    {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'},
    {'M', 'K'}, {'R', 'Y'}, {'W', 'W'}, {'S', 'S'},
    {'Y', 'R'}, {'K', 'M'}, {'V', 'B'}, {'H', 'D'},
    {'D', 'H'}, {'B', 'V'}, {'N', 'N'}
};

struct Codon {
    std::string orf_name;
    std::string nucleotides;
    int codon_number;
    int start_position;
    char protein;

    // Translate codon to amino acid, allowing for ambiguous codons
    inline char translate_codon(std::string nt) {
        auto it = translation_map.find(nt);
        if (it == translation_map.end()) {
            return 'X'; // ambiguous, couldn't resolve aa
        } else {
            return it->second;
        }
    }

    inline void mutate(int nuc_pos, char mutated_nuc) {
        // The nt to mutate is the difference between the
        // genomic coordinate of the mutated nt and the
        // starting coordinate of the codon
        nucleotides[abs(nuc_pos-start_position)] = mutated_nuc;
        protein = translate_codon(nucleotides);
    }

    Codon (std::string _orf_name, int _codon_number, int _start_position, char nt[3]) {
        orf_name = _orf_name;
        start_position = _start_position;
        codon_number = _codon_number;
        nucleotides = "";
        nucleotides += nt[0];
        nucleotides += nt[1];
        nucleotides += nt[2];
        protein = translate_codon(nt);
    }

    inline std::string get_string() const {
        return std::to_string(start_position) + ':'
               + nucleotides[0]
               + nucleotides[1]
               + nucleotides[2]
               + '=' + protein;
    }
};

std::string do_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map, bool taxodium_format);
void translate_main(MAT::Tree *T, std::string output_filename, std::string gff_filename, std::string fasta_filename);
void translate_and_populate_node_data(MAT::Tree *T, std::string gtf_filename, std::string fasta_filename, Taxodium::AllNodeData *node_data, Taxodium::AllData *all_data, std::unordered_map<std::string, std::vector<std::string>> &metadata, MetaColumns fixed_columns, std::vector<GenericMetadata> &generic_metadata, float x_scale, bool include_nt);
void cleanup_codon_map(std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map);
void undo_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map);
char complement(char nt);
void save_taxodium_tree (MAT::Tree &tree, std::string out_filename, std::vector<std::string> meta_filenames, std::string gtf_filename, std::string fasta_filename, std::string title, std::string description, std::vector<std::string> additional_meta_fields, float x_scale, bool include_nt);
std::unordered_map<std::string, std::vector<std::string>> read_metafiles_tax(std::vector<std::string> filenames, Taxodium::AllData &all_data, Taxodium::AllNodeData *node_data, MetaColumns &columns, std::vector<GenericMetadata> &generic_metadata, std::vector<std::string> additional_meta_fields);
