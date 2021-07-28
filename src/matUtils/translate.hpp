#include "common.hpp"



struct Codon {
    int startPosition;
    char nucleotides[3];
    char protein;

    Codon (std::string orf, int pos, std::string nt, char prot) {
        startPosition = pos;
        nucleotides[0] = nt[0];
        nucleotides[1] = nt[1];
        nucleotides[2] = nt[2];
        protein = prot;
    }

    inline std::string get_string() const {
            return std::to_string(startPosition) + ':'
                + nucleotides[0]
                + nucleotides[1]
                + nucleotides[2] 
                + '=' + protein;
        }
};

char translate_codon(std::string nt);
void translate_main(po::parsed_options parsed);