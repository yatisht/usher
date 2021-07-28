#include "common.hpp"

struct Codon {
    int startPosition;
    char nt[3];
    char prot;

    Codon (std::string orf, int pos, char nt1, char nt2, char nt3) {
        startPosition = pos;
        nt[0] = nt1;
        nt[1] = nt2;
        nt[2] = nt3;
        prot = 'K';
    }

    inline std::string get_string() const {
            return std::to_string(startPosition) + ':' + nt[0] + nt[1] + nt[2];
        }
};


void translate_main(po::parsed_options parsed);