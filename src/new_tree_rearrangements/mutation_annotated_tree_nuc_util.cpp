#include "mutation_annotated_tree.hpp"
int8_t Mutation_Annotated_Tree::get_nuc_id(char nuc) {
    int8_t ret = 0b1111;
    switch(nuc) {
        case 'a':
        case 'A': ret = 0b1;
                  break;
        case 'c':
        case 'C': ret = 0b10;
                  break;
        case 'g':
        case 'G': ret = 0b100; 
                  break;
        case 't':
        case 'T': ret = 0b1000; 
                  break;
        case 'R': ret = 0b101;
                  break;
        case 'Y': ret = 0b1010;
                  break;
        case 'S': ret = 0b110;
                  break;
        case 'W': ret = 0b1001;
                  break;
        case 'K': ret = 0b1100;
                  break;
        case 'M': ret = 0b11;
                  break;
        case 'B': ret = 0b1110;
                  break;
        case 'D': ret = 0b1101;
                  break;
        case 'H': ret = 0b1011;
                  break;
        case 'V': ret = 0b111;
        case 'n':
        case 'N': 
        default: ret = 0b1111;
                 break;
    }
    return ret;
}

// Sets bits at positions specified by nuc_vec to 1 in int8
int8_t Mutation_Annotated_Tree::get_nuc_id (std::vector<int8_t> nuc_vec) {
    int8_t ret = 0;
    int8_t one = 1;
    for (auto nuc: nuc_vec) {
        assert((nuc >= 0) && (nuc <=3));
        ret += (one << nuc);
    }
    return ret;
}

// Convert nuc_id back to IUPAC base 
char Mutation_Annotated_Tree::get_nuc (int8_t nuc_id) {
    char ret = 'N';
    //assert ((nuc_id >= 1) && (nuc_id <= 15));
    switch(nuc_id) {
        case 1: ret = 'A';
                break;
        case 2: ret = 'C';
                break;
        case 3: ret = 'M';
                break;
        case 4: ret = 'G';
                break;
        case 5: ret = 'R';
                break;
        case 6: ret = 'S';
                break;
        case 7: ret = 'V';
                break;
        case 8: ret = 'T';
                break;
        case 9: ret = 'W';
                break;
        case 10: ret = 'Y';
                break;
        case 11: ret = 'H';
                break;
        case 12: ret = 'K';
                break;
        case 13: ret = 'D';
                break;
        case 14: ret = 'B';
                break;
        default: ret = 'N';
                 break;
    }
    return ret;
}

// A:0, C:1, G:2, T:3 
int8_t Mutation_Annotated_Tree::get_nt (int8_t nuc_id) {
    int8_t ret = 0;
    switch(nuc_id) {
        case 1: ret = 0;
                break;
        case 2: ret = 1;
                break;
        case 4: ret = 2;
                break;
        case 8: ret = 3;
                break;
        default: ret = -1;
                 break;
    }
    return ret;
}

std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec (char c) {
    switch (c) {
        case 'a':
        case 'A': return std::vector<int8_t>{0};
        case 'c':
        case 'C': return std::vector<int8_t>{1};
        case 'g':
        case 'G': return std::vector<int8_t>{2};
        case 't':
        case 'T': return std::vector<int8_t>{3};
        case 'R': return std::vector<int8_t>{0,2};
        case 'Y': return std::vector<int8_t>{1,3};
        case 'S': return std::vector<int8_t>{1,2};
        case 'W': return std::vector<int8_t>{0,3};
        case 'K': return std::vector<int8_t>{2,3};
        case 'M': return std::vector<int8_t>{0,1};
        case 'B': return std::vector<int8_t>{1,2,3};
        case 'D': return std::vector<int8_t>{0,2,3};
        case 'H': return std::vector<int8_t>{0,1,3};
        case 'V': return std::vector<int8_t>{0,1,2};
        case 'n':
        case 'N': return std::vector<int8_t>{0,1,2,3};
        default: return std::vector<int8_t>{0,1,2,3};
    }
}
