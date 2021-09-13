#include "transpose_vcf.hpp"
#include <cstdio>
struct Sample_Pos_Mut_Wrap {
    void add_Not_N(int position, uint8_t allele) {
    }
    void add_N(int first, int second) {
    }
};

struct All_Sample_Appender {
    Sample_Pos_Mut_Wrap set_name(std::string &&name) {
        puts(name.c_str());
        return Sample_Pos_Mut_Wrap{};
    }
};
int main(int argc, char **argv) {
    All_Sample_Appender appender;
    load_mutations(argv[1], 80, appender);
}
