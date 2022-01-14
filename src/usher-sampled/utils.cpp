#include "usher.hpp"
#include <cstdint>
#include <vector>
void convert_mut_type(std::vector<Sampled_Tree_Mutation> &in,
                      std::vector<To_Place_Sample_Mutation> &out) {
    out.reserve(in.size());
    for (const auto &mut : in) {
        if (mut.mut_nuc != 0xf) {
            out.push_back(To_Place_Sample_Mutation(mut.position,mut.chrom_idx,mut.mut_nuc,mut.par_nuc));
        } else {
            if (out.empty() || // first
                out.back().mut_nuc!=0xf|| //last mutation not N
                (out.back().get_end_range() + 1) !=
                    mut.position ||           // not contiguous
                out.back().range == UINT8_MAX // overflow
            ) {
                out.push_back(To_Place_Sample_Mutation(
                    mut.position, mut.chrom_idx, mut.mut_nuc));
            } else {
                out.back().range++;
            }
        }
    }
}