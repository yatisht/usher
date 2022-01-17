#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "usher.hpp"
#include <cstdint>
#include <vector>
void convert_mut_type(std::vector<MAT::Mutation> &in,
                      std::vector<To_Place_Sample_Mutation> &out) {
    out.reserve(in.size());
    for (const auto &mut : in) {
        if (mut.get_mut_one_hot() != 0xf) {
            out.push_back(To_Place_Sample_Mutation(mut.get_position(),mut.get_chromIdx(),mut.get_mut_one_hot(),mut.get_par_one_hot()));
        } else {
            if (out.empty() || // first
                out.back().mut_nuc!=0xf|| //last mutation not N
                (out.back().get_end_range() + 1) !=
                    mut.get_position() ||           // not contiguous
                out.back().range == UINT8_MAX // overflow
            ) {
                out.push_back(To_Place_Sample_Mutation(
                    mut.get_position(), mut.get_chromIdx(), mut.get_mut_one_hot()));
            } else {
                out.back().range++;
            }
        }
    }
}