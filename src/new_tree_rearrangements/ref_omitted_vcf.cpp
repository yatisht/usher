#include "src/new_tree_rearrangements/check_samples.hpp"
#include "vcf.pb.h"
#include "tree_rearrangement_internal.hpp"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>
#include <fstream>
void save_variant_samples(const Original_State_t& origin_state,const std::string& out_name){
    std::unordered_map<std::string,uint32_t> samples;
    VCF::data data;
    uint32_t sample_idx=0;
    for(auto sample:origin_state){
        data.add_samples(sample.first);
        samples.emplace(sample.first,sample_idx++);
    }
    for(const auto& pos: mutated_positions){
        data.add_pos(pos.first.get_position());
        auto genotype=data.add_genotypes();
        for(const auto& sample_mut:*pos.second){
            genotype->add_sample_idx(samples[sample_mut.first]);
            genotype->add_allele(sample_mut.second);
        }
    }
    
}