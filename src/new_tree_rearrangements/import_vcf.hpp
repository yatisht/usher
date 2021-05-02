#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <tbb/pipeline.h>
#include <vector>
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <utility>
struct Parsed_VCF_Line{
    Mutation_Annotated_Tree::Mutation mutation;
    //std::unordered_map<std::string, nuc_one_hot> *mutated;
    std::unordered_map<std::string, nuc_one_hot> mutated;
};
struct VCF_Reader{
    static const int CHROM_IDX=0;
    static const int POS_IDX=1;
    static const int REF_IDX=3;
    static const int ALT_IDX=4;
    static const int SAMPLE_START_IDX=9;
    boost::iostreams::filtering_istream& instream;
    std::vector<std::string> header;
    VCF_Reader(boost::iostreams::filtering_istream& instream);
    Parsed_VCF_Line operator()(tbb::flow_control& )const;
};