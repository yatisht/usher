#pragma once

#include "recombinant.hpp"
#include <fstream>

namespace ripples::filtration {

class recomb_writer {
  public:
    using string = std::string;
    using string_view = std::string_view;
    using Status = ripples::server::Status;

    recomb_writer(MAT::Tree &tree, const string &outfilepath);
    Status write(const candidate_recombinants &candidates);

  private:
    void write(const recombinant &recomb);

    static constexpr string_view HEADER{
        "recomb_node_id\tdonor_node_id\tacceptor_node_id\trecombinant_num_"
        "desc\tdonor_num_desc\tacceptor_num_desc\tbreakpoint interval "
        "1\tbreakpoint interval 2\trecombinant clade\trecombinant "
        "lineage\tdonor clade\tdonor lineage\tacceptor clade\tacceptor "
        "lineage\trepresentative descendant\toriginal parsimony "
        "score\tparsimony score improvement\trecombinant mutations\t"
        "donor mutations\tacceptor mutations\n"};

    MAT::Tree &tree_;
    std::ofstream log_;
};

}; // namespace ripples::filtration

