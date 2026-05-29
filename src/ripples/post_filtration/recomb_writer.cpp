#include "recomb_writer.hpp"

ripples::filtration::recomb_writer::recomb_writer(
    MAT::Tree &tree, const std::string &outfilepath)
    : tree_(tree), log_(outfilepath) {}

void ripples::filtration::recomb_writer::write(const recombinant &recomb) {
    const auto &[recomb_clade, recomb_lineage] =
        recomb.get_clade_assignments(tree_);
    const auto &[donor_clade, donor_lineage] =
        recomb.get_donor_clade_assignments(tree_);
    const auto &[acceptor_clade, acceptor_lineage] =
        recomb.get_acceptor_clade_assignments(tree_);

    log_ << recomb.id() << "\t" << recomb.donor_id() << "\t"
         << recomb.acceptor_id() << "\t" << recomb.num_desc(tree_) << "\t"
         << recomb.donor_num_desc(tree_) << "\t"
         << recomb.acceptor_num_desc(tree_) << "\t" << recomb.bp1() << "\t"
         << recomb.bp2() << "\t" << recomb_clade << "\t" << recomb_lineage
         << "\t" << donor_clade << "\t" << donor_lineage << "\t"
         << acceptor_clade << "\t" << acceptor_lineage << "\t"
         << recomb.find_representative_sample(tree_) << "\t"
         << recomb.parsimony_original() << "\t"
         << recomb.parsimony_improvement() << "\t"
         << recomb.get_collapsed_mutations(tree_, recomb.id()) << "\t"
         << recomb.get_collapsed_mutations(tree_, recomb.donor_id()) << "\t"
         << recomb.get_collapsed_mutations(tree_, recomb.acceptor_id()) << "\n";
}

ripples::server::Status ripples::filtration::recomb_writer::write(
    const candidate_recombinants &candidates) {
    if (!log_.is_open()) {
        return Status{ripples::server::error_t::CANNOT_OPEN_LOG};
    }
    log_ << HEADER;
    for (const auto &recomb : candidates) {
        write(recomb);
    }
    return Status{ripples::server::error_t::NO_ERROR};
}

