#include "recombinant.hpp"

ripples::filtration::parsimony_scores::parsimony_scores(std::string_view orig,
                                                        std::string_view recomb)
    : original(std::stoi(std::string(orig))),
      recomb(std::stoi(std::string(recomb))) {}

ripples::filtration::recombinant::breakpoints::breakpoints(std::string_view bp1,
                                                           std::string_view bp2)
    : bp1(bp1), bp2(bp2) {}

ripples::filtration::recombinant::recombinant(const recomb_results_row &row)
    : breakpoints_(row.bp1, row.bp2),
      scores_(row.original_score, row.recomb_score), node_(row.recomb),
      donor_(row.donor), acceptor_(row.acceptor) {}

std::string_view ripples::filtration::recombinant::id() const {
    return node_->identifier;
}

std::string_view ripples::filtration::recombinant::donor_id() const {
    return donor_.node->identifier;
}

std::string_view ripples::filtration::recombinant::acceptor_id() const {
    return acceptor_.node->identifier;
}

std::string_view ripples::filtration::recombinant::bp1() const {
    return breakpoints_.bp1;
}

std::string_view ripples::filtration::recombinant::bp2() const {
    return breakpoints_.bp2;
}

size_t ripples::filtration::recombinant::num_desc(MAT::Tree &tree) const {
    return tree.get_num_leaves(node_);
}

size_t ripples::filtration::recombinant::donor_num_desc(MAT::Tree &tree) const {
    return tree.get_num_leaves(donor_.node);
}

size_t
ripples::filtration::recombinant::acceptor_num_desc(MAT::Tree &tree) const {
    return tree.get_num_leaves(acceptor_.node);
}

int ripples::filtration::recombinant::parsimony_original() const {
    return scores_.original;
}

int ripples::filtration::recombinant::parsimony_improvement() const {
    return scores_.original - scores_.recomb;
}

void ripples::filtration::recombinant::set_bp1(std::string_view bp) {
    breakpoints_.bp1 = bp;
}

void ripples::filtration::recombinant::set_bp2(std::string_view bp) {
    breakpoints_.bp2 = bp;
}

ripples::filtration::CladeLineagePair
ripples::filtration::recombinant::get_clade_assignments(MAT::Tree &tree) const {
    return std::make_pair(tree.get_clade_assignment(node_, 0),
                          tree.get_clade_assignment(node_, 1));
}

ripples::filtration::CladeLineagePair
ripples::filtration::recombinant::get_donor_clade_assignments(
    MAT::Tree &tree) const {
    return std::make_pair(tree.get_clade_assignment(donor_.node, 0),
                          tree.get_clade_assignment(donor_.node, 1));
}

ripples::filtration::CladeLineagePair
ripples::filtration::recombinant::get_acceptor_clade_assignments(
    MAT::Tree &tree) const {
    return std::make_pair(tree.get_clade_assignment(acceptor_.node, 0),
                          tree.get_clade_assignment(acceptor_.node, 1));
}

std::string ripples::filtration::recombinant::find_representative_sample(
    MAT::Tree &tree) const {
    // Get all the descendant nodes for internal node
    auto desc_nodes_vec = tree.get_leaves(node_->identifier);
    if (desc_nodes_vec.empty()) {
        return "None";
    }
    std::vector<RepresentativeSample> rep_samples;
    rep_samples.reserve(desc_nodes_vec.size());
    // Go through all of the internal_node descendants and find a sample
    // with no or fewest additional mutations
    for (const auto *node : desc_nodes_vec) {
        rep_samples.emplace_back(node->mutations.size(), node->identifier);
    }
    // Sort representative samples by fewest
    // additional mutations wrt given internal node
    std::sort(rep_samples.begin(), rep_samples.end(),
              [](const RepresentativeSample &a, const RepresentativeSample &b) {
                  return a.num_novel_mutations < b.num_novel_mutations;
              });
    // Return the sample with fewest additional mutations
    return rep_samples[0].node_id;
}

std::string ripples::filtration::recombinant::get_collapsed_mutations(MAT::Tree &tree, std::string_view node_id) const {
    // Get the path from node to root, including node, then reverse it to get root to node.
    std::vector<MAT::Node*> path = tree.rsearch(std::string(node_id), true);
    std::reverse(path.begin(), path.end());
    // Make a new node to accumulate all mutations on the path from root to node, collapsing multiple changes at the
    // same position including reversions.
    MAT::Node n;
    for (auto node: path) {
        for (auto mut: node->mutations) {
            n.add_mutation(mut);
        }
    }
    bool is_first = true;
    std::string mut_string;
    for (auto mut: n.mutations) {
        if (is_first) {
            is_first = false;
        } else {
            mut_string += ",";
        }
        mut_string += mut.get_string();
    }
    return mut_string;
}
