#pragma once

#include "src/ripples/util/status.hpp"
#include "src/usher_graph.hpp"
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace ripples::filtration {

using CladeLineagePair = std::pair<std::string, std::string>;

struct recomb_results_row {
    using string_view = std::string_view;
    MAT::Node *recomb;
    MAT::Node *donor;
    MAT::Node *acceptor;
    string_view bp1;
    string_view bp2;
    string_view original_score;
    string_view recomb_score;
};

struct parent {
    parent(MAT::Node *node) : node(node) {}
    MAT::Node *node;
};

struct parsimony_scores {
    using string_view = std::string_view;
    parsimony_scores(string_view orig, string_view recomb);
    int original;
    int recomb;
};

class recombinant {
    using string = std::string;
    using string_view = std::string_view;
    using Mutation = MAT::Mutation;
    using MutationVec = std::vector<Mutation>;

    struct breakpoints {
        breakpoints(string_view bp1, string_view bp2);
        string bp1;
        string bp2;
    };

    struct RepresentativeSample {
        RepresentativeSample(int num_mutations, string_view id)
            : num_novel_mutations(num_mutations), node_id(id) {}
        int num_novel_mutations;
        string node_id;
    };

  public:
    recombinant(const recomb_results_row &row);

    // Getters
    string_view id() const;
    string_view donor_id() const;
    string_view acceptor_id() const;
    string_view bp1() const;
    string_view bp2() const;
    size_t num_desc(MAT::Tree &tree) const;
    size_t donor_num_desc(MAT::Tree &tree) const;
    size_t acceptor_num_desc(MAT::Tree &tree) const;
    CladeLineagePair get_clade_assignments(MAT::Tree &tree) const;
    CladeLineagePair get_donor_clade_assignments(MAT::Tree &tree) const;
    CladeLineagePair get_acceptor_clade_assignments(MAT::Tree &tree) const;
    int parsimony_original() const;
    int parsimony_improvement() const;

    // Setters
    void set_bp1(string_view bp);
    void set_bp2(string_view bp);

    string find_representative_sample(MAT::Tree &tree) const;
    string get_collapsed_mutations(MAT::Tree &tree, std::string_view node_id) const;

  private:
    breakpoints breakpoints_;
    parsimony_scores scores_;
    MAT::Node *node_;
    parent donor_;
    parent acceptor_;
};

using candidate_recombinants = std::vector<recombinant>;

}; // namespace ripples::filtration

