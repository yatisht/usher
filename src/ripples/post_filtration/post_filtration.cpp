#include "post_filtration.hpp"
#include "recomb_writer.hpp"
#include "src/ripples/util/text_parser.hpp"
#include <stdexcept>
#include <charconv>
#include <limits>

namespace {

// For parsing breakpoint ranges to see if they overlap and can be merged.

struct IntRange {
    unsigned int start;
    unsigned int end;
};

// Parse a range string of the form "(uint,uint)" or "(uint,GENOME_SIZE)".
// genome_size is substituted when "GENOME_SIZE" appears as the second value.
// Throws std::invalid_argument if string is malformed.
IntRange parse_breakpoint_range(std::string_view s, unsigned int genome_size) {
    // Must start with '(' and end with ')'
    if (s.size() < 5 || s.front() != '(' || s.back() != ')')
        throw std::invalid_argument("Breakpoint range must start and end with parentheses");
    s = s.substr(1, s.size() - 2);  // strip parens -> "uint,uint" or "uint,GENOME_SIZE"
    // Split on comma
    auto comma = s.find(',');
    if (comma == std::string_view::npos)
        throw std::invalid_argument("Breakpoint range is missing comma");
    std::string_view lhs = s.substr(0, comma);
    std::string_view rhs = s.substr(comma + 1);
    // Parse left side (must be a plain uint)
    IntRange range;
    auto [ptr1, ec1] = std::from_chars(lhs.data(), lhs.data() + lhs.size(), range.start);
    if (ec1 != std::errc{} || ptr1 != lhs.data() + lhs.size())
        throw std::invalid_argument("Breakpoint range: invalid start coordinate");
    // Parse right side: GENOME_SIZE or plain uint
    if (rhs == "GENOME_SIZE") {
        range.end = genome_size;
    } else {
        auto [ptr2, ec2] = std::from_chars(rhs.data(), rhs.data() + rhs.size(), range.end);
        if (ec2 != std::errc{} || ptr2 != rhs.data() + rhs.size())
            throw std::invalid_argument("Breakpoint range: invalid end coordinate");
    }
    return range;
}

// Return true if r1 completely contains, or is the same as, r2.
bool range_contains(IntRange& r1, IntRange& r2) {
    return (r2.start >= r1.start && r2.end <= r1.end);
}

// Convert IntRange back to a breakpoint range string with symbolic "GENOME_SIZE"
std::string format_breakpoint_range(const IntRange& r, unsigned int genome_size) {
  if (r.end == genome_size) {
      return '(' + std::to_string(r.start) + ",GENOME_SIZE)";
  } else {
      return '(' + std::to_string(r.start) + ',' + std::to_string(r.end) + ')';
  }
}

} // namespace

ripples::filtration::post_filtration::post_filtration(
    MAT::Tree &tree, const std::string &recomb_file)
    : tree_(tree) {
    // Column indices in ripples-fast 'recombination.tsv' file
    static constexpr int RECOMB_ID_COL{0};
    static constexpr int BREAKPOINT_1_COL{1};
    static constexpr int BREAKPOINT_2_COL{2};
    static constexpr int DONOR_ID_COL{3};
    static constexpr int ACCEPTOR_ID_COL{6};
    static constexpr int ORIG_PARSIMONY_COL{9};
    static constexpr int RECOMB_PARSIMONY_COL{11};

    // Large integer for the symbolic "GENOME_SIZE" used in breakpoint range strings
    static constexpr unsigned int GENOME_SIZE = std::numeric_limits<unsigned int>::max();

    text_parser parser(recomb_file);
    // Skip over header
    parser.next_line();

    auto process_line = [&](const text_parser &parser) -> recombinant {
        auto recomb_id = parser.get_value(RECOMB_ID_COL);
        auto bp1 = parser.get_value(BREAKPOINT_1_COL);
        auto bp2 = parser.get_value(BREAKPOINT_2_COL);
        auto orig_score = parser.get_value(ORIG_PARSIMONY_COL);
        auto recomb_score = parser.get_value(RECOMB_PARSIMONY_COL);
        auto donor_id = parser.get_value(DONOR_ID_COL);
        auto acceptor_id = parser.get_value(ACCEPTOR_ID_COL);
        MAT::Node *recomb = tree_.get_node(string{recomb_id});
        MAT::Node *donor = tree_.get_node(string{donor_id});
        MAT::Node *acceptor = tree_.get_node(string{acceptor_id});

        recomb_results_row row{recomb, donor,      acceptor,    bp1,
                               bp2,    orig_score, recomb_score};
        return recombinant(row);
    };

    std::vector<recombinant> all_recombinants;
    for (; !parser.done(); parser.next_line()) {
        all_recombinants.emplace_back(process_line(parser));
    }

    // Merge recombinants that have the same recomb_id, donor_id and acceptor_id and overlapping breakpoint ranges,
    // and filter to keep only those with the highest parsimony score improvement per node triplet.
    // Sort the merged and filtered recombinants by parsimony improvement (descending).
    auto make_triplet_key = [&](const recombinant& r) -> std::string {
        return std::string(r.id()) + '\t' +
               std::string(r.donor_id()) + '\t' +
               std::string(r.acceptor_id());
    };
    auto merge_recombinant_breakpoints = [&](recombinant& old_r, recombinant& r) -> bool {
        // If both breakpoint ranges meet the condition of old containing new or new containing old (i.e. overlap,
        // but not incomplete overlap), then merge both breakpoint ranges into old_r's.
        bool did_merge = false;
        IntRange r_bp1 = parse_breakpoint_range(r.bp1(), GENOME_SIZE);
        IntRange old_r_bp1 = parse_breakpoint_range(old_r.bp1(), GENOME_SIZE);
        if (range_contains(r_bp1, old_r_bp1) || range_contains(old_r_bp1, r_bp1)) {
            IntRange r_bp2 = parse_breakpoint_range(r.bp2(), GENOME_SIZE);
            IntRange old_r_bp2 = parse_breakpoint_range(old_r.bp2(), GENOME_SIZE);
            if (range_contains(r_bp2, old_r_bp2) || range_contains(old_r_bp2, r_bp2)) {
                IntRange merged_bp1 = { std::min(r_bp1.start, old_r_bp1.start), std::max(r_bp1.end, old_r_bp1.end) };
                IntRange merged_bp2 = { std::min(r_bp2.start, old_r_bp2.start), std::max(r_bp2.end, old_r_bp2.end) };
                old_r.set_bp1(format_breakpoint_range(merged_bp1, GENOME_SIZE));
                old_r.set_bp2(format_breakpoint_range(merged_bp2, GENOME_SIZE));
                did_merge = true;
            }
        }
        return did_merge;
    };
    auto merge_recombinants = [&](std::vector<recombinant>&& all_recombinants) -> std::vector<recombinant> {
        std::unordered_map<std::string,std::vector<recombinant>> recomb_map;
        recomb_map.reserve(all_recombinants.size()); // this is a big overestimate, I expect mostly discards and merges
        for (recombinant& r : all_recombinants) {
            std::string key = make_triplet_key(r);
            auto it = recomb_map.find(key);
            if (it == recomb_map.end()) {
                // First recombinant encountered for this node triplet
                std::vector<recombinant> vec;
                vec.push_back(std::move(r));
                recomb_map.emplace(std::move(key), std::move(vec));
            } else {
                // At least one recombinant for this node triplet has already been found; depending on parsimony
                // improvement and breakpoint ranges, replace, discard, merge, or add.
                bool add_unmergeable = false;
                for (recombinant& old_r: it->second) {
                    if (r.parsimony_improvement() > old_r.parsimony_improvement()) {
                        // New recombinant has a greater parsimony improvement; replace the old vector.
                        it->second = {};
                        it->second.push_back(std::move(r));
                        break;
                    } else if (r.parsimony_improvement() < old_r.parsimony_improvement()) {
                        // New recombinant has a smaller improvement, discard it.
                        break;
                    } else {
                        // Equal parsimony improvement; merge breakpoint ranges if possible.
                        if (merge_recombinant_breakpoints(old_r, r)) {
                            add_unmergeable = false;
                            break;
                        } else {
                            add_unmergeable = true;
                        }
                    }
                }
                if (add_unmergeable) {
                    it->second.push_back(std::move(r));
                }
            }
        }
        // Extract merged recombinants from recomb_map and sort by parsimony improvement (descending)
        std::vector<recombinant> merged_recombinants;
        merged_recombinants.reserve(recomb_map.size() * 2);
        for (auto& [key, vec] : recomb_map) {
            for (recombinant& r: vec) {
                merged_recombinants.push_back(std::move(r));
            }
        }
        std::sort(merged_recombinants.begin(), merged_recombinants.end(),
                  [](const recombinant& a, const recombinant& b) {
                      return a.parsimony_improvement() > b.parsimony_improvement(); // > for descending
                  });
        return merged_recombinants;
    };
    std::vector<recombinant> merged_recombinants = merge_recombinants(std::move(all_recombinants));

    // Filter the sorted recombinants to keep only the top n (including ties) for each recomb_id
    auto filter_recombinants = [&](std::vector<recombinant>&& recombinants, int top_n) -> std::vector<recombinant> {
        std::unordered_map<std::string_view,int> result_counts;
        result_counts.reserve(recombinants.size());
        std::unordered_map<std::string_view,int> min_scores;
        min_scores.reserve(recombinants.size());
        std::vector<recombinant> filtered_recombinants;
        filtered_recombinants.reserve(recombinants.size());
        for (recombinant& r: recombinants) {
            bool keep_this = true;
            const std::string_view recomb_id = r.id();
            int count = ++result_counts[recomb_id];
            if (count == top_n) {
                // This is top nth place, note the score in case there are ties
                min_scores[recomb_id] = r.parsimony_improvement();
            } else if (count > top_n && r.parsimony_improvement() < min_scores[recomb_id]) {
                keep_this = false;
            }
            if (keep_this) {
                filtered_recombinants.push_back(std::move(r));
            }
        }
        return filtered_recombinants;
    };
    recombinants_ = filter_recombinants(std::move(merged_recombinants), 3);
}

ripples::server::Status ripples::filtration::post_filtration::write(
    const std::string &outfilepath) {
    recomb_writer writer(tree_, outfilepath);
    return writer.write(recombinants_);
}

