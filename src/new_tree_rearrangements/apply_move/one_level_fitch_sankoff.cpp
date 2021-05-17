#include "apply_move.hpp"
#include "src/new_tree_rearrangements/Fitch_Sankoff.hpp"
#include <cstddef>
typedef MAT::Mutations_Collection::const_iterator mut_iter;

struct Allele_Count_t {
    MAT::Mutation base;
    std::array<int, 4> count;
    int node_cnt;
    int get_position() const { return base.get_position(); }

    void operator+=(const Allele_Count_t &to_add) {
        assert(base.get_position() == to_add.base.get_position());
        for (int i = 0; i < 4; i++) {
            count[i] += to_add.count[i];
        }
        node_cnt += to_add.node_cnt;
    }

    void operator+=(const MAT::Mutation &to_add) {
        assert(base.get_position() == to_add.get_position());
        for (int i = 0; i < 4; i++) {
            if ((1 << i) & to_add.get_all_major_allele()) {
                count[i] += 1;
            }
        }
        node_cnt++;
    }

    Allele_Count_t(){};
    Allele_Count_t(const MAT::Mutation &in)
        : base(in), count({0, 0, 0, 0}), node_cnt(0) {
        (*this) += in;
    };
};
template <typename type1, typename type2>
static void merge_allele_count(const type1 &in1, const type2 &in2,
                               std::vector<Allele_Count_t> &out) {
    auto in1_iter = in1.begin();
    auto in1_end = in1.end();
#ifdef MERGE_ALLELE_CHECK_ORDER
    size_t last_in_pos = 0;
#endif
    for (const auto &other_mut : in2) {
        while (in1_iter != in1_end &&
               in1_iter->get_position() < other_mut.get_position()) {
            out.emplace_back(*in1_iter);
#ifdef MERGE_ALLELE_CHECK_ORDER
            assert(last_in_pos < in1_iter->get_position());
            last_in_pos = in1_iter->get_position();
#endif
            in1_iter++;
        }
        out.emplace_back(other_mut);
#ifdef MERGE_ALLELE_CHECK_ORDER
        assert(last_in_pos < other_mut.get_position());
        last_in_pos = other_mut.get_position();
#endif
        if (in1_iter != in1_end &&
            in1_iter->get_position() == other_mut.get_position()) {
            out.back() += *in1_iter;
            in1_iter++;
        }
    }

    while (in1_iter != in1_end) {
        out.emplace_back(*in1_iter);
#ifdef MERGE_ALLELE_CHECK_ORDER
        assert(last_in_pos < in1_iter->get_position());
        last_in_pos = in1_iter->get_position();
#endif
        in1_iter++;
    }
}

typedef std::vector<MAT::Node *>::const_iterator node_iter;
static void get_mut_count(const node_iter begin, const node_iter end,
                          std::vector<Allele_Count_t> &out) {
    if (end - begin == 3) {
        std::vector<Allele_Count_t> first;
        merge_allele_count((*begin)->mutations, (*(begin + 1))->mutations,
                           first);
        merge_allele_count((*(begin + 2))->mutations, first, out);
    } else if (end - begin == 2) {
        merge_allele_count((*begin)->mutations, (*(begin + 1))->mutations, out);
    } else if (begin == end - 1) {
        for (const auto &mut : (*begin)->mutations) {
            out.emplace_back(mut);
        }
    } else {
        assert(end - begin >= 2);
        std::vector<Allele_Count_t> first;
        std::vector<Allele_Count_t> second;
        int offset = (end - begin) / 2;
        get_mut_count(begin, begin + offset, first);
        get_mut_count(begin + offset, end, second);
        merge_allele_count(first, second, out);
    }
}

static void
rewind_ori_mut_ploytomy(MAT::Mutations_Collection &new_major_alleles_out,
                        const mut_iter &iter, bool &changed) {
    if (iter->get_all_major_allele() != iter->get_mut_one_hot() ||
        iter->get_boundary1_one_hot()) {
        changed = true;
    }
    if (iter->is_valid()) {
        new_major_alleles_out.push_back(*iter);
        new_major_alleles_out.back().set_auxillary(iter->get_mut_one_hot(), 0);
    }
}

bool get_major_allele_polytomy(MAT::Node *node,
                               MAT::Mutations_Collection &new_major_alleles_out,
                               State_Change_Collection &changed_positions) {
    std::vector<Allele_Count_t> allele_count;
    get_mut_count(node->children.begin(), node->children.end(), allele_count);
    for (auto &cnt : allele_count) {
        auto par_nuc = cnt.base.get_par_one_hot();
        int par_idx = one_hot_to_two_bit(par_nuc);
        cnt.count[par_idx] += (node->children.size() - cnt.node_cnt);
    }
    auto iter = node->mutations.begin();
    auto end = node->mutations.end();
    bool changed = false;
    std::vector<int> backward_mut;
    for (const auto &allele : allele_count) {
        while (iter != end && iter->get_position() < allele.get_position()) {
            rewind_ori_mut_ploytomy(new_major_alleles_out, iter, changed);
            iter++;
        }
        uint8_t boundary1_major_allele;
        set_state_from_cnt(allele.count, boundary1_major_allele);
        nuc_one_hot major_alleles = boundary1_major_allele & 0xf;
        if (iter != end && iter->get_position() == allele.get_position()) {
            MAT::Mutation altered = *iter;
            if (major_alleles != iter->get_all_major_allele()) {
                changed = true;
                nuc_one_hot ori_mut_nuc = altered.get_mut_one_hot();
                nuc_one_hot par_nuc = altered.get_par_one_hot();
                if (major_alleles & par_nuc) {
                    altered.set_mut_one_hot(major_alleles & par_nuc);
                } else if (!(major_alleles & ori_mut_nuc)) {
                    altered.set_mut_one_hot(major_alleles.choose_first());
                }
                if (altered.get_mut_one_hot() != ori_mut_nuc) {
                    changed_positions.emplace_back(altered, ori_mut_nuc);
                }
            }
            altered.set_auxillary(major_alleles, boundary1_major_allele >> 4);
            if ((altered.get_par_one_hot() != altered.get_all_major_allele()) ||
                (altered.get_boundary1_one_hot() &&
                 node->children.size() > 1)) {
                new_major_alleles_out.push_back(altered);
            }
            iter++;
        } else {
            nuc_one_hot par_nuc = allele.base.get_par_one_hot();
            nuc_one_hot boundary1_allele = boundary1_major_allele >> 4;
            changed |= (major_alleles != par_nuc);
            if (boundary1_allele || major_alleles != par_nuc) {
                new_major_alleles_out.push_back(allele.base);
                MAT::Mutation &altered = new_major_alleles_out.back();
                if (par_nuc & major_alleles) {
                    altered.set_mut_one_hot(par_nuc & major_alleles);
                    altered.set_auxillary(major_alleles, boundary1_allele);
                } else {
                    altered.set_mut_one_hot(major_alleles.choose_first());
                    altered.set_auxillary(major_alleles, boundary1_allele);
                    changed_positions.emplace_back(altered, par_nuc);
                }
            }
        }
    }
    while (iter != end) {
        if (iter->get_all_major_allele() != iter->get_mut_one_hot()) {
            changed = true;
        }
        if (iter->is_valid()) {
            new_major_alleles_out.push_back(*iter);
            new_major_alleles_out.back().set_auxillary(iter->get_mut_one_hot(),
                                                       0);
        }
        iter++;
    }
    return changed;
}
