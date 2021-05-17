#include "apply_move.hpp"
typedef MAT::Mutations_Collection::const_iterator mut_iter;
bool get_new_mut_binary(MAT::Mutation &base, nuc_one_hot left_branch,
                        nuc_one_hot right_branch) {
    nuc_one_hot major_allele = left_branch & right_branch;
    bool have_shared = major_allele;
    nuc_one_hot boundary1_alleles = 0;
    if (!major_allele) {
        major_allele = left_branch | right_branch;
        boundary1_alleles = (~major_allele) & 0xf;
    }

    if ((left_branch | right_branch) != major_allele) {
        boundary1_alleles = (left_branch | right_branch) & (~major_allele);
    }
    base.set_children(boundary1_alleles, major_allele, left_branch,
                      right_branch);
    return have_shared;
}
void place_ori_mutation(mut_iter &iter,
                        MAT::Mutations_Collection &major_alleles_out,
                        bool &changed) {
    if (iter->get_all_major_allele() != iter->get_mut_one_hot()) {
        // Originally there was some other alleles,now there isn't
        changed = true;
    }
    if (iter->is_valid()) {
        major_alleles_out.push_back(*iter);
        nuc_one_hot major_allele = iter->get_mut_one_hot();
        major_alleles_out.back().set_children(0, major_allele, major_allele,
                                              major_allele);
    }
}
static void set_major_allele_state_binary(
    MAT::Mutation &to_set, mut_iter &iter, const mut_iter &end,
    State_Change_Collection &changed_states,
    MAT::Mutations_Collection &major_alleles_out, bool &changed) {
    while (iter != end && iter->get_position() < to_set.get_position()) {
        place_ori_mutation(iter, major_alleles_out, changed);
        iter++;
    }
    nuc_one_hot par_nuc = to_set.get_par_one_hot();
    nuc_one_hot old_state = par_nuc;
    nuc_one_hot major_allele = to_set.get_all_major_allele();
    if (iter != end && iter->get_position() == to_set.get_position()) {
        par_nuc = iter->get_par_one_hot();
        old_state = iter->get_mut_one_hot();
        if (iter->get_all_major_allele() != major_allele) {
            changed = true;
        }
        iter++;
    } else if (par_nuc != major_allele) {
        changed = true;
    }
    nuc_one_hot new_state = major_allele & par_nuc;
    if (!new_state) {
        if (major_allele & old_state) {
            new_state = old_state;
        } else {
            new_state = major_allele.choose_first();
        }
    }
    to_set.set_par_mut(par_nuc, new_state);
    if (new_state != par_nuc || major_allele != new_state ||
        to_set.get_boundary1_one_hot()) {
        major_alleles_out.push_back(to_set);
    }
    if (new_state != old_state) {
        changed_states.emplace_back(to_set, old_state);
    }
}
bool get_major_allele_binary(const MAT::Node *node,
                             MAT::Mutations_Collection &out,
                             State_Change_Collection &changed_states) {
    auto iter = node->children[0]->mutations.begin();
    auto end = node->children[0]->mutations.end();
    mut_iter ori_mut_iter = node->mutations.begin();
    mut_iter ori_mut_end = node->mutations.end();
    bool changed = false;
    for (const auto &right_mut : node->children[1]->mutations) {
        while (iter != end && iter->get_position() < right_mut.get_position()) {
            MAT::Mutation mut = *iter;
            get_new_mut_binary(mut, iter->get_all_major_allele(),
                               iter->get_par_one_hot());
            set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                          changed_states, out, changed);
            iter++;
        }
        MAT::Mutation mut = right_mut;
        if (iter != end && iter->get_position() == right_mut.get_position()) {
            get_new_mut_binary(mut, iter->get_all_major_allele(),
                               right_mut.get_all_major_allele());
            iter++;
        } else {
            get_new_mut_binary(mut, right_mut.get_par_one_hot(),
                               right_mut.get_all_major_allele());
        }
        set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                      changed_states, out, changed);
    }
    while (iter != end) {
        MAT::Mutation mut = *iter;
        get_new_mut_binary(mut, iter->get_all_major_allele(),
                           iter->get_par_one_hot());
        set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                      changed_states, out, changed);
        iter++;
    }
    while (ori_mut_iter != ori_mut_end) {
        place_ori_mutation(ori_mut_iter, out, changed);
        ori_mut_iter++;
    }
    return changed;
}
