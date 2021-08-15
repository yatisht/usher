#include "apply_move.hpp"
typedef MAT::Mutations_Collection::const_iterator mut_iter;
//set the major and boundary one allele given state of its children as in Fitch Sankoff
bool get_new_mut_binary(MAT::Mutation &base, nuc_one_hot left_branch,
                        nuc_one_hot right_branch) {
    //Share any allele?
    nuc_one_hot major_allele = left_branch & right_branch;
    bool have_shared = major_allele;
    nuc_one_hot boundary1_alleles = 0;
    if (!major_allele) {
        //if not, all alleles have count of 1
        major_allele = left_branch | right_branch;
        boundary1_alleles = (~major_allele) & 0xf;
    }

    if ((left_branch | right_branch) != major_allele) {
        //share some of the allele(allele count =2), some other  not shared (allele count =1)
        boundary1_alleles = (left_branch | right_branch) & (~major_allele);
    }
    //set anxillary states
    base.set_auxillary(major_allele,boundary1_alleles);
    return have_shared;
}
//update original mutations not present in mutation vector of either children (they are shared by both)
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
        //they didn't come up in children mutation vector, so no anxillary states
        major_alleles_out.back().set_auxillary(major_allele,0);
    }
}
/**
 * @brief Set par_nuc and mut_nuc of an allele that present in either child of a node and output it
    , and copy over all valid old mutations preceeding it
 * @param[in] to_set The mutation present in either child, with all_major_allele and boundary1_allele set
 * @param[inout] iter original mutation vector iterator from last invocation
 * @param[in] end end iterator of original mutation vector
 * @param[out] changed_states changed loci with mut_nuc
 * @param[out] major_alleles_out new mutation vector to append
 * @param[out] changed whether major allele set have changed (so whether fitch sankoff on its parent node is necessary)
 * @return (void)
 */
static void set_major_allele_state_binary(
    MAT::Mutation &to_set, mut_iter &iter, const mut_iter &end,
    State_Change_Collection &changed_states,
    MAT::Mutations_Collection &major_alleles_out, bool &changed) {
    while (iter != end && iter->get_position() < to_set.get_position()) {
        //original allele not present in both children
        place_ori_mutation(iter, major_alleles_out, changed);
        iter++;
    }
    //if originally there were no mutation at this loci, then the input to_set have valid par_nuc
    nuc_one_hot par_nuc = to_set.get_par_one_hot();
    nuc_one_hot old_state = par_nuc;
    nuc_one_hot major_allele = to_set.get_all_major_allele();
    if (iter != end && iter->get_position() == to_set.get_position()) {
        //there were mutations at this loci, so the original mutation have the correct par_nuc
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
        //if cannot follow parent then prefer old state to minimize changes to te tree
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
//specialized backward patching routine for binary nodes (they do not need two step count + find max)
bool get_major_allele_binary(const MAT::Node *node,
                             MAT::Mutations_Collection &out,
                             State_Change_Collection &changed_states) {
    //left node iterator
    auto iter = node->children[0]->mutations.begin();
    auto end = node->children[0]->mutations.end();
    //original mutation vector iterator
    mut_iter ori_mut_iter = node->mutations.begin();
    mut_iter ori_mut_end = node->mutations.end();
    bool changed = false;
    //iterating on right child mutation
    for (const auto &right_mut : node->children[1]->mutations) {
        while (iter != end && iter->get_position() < right_mut.get_position()) {
            //only present in left child, right is par_nuc
            MAT::Mutation mut = *iter;
            get_new_mut_binary(mut, iter->get_all_major_allele(),
                               iter->get_par_one_hot());
            set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                          changed_states, out, changed);
            iter++;
        }
        MAT::Mutation mut = right_mut;
        if (iter != end && iter->get_position() == right_mut.get_position()) {
            //left right children match
            get_new_mut_binary(mut, iter->get_all_major_allele(),
                               right_mut.get_all_major_allele());
            iter++;
        } else {
            //only right, left parental allele
            get_new_mut_binary(mut, right_mut.get_par_one_hot(),
                               right_mut.get_all_major_allele());
        }
        set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                      changed_states, out, changed);
    }
    while (iter != end) {
        //left only
        MAT::Mutation mut = *iter;
        get_new_mut_binary(mut, iter->get_all_major_allele(),
                           iter->get_par_one_hot());
        set_major_allele_state_binary(mut, ori_mut_iter, ori_mut_end,
                                      changed_states, out, changed);
        iter++;
    }
    while (ori_mut_iter != ori_mut_end) {
        //continue copying over original mutations
        place_ori_mutation(ori_mut_iter, out, changed);
        ori_mut_iter++;
    }
    return changed;
}
