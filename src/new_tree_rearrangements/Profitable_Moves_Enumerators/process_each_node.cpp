#include "process_individual_mutation.hpp"
#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <vector>

void consume_parent_mutations(
    MAT::Mutations_Collection::const_iterator &parent_mutation_iter,
    const MAT::Mutations_Collection::const_iterator &parent_mutation_end,
    const Mutation_Count_Change &this_mut
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,dbg_iter &debug_iter,dbg_iter &debug_end
#endif
    ) {
    while (parent_mutation_iter != parent_mutation_end &&
           (*parent_mutation_iter) < this_mut) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        update_par_cur_nuc(parent_mutation_iter, debug_iter, debug_end);
        #endif
        parent_mutation_iter++;
    }
}

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
void intermediate_node_tail_debug_update(
    std::vector<state_change_hist_dbg> &debug,
    MAT::Mutations_Collection::const_iterator &mutation_iter,
    MAT::Mutations_Collection::const_iterator &mutation_end,
    std::vector<state_change_hist_dbg>::iterator &debug_iter,
    std::vector<state_change_hist_dbg>::iterator &debug_end) {
    while (mutation_iter < mutation_end) {
        update_par_cur_nuc(mutation_iter, debug_iter, debug_end);
        mutation_iter++;
    }
    while (debug_iter < debug.end()) {
        nuc_one_hot parent_nuc = debug_iter->par_nuc.back();
        debug_iter->par_nuc.push_back(parent_nuc);
        debug_iter->major_allele.push_back(parent_nuc);
        debug_iter->mutation_score_change.push_back(0);
        debug_iter->count_change.push_back(debug_iter->count_change.back());
        debug_iter->count_change.back().set_change(0, 0, 0, true);
        debug_iter++;
    }
}
#endif

void get_intermediate_nodes_mutations(
    const MAT::Node *node,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<MAT::Node *> &node_stack
#endif
) {
    MAT::Mutations_Collection::const_iterator mutation_iter =
        node->mutations.begin();
    MAT::Mutations_Collection::const_iterator mutation_end =
        node->mutations.end();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    auto debug_iter = debug.begin();
    auto debug_end=debug.end();
#endif
    for (const auto &this_mut : this_node_mutation_count_change) {
        /*if(this_mut.get_position()==241){
            fputc('a', stderr);
        }*/
        consume_parent_mutations(
            mutation_iter, mutation_end, this_mut
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_iter,debug_end
#endif
            );
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        int old_score_change = parent_parsimony_score_change;
#endif
        if (mutation_iter != mutation_end &&
            mutation_iter->get_position() == this_mut.get_position()) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            int old_score_change = parent_parsimony_score_change;
#endif
            nuc_one_hot major_alleles = decrement_increment_mutation_count(
                *mutation_iter, this_mut, parent_node_mutation_count_change,
                parent_parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_count_out(node, *mutation_iter, major_alleles, 0,
                                  parent_node_mutation_count_change,
                                  parent_parsimony_score_change -
                                      old_score_change,
                                  debug_iter,debug_end, node_stack);
#endif
            mutation_iter++;
        } else {
            if (node->children.size()<=1) {
                MAT::Mutation this_node_mutation(this_mut.get_position());
                nuc_one_hot parent_state=this_mut.get_par_state();
                this_node_mutation.set_par_mut(parent_state, parent_state);
                this_node_mutation.set_auxillary(parent_state, (~parent_state)&0xf);
                            nuc_one_hot major_alleles = decrement_increment_mutation_count(
                this_node_mutation, this_mut, parent_node_mutation_count_change,
                parent_parsimony_score_change);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_count_out(node, this_node_mutation, major_alleles, 0,
                                  parent_node_mutation_count_change,
                                  parent_parsimony_score_change -
                                      old_score_change,
                                  debug_iter,debug_end, node_stack);
#endif
            }else{
            parent_parsimony_score_change += this_mut.get_default_change_internal();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            update_dbg_vector_score_only(this_mut.get_position(), debug_iter,debug_end,
                                         parent_parsimony_score_change -
                                             old_score_change);
#endif
}
        }
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    intermediate_node_tail_debug_update(debug, mutation_iter, mutation_end, debug_iter, debug_end);
#endif
}

bool get_parsimony_score_change_from_add(
    MAT::Node *node,
    const Mutation_Count_Change_Collection &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<MAT::Node *> &node_stack
#endif
) {
    bool have_not_shared=false;
    auto sibling_addable_iter = node->mutations.begin();
    auto sibling_addable_end = node->mutations.end();
    for (const auto &added_child_mutation : children_added_mutations) {
        /*if(added_child_mutation.get_position()==23635){
            fputc('ab', stderr);
        }*/
        while (sibling_addable_end != sibling_addable_iter &&
               sibling_addable_iter->get_position() < added_child_mutation.get_position()) {
                // Adding a children with parent state
                nuc_one_hot major_alleles = sibling_addable_iter->get_all_major_allele()&sibling_addable_iter->get_par_one_hot();
                if (!major_alleles) {
                    major_alleles=sibling_addable_iter->get_all_major_allele()|sibling_addable_iter->get_par_one_hot();
                    have_not_shared=true;
                    parsimony_score_change++;
                }
                register_change_from_new_state(parent_added_mutations, 0, *sibling_addable_iter, major_alleles);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                test_allele_out_init(node, 0, *sibling_addable_iter,
                                     0,
                                     major_alleles,
                                     sibling_addable_iter->get_mut_one_hot(),
                                     parent_added_mutations, debug);

#endif
            
            sibling_addable_iter++;
        }
int old_parsimony_score=parsimony_score_change;
        if (sibling_addable_end != sibling_addable_iter &&
            sibling_addable_iter->get_position() ==
                added_child_mutation.get_position()) {

            nuc_one_hot major_alleles;
            int new_score=get_new_major_allele_binary_node(sibling_addable_iter->get_all_major_allele(), added_child_mutation.get_incremented(), major_alleles);
            have_not_shared=have_not_shared||(new_score);
            parsimony_score_change+=new_score;
            register_change_from_new_state(parent_added_mutations, new_score, *sibling_addable_iter, major_alleles);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(node, 0, *sibling_addable_iter,
                                 parsimony_score_change-old_parsimony_score,
                                 major_alleles,
                                 added_child_mutation.get_incremented(),
                                 parent_added_mutations, debug);

#endif
            sibling_addable_iter++;
        } else {
            nuc_one_hot major_allele;
            int new_score=get_new_major_allele_binary_node(added_child_mutation.get_par_state(), added_child_mutation.get_incremented(), major_allele);
            have_not_shared=have_not_shared||(new_score);
            parsimony_score_change+=new_score;
            if (major_allele!=added_child_mutation.get_par_state()) {
                parent_added_mutations.emplace_back(added_child_mutation,major_allele);
            }

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            debug.emplace_back(added_child_mutation.get_position(),added_child_mutation.get_par_state(),major_allele,parsimony_score_change-old_parsimony_score,added_child_mutation);
#endif

        }
    }
    while (sibling_addable_iter != sibling_addable_end) {

                nuc_one_hot major_alleles = sibling_addable_iter->get_all_major_allele()&sibling_addable_iter->get_par_one_hot();
                if (!major_alleles) {
                    major_alleles=sibling_addable_iter->get_all_major_allele()|sibling_addable_iter->get_par_one_hot();
                    have_not_shared=true;
                    parsimony_score_change++;
                }
                register_change_from_new_state(parent_added_mutations, 0, *sibling_addable_iter, major_alleles);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                test_allele_out_init(node, 0, *sibling_addable_iter,
                                     0,
                                     major_alleles,
                                     sibling_addable_iter->get_mut_one_hot(),
                                     parent_added_mutations, debug);

#endif

        sibling_addable_iter++;
    }
    return have_not_shared;
}

/*
parent of parent (not involved in this function, but the output
parent_mutation_count_change is used by it to adjust its mutation list)
|
parent (parent_mutation_count_change reports what happens to be mutation count
at each loci for this node)
|
src (being removed, amount to decrement count of all its mutations)
*/
void get_parent_altered_remove(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    //,New_Tie_Collection_t &decrease_if_match_parent_state,
    const MAT::Node *src, int &parent_parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<MAT::Node *> &node_stack
#endif
) {
    auto parent_iter = src->parent->mutations.begin();
    auto parent_end = src->parent->mutations.end();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
#endif
    const MAT::Mutations_Collection &this_node_altered = src->mutations;
    for (const auto &src_mut : this_node_altered) {
        /*if (src_mut.get_position()==186) {
            fputc('a',stderr);
        }*/
        while (parent_iter != parent_end && *parent_iter < src_mut) {
// This is removing major/tie mutations that are present in parent
// node
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            int old_score_change = parent_parsimony_score_change;
#endif
            Mutation_Count_Change temp(*parent_iter);
            temp.set_change(parent_iter->get_mut_one_hot(), 0,0,true);
            nuc_one_hot major_alleles = decrement_mutation_count(
                parent_mutation_count_change_out, *parent_iter, temp,
                parent_parsimony_score_change);
                parent_parsimony_score_change--;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(
                src->parent, src, *parent_iter,
                parent_parsimony_score_change - old_score_change, major_alleles,
                0, parent_mutation_count_change_out, debug);

#endif
            parent_iter++;
        }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            int old_score_change = parent_parsimony_score_change;
#endif
        if (parent_iter != parent_end&&parent_iter->get_position() == src_mut.get_position()) {
            Mutation_Count_Change temp(src_mut);
            nuc_one_hot major_alleles;
            temp.set_change(src_mut.get_all_major_allele(), 0,0);
            major_alleles = decrement_mutation_count(
                parent_mutation_count_change_out, *parent_iter, temp,
                parent_parsimony_score_change);
                parent_parsimony_score_change--;
            
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(
                src->parent, src, *parent_iter,
                parent_parsimony_score_change - old_score_change, major_alleles,
                0, parent_mutation_count_change_out, debug);

#endif
            parent_iter++;
        }else{
            assert(parent_iter == parent_end||parent_iter->get_position() > src_mut.get_position());
            if (src_mut.is_valid()) {
                parent_parsimony_score_change--;
            }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            debug.emplace_back(src_mut.get_position(),src_mut.get_par_one_hot(),src_mut.get_par_one_hot(),parent_parsimony_score_change - old_score_change,Mutation_Count_Change());
#endif
        }

        assert(parent_iter==parent_end||src_mut < (*parent_iter));
        // This is removing minor alleles, not interesting
        // They are not not in the state of parent,
        // as they needs to be mutated to reach that state in this_node
    }
    while (parent_iter != parent_end) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        int old_score_change = parent_parsimony_score_change;
#endif
            Mutation_Count_Change temp(*parent_iter);
            temp.set_change(parent_iter->get_mut_one_hot(), 0,0,true);
        nuc_one_hot major_alleles = decrement_mutation_count(
            parent_mutation_count_change_out, *parent_iter, temp,
            parent_parsimony_score_change);
            parent_parsimony_score_change--;

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        test_allele_out_init(src->parent, src, *parent_iter,
                             parent_parsimony_score_change - old_score_change,
                             major_alleles, 0, parent_mutation_count_change_out,
                             debug);

#endif
        // register_new_tied_mutations(decrease_if_match_parent_state,major_alleles,
        // *parent_iter,parent_parsimony_score_change);
        parent_iter++;
    }
    // return parent_parsimony_score_change;
}
