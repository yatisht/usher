#include "process_individual_mutation.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <vector>
extern int debug_mask;

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
static void update_par_cur_nuc(MAT::Mutations_Collection::const_iterator parent_mutation_iter,dbg_iter& debug_iter,dbg_iter& debug_end){
    rewind_mutations(parent_mutation_iter->get_position(), debug_iter,debug_end);
    if (debug_iter!=debug_end&&debug_iter->position==parent_mutation_iter->get_position()) {
                debug_iter->par_nuc.push_back(parent_mutation_iter->get_par_one_hot());
                debug_iter->major_allele.push_back(parent_mutation_iter->get_mut_one_hot()|parent_mutation_iter->get_tie_one_hot());
                debug_iter->mutation_score_change.push_back(0);
                debug_iter->count_change.emplace_back(*parent_mutation_iter);
                debug_iter->count_change.back().set_change(0, 0);
                debug_iter++;
    }
}
#endif
static void consume_parent_mutations(
    MAT::Mutations_Collection::const_iterator &parent_mutation_iter,
    const MAT::Mutations_Collection::const_iterator &parent_mutation_end,
    const Mutation_Count_Change &this_mut,int &parsimony_score_change
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
static void update_dbg_vector_score_only(
    int position, dbg_iter &debug_iter,dbg_iter &debug_end,
    int par_score_change) {
    rewind_mutations(position, debug_iter,debug_end);
    assert(debug_iter->position == position);
    nuc_one_hot parent_nuc = debug_iter->par_nuc.back();
    debug_iter->par_nuc.push_back(parent_nuc);
    debug_iter->major_allele.push_back(parent_nuc);
    debug_iter->mutation_score_change.push_back(par_score_change);
    debug_iter->count_change.push_back(debug_iter->count_change.back());
    debug_iter->count_change.back().set_change(0, 0);
    debug_iter++;
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
            mutation_iter, mutation_end, this_mut,parent_parsimony_score_change
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
                this_node_mutation.set_auxillary(0, (~parent_state)&0xf, 0);
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
    while (mutation_iter<mutation_end) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        update_par_cur_nuc(mutation_iter, debug_iter, debug_end);
#endif
        mutation_iter++;
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    while (debug_iter < debug.end()) {
        nuc_one_hot parent_nuc = debug_iter->par_nuc.back();
        debug_iter->par_nuc.push_back(parent_nuc);
        debug_iter->major_allele.push_back(parent_nuc);
        debug_iter->mutation_score_change.push_back(0);
    debug_iter->count_change.push_back(debug_iter->count_change.back());
    debug_iter->count_change.back().set_change(0, 0);
        debug_iter++;
    }
#endif
}

void get_parsimony_score_change_from_add(
    MAT::Node *node,
    const Mutation_Count_Change_Collection &children_added_mutations,
    Mutation_Count_Change_Collection &parent_added_mutations, bool is_terminal,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<MAT::Node *> &node_stack
#endif
) {
    auto parent_addable_iter = node->mutations.begin();
    auto parent_addable_end = node->mutations.end();
    for (const auto &added_child_mutation : children_added_mutations) {
        /*if(added_child_mutation.get_position()==22168){
            fputc('a', stderr);
        }*/
        while (parent_addable_end != parent_addable_iter &&
               parent_addable_iter->get_position() < added_child_mutation.get_position()) {
            if (is_terminal) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                int old_score = parsimony_score_change;
#endif
                // Adding a children with parent state
                Mutation_Count_Change temp(*parent_addable_iter);
                temp.set_change(0, parent_addable_iter->get_mut_one_hot());
                nuc_one_hot major_alleles = increment_mutation_count(
                    parent_added_mutations, *parent_addable_iter, temp,
                    parsimony_score_change, is_terminal);
                if (node->is_leaf()) {
                    assert(parent_added_mutations.empty()||parent_added_mutations.back().get_position()!=parent_addable_iter->get_position());
                }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                test_allele_out_init(node, 0, *parent_addable_iter,
                                     parsimony_score_change - old_score,
                                     major_alleles,
                                     parent_addable_iter->get_mut_one_hot(),
                                     parent_added_mutations, debug);

#endif
            }
            parent_addable_iter++;
        }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            int old_score = parsimony_score_change;
#endif
        if (parent_addable_end != parent_addable_iter &&
            parent_addable_iter->get_position() ==
                added_child_mutation.get_position()) {
            /*MAT::Mutation parent_mut_to_mod=*parent_addable_iter;
            if (node->children.size()<=1) {
                assert(!parent_mut_to_mod.get_boundary1_one_hot());
                parent_mut_to_mod.set_boundary_one_hot((~(parent_mut_to_mod.get_mut_one_hot()|parent_mut_to_mod.get_tie_one_hot())&0xf));
            }*/
            nuc_one_hot major_alleles;
            major_alleles = increment_mutation_count(
                parent_added_mutations, *parent_addable_iter,
                added_child_mutation, parsimony_score_change, is_terminal);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(node, 0, *parent_addable_iter,
                                 parsimony_score_change - old_score,
                                 major_alleles,
                                 added_child_mutation.get_incremented(),
                                 parent_added_mutations, debug);

#endif
            parent_addable_iter++;
        } else if (is_terminal) {
            if (node->children.size()<=1) {
                MAT::Mutation parent_mutation(added_child_mutation.get_position());
                //nuc_one_hot parent_state=get_parent_state(node, added_child_mutation.get_position());
                nuc_one_hot parent_state=added_child_mutation.get_par_state();
                parent_mutation.set_par_mut(parent_state, parent_state);
                parent_mutation.set_auxillary(0, (~parent_state)&0xf, 0);
                nuc_one_hot major_alleles = increment_mutation_count(
                parent_added_mutations, parent_mutation,
                added_child_mutation, parsimony_score_change, is_terminal);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(node, 0, parent_mutation,
                                 parsimony_score_change - old_score,
                                 major_alleles,
                                 added_child_mutation.get_incremented(),
                                 parent_added_mutations, debug);
#endif
            }else{
            if (!(added_child_mutation.get_par_state()&added_child_mutation.get_incremented()) ){
                parsimony_score_change++;
            }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            debug.emplace_back(added_child_mutation.get_position(),added_child_mutation.get_par_state(),added_child_mutation.get_par_state(),parsimony_score_change-old_score,added_child_mutation);
#endif
            }
        }
    }
    while (parent_addable_iter != parent_addable_end) {
        if (is_terminal) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            int old_score = parsimony_score_change;
#endif
            // Adding a children with parent state
            Mutation_Count_Change temp(*parent_addable_iter);
            temp.set_change(0, parent_addable_iter->get_mut_one_hot());
            nuc_one_hot major_alleles = increment_mutation_count(
                parent_added_mutations, *parent_addable_iter, temp,
                parsimony_score_change, is_terminal);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(node, node_stack.back(), *parent_addable_iter,
                                 parsimony_score_change - old_score,
                                 major_alleles,
                                 parent_addable_iter->get_mut_one_hot(),
                                 parent_added_mutations, debug);

#endif
        }
        parent_addable_iter++;
    }
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
            temp.set_change(parent_iter->get_mut_one_hot(), 0);
            nuc_one_hot major_alleles = decrement_mutation_count(
                parent_mutation_count_change_out, *parent_iter, temp,
                parent_parsimony_score_change, true);
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
            temp.set_change(src_mut.get_mut_one_hot()|src_mut.get_tie_one_hot(), 0);
            major_alleles = decrement_mutation_count(
                parent_mutation_count_change_out, *parent_iter, temp,
                parent_parsimony_score_change, true);
            
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
            temp.set_change(parent_iter->get_mut_one_hot(), 0);
        nuc_one_hot major_alleles = decrement_mutation_count(
            parent_mutation_count_change_out, *parent_iter, temp,
            parent_parsimony_score_change, true);

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
