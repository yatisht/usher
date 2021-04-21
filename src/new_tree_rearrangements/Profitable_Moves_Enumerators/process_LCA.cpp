#include "process_individual_mutation.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT

static void set_LCA_par_score_change(
    const Mutation_Count_Change &change,
    std::vector<state_change_hist_dbg> &debug,
    std::vector<LCA_merged_states>::const_iterator &in,
    const std::vector<LCA_merged_states>::const_iterator &begin,
    const std::vector<LCA_merged_states>::const_iterator &end,
    int score_change) {
    int position = change.get_position();
    while (in != end && in->position < position) {
        assert(debug.empty() || debug.back().position < in->position);
        debug.emplace_back(in->position, in->par_allele,in->par_allele, 0,Mutation_Count_Change());
        in++;
    }
    assert((in==begin|| (in-1)->position <= position)&&(debug.empty()||debug.back().position<position));
    if (in!=end&&in->position==position) {
        assert(change.get_par_state()==in->par_allele);
        in++;
    }
    debug.emplace_back(position,change.get_par_state(),change.get_par_state(),score_change,Mutation_Count_Change());
}
#endif
static void LCA_no_match(
    const MAT::Node *LCA, bool is_src_terminal, bool is_dst_terminal,
    MAT::Mutations_Collection::const_iterator &LCA_mutation_iter,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
    const std::vector<LCA_merged_states>::const_iterator &end,
    const MAT::Node *src_branch, const MAT::Node *dst_branch
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    int old_score = parsimony_score_change;
#endif
/*if (LCA_mutation_iter->get_position()==241) {
    fputc('a', stderr);
}*/
    nuc_one_hot major_allele=LCA_mutation_iter->get_mut_one_hot()|LCA_mutation_iter->get_tie_one_hot();
    Mutation_Count_Change temp(*LCA_mutation_iter);
    temp.set_change(LCA_mutation_iter->get_mut_one_hot(),
                    LCA_mutation_iter->get_mut_one_hot());
    if (is_dst_terminal) {
        major_allele = increment_mutation_count(
            LCA_parent_mutation_count_change_out, *LCA_mutation_iter, temp,
            parsimony_score_change, true);
    } else if (is_src_terminal) {
        major_allele = decrement_mutation_count(
            LCA_parent_mutation_count_change_out, *LCA_mutation_iter, temp,
            parsimony_score_change, true);
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    test_allele_count_out_LCA(
        LCA, src_branch, is_src_terminal, dst_branch,
        is_dst_terminal ? (uint8_t)LCA_mutation_iter->get_mut_one_hot()
                        : (uint8_t)0,
        *LCA_mutation_iter, major_allele, parsimony_score_change - old_score,
        LCA_parent_mutation_count_change_out, debug, debug_in_iter, end);
#endif
    LCA_mutation_iter++;
}
static void LCA_dst_match(const MAT::Mutation* LCA_mutation_iter,
        const Mutation_Count_Change_Collection::const_iterator &dst_add_count_iter,
        const MAT::Node *LCA, bool is_src_terminal, bool is_dst_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end,
        const MAT::Node *src_branch, const MAT::Node *dst_branch,int old_score
#endif
){
            nuc_one_hot major_allele;
        if (is_dst_terminal) {
            major_allele = increment_mutation_count(
                LCA_parent_mutation_count_change_out, *LCA_mutation_iter,
                *dst_add_count_iter, parsimony_score_change, true);
        } else if (is_src_terminal) {
            Mutation_Count_Change temp(*LCA_mutation_iter);
            temp.set_change(LCA_mutation_iter->get_mut_one_hot(),
                            0);
            major_allele = dbl_inc_dec_mutations(
                *dst_add_count_iter, false, temp, true, *LCA_mutation_iter,
                parsimony_score_change, LCA_parent_mutation_count_change_out);
        } else {
            major_allele = decrement_increment_mutation_count(
                *LCA_mutation_iter, *dst_add_count_iter,
                LCA_parent_mutation_count_change_out, parsimony_score_change);
        }

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        test_allele_count_out_LCA(
            LCA, src_branch, is_src_terminal, dst_branch, is_dst_terminal?dst_add_count_iter->get_incremented():(nuc_one_hot)0, *LCA_mutation_iter,
            major_allele, parsimony_score_change - old_score,
            LCA_parent_mutation_count_change_out, debug, debug_in_iter, end);
#endif
}
static void
LCA_dst(MAT::Mutations_Collection::const_iterator &LCA_mutation_iter,
        MAT::Mutations_Collection::const_iterator LCA_mutation_end,
        const Mutation_Count_Change_Collection::const_iterator &dst_add_count_iter,
        const MAT::Node *LCA, bool is_src_terminal, bool is_dst_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end,
        const MAT::Node *src_branch, const MAT::Node *dst_branch
#endif
) {
    while (LCA_mutation_iter != LCA_mutation_end&&LCA_mutation_iter->get_position() <
               dst_add_count_iter->get_position()
           ) {
        LCA_no_match(LCA, is_src_terminal, is_dst_terminal, LCA_mutation_iter,
                     LCA_parent_mutation_count_change_out,
                     parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                     , debug, debug_in_iter, end,
                     src_branch, dst_branch
#endif
                    );
    }

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    int old_score = parsimony_score_change;
#endif
/*if (LCA_mutation_iter->get_position()==241) {
    fputc('a',stderr);
}*/
    if (LCA_mutation_iter != LCA_mutation_end &&
        dst_add_count_iter->get_position() ==
            LCA_mutation_iter->get_position()) {
        LCA_dst_match(&(*LCA_mutation_iter), dst_add_count_iter, LCA, is_src_terminal, is_dst_terminal, LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        , debug, debug_in_iter, begin, end, src_branch, dst_branch, old_score
#endif
        );
        LCA_mutation_iter++;
    } else if(LCA->children.size()<=1){
        MAT::Mutation temp(dst_add_count_iter->get_position());
        nuc_one_hot par_state=dst_add_count_iter->get_par_state();
        temp.set_par_mut(par_state, par_state);
        temp.set_auxillary(0, (~par_state)&0xf, 0);
        LCA_dst_match(&temp, dst_add_count_iter, LCA, is_src_terminal, is_dst_terminal, LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        , debug, debug_in_iter, begin, end, src_branch, dst_branch, old_score
#endif
        );
    }else {
        parsimony_score_change += is_dst_terminal?dst_add_count_iter->get_default_change_terminal():dst_add_count_iter->get_default_change_internal();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        set_LCA_par_score_change(*dst_add_count_iter, debug, debug_in_iter, begin,end,
                                 parsimony_score_change - old_score);
#endif
    }
}
static void single_child(const MAT::Mutation& LCA_mutation,nuc_one_hot new_child_state,const MAT::Node* LCA,const Mutation_Count_Change& src_count_change,nuc_one_hot& major_allele,int& parsimony_score_change,Mutation_Count_Change_Collection& LCA_parent_mutation_count_change_out){
                nuc_one_hot ori_child_ori_state=LCA_mutation.get_mut_one_hot()|LCA_mutation.get_tie_one_hot();
            assert(LCA->children.size()==1);
            nuc_one_hot ori_child_state=(ori_child_ori_state&(~src_count_change.get_decremented()))|src_count_change.get_incremented();
            major_allele=ori_child_state&new_child_state;
            if(!major_allele){
                parsimony_score_change++;
                major_allele=ori_child_state|new_child_state;
            }
            if(major_allele!=ori_child_ori_state){
                Mutation_Count_Change change(LCA_mutation);
                change.set_change(ori_child_ori_state&(~major_allele),major_allele&(~ori_child_ori_state));
                LCA_parent_mutation_count_change_out.push_back(change);
            }
}
static void LCA_src_match(Mutation_Count_Change_Collection::const_iterator& dst_add_count_iter,Mutation_Count_Change_Collection::const_iterator& dst_add_count_end,const Mutation_Count_Change&src_count_change, nuc_one_hot& major_allele,bool is_dst_terminal,bool is_src_terminal,const MAT::Mutation& LCA_mutation,Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,int &parsimony_score_change,const MAT::Node* LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end,
        const MAT::Node *src_branch, const MAT::Node *dst_branch,int old_score
#endif
) {
    /*if(LCA_mutation.get_position()==28654){
        fputc('a', stderr);
    }*/
    nuc_one_hot dst_added=0;
    if (dst_add_count_iter != dst_add_count_end &&
        dst_add_count_iter->get_position() == src_count_change.get_position()) {
            if (LCA->children.size()<=1) {
                single_child(LCA_mutation, dst_add_count_iter->get_incremented(), LCA, src_count_change, major_allele, parsimony_score_change, LCA_parent_mutation_count_change_out);
            }else{
            major_allele = dbl_inc_dec_mutations(
            *dst_add_count_iter, is_dst_terminal, src_count_change,
            is_src_terminal, LCA_mutation, parsimony_score_change,
            LCA_parent_mutation_count_change_out);
            }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        if (is_dst_terminal) {
            dst_added = dst_add_count_iter->get_incremented();
        }
#endif
        dst_add_count_iter++;
    } else if (is_src_terminal) {
        major_allele = decrement_mutation_count(
            LCA_parent_mutation_count_change_out, LCA_mutation,
            src_count_change, parsimony_score_change, true);
    } else if (is_dst_terminal) {
        if(LCA->children.size()<=1){
            single_child(LCA_mutation, LCA_mutation.get_mut_one_hot(), LCA, src_count_change, major_allele, parsimony_score_change, LCA_parent_mutation_count_change_out);
        }
        else{Mutation_Count_Change temp(LCA_mutation);
        temp.set_change(0, LCA_mutation.get_mut_one_hot());
        major_allele = dbl_inc_dec_mutations(
            temp, true, src_count_change, false, LCA_mutation,
            parsimony_score_change, LCA_parent_mutation_count_change_out);}
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        if (is_dst_terminal) {
            dst_added = LCA_mutation.get_mut_one_hot();
        }
#endif
    } else {
        major_allele = decrement_increment_mutation_count(
            LCA_mutation, src_count_change,
            LCA_parent_mutation_count_change_out, parsimony_score_change);
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    test_allele_count_out_LCA(LCA, src_branch, is_src_terminal, dst_branch,
                              dst_added, LCA_mutation, major_allele,
                              parsimony_score_change - old_score,
                              LCA_parent_mutation_count_change_out,
                              debug, debug_in_iter, end);
#endif
}
void get_LCA_mutation(
    const MAT::Node *LCA, MAT::Node* src_terminal, MAT::Node* dst_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<LCA_merged_states> &dbg_states_in,
    std::vector<MAT::Node *> node_stack, const MAT::Node *src_branch,
    const MAT::Node *dst_branch
#endif
) {
    bool is_src_terminal=src_terminal!=nullptr;
    bool is_dst_terminal=dst_terminal!=nullptr;
    MAT::Mutations_Collection::const_iterator LCA_mutation_iter =
        LCA->mutations.begin();
    MAT::Mutations_Collection::const_iterator LCA_mutation_end =
        LCA->mutations.end();
    auto dst_add_count_iter = from_dst_add.begin();
    auto dst_add_count_end = from_dst_add.end();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    auto debug_in_iter = dbg_states_in.begin();
    auto end = dbg_states_in.end();
    auto debug_begin = dbg_states_in.begin();
#endif
    assert(!(is_src_terminal & is_dst_terminal));
    int mut_size =
        is_src_terminal ? src_terminal->mutations.size() : from_src_remove.size();
    for (size_t idx = 0; idx < mut_size; idx++) {
        Mutation_Count_Change src_count_change =
            is_src_terminal ? src_terminal->mutations[idx] : from_src_remove[idx];
        if (is_src_terminal) {
            const auto& src_mut=src_terminal->mutations[idx];
            src_count_change.set_change(
                src_mut.get_mut_one_hot()|src_mut.get_tie_one_hot(), 0);
        }
        /*if (src_count_change.get_position()==241) {
            fputc('a', stderr);
        }*/
        while ((dst_add_count_iter != dst_add_count_end )&&
               dst_add_count_iter->get_position() < src_count_change.get_position()) {
            LCA_dst(LCA_mutation_iter, LCA_mutation_end, dst_add_count_iter,
                    LCA, is_src_terminal, is_dst_terminal,
                    LCA_parent_mutation_count_change_out,
                    parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    , debug, debug_in_iter,debug_begin ,end,
                    src_branch, dst_branch
#endif
                );
            dst_add_count_iter++;
        }
        while (LCA_mutation_iter != LCA_mutation_end&&LCA_mutation_iter->get_position() <
                   src_count_change.get_position()
                   ) {
            LCA_no_match(
                LCA, is_src_terminal, is_dst_terminal, LCA_mutation_iter,
                LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,debug, debug_in_iter, end, src_branch, dst_branch
            #endif
            );
        }
        nuc_one_hot major_allele;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        int old_score = parsimony_score_change;
#endif
        /*if (src_count_change.get_position() ==241) {
            fputc('a', stderr);
        }*/
        if (LCA_mutation_iter != LCA_mutation_end &&
            LCA_mutation_iter->get_position() ==
                src_count_change.get_position()) {
                    LCA_src_match(dst_add_count_iter, dst_add_count_end, src_count_change, major_allele, is_dst_terminal, is_src_terminal, *LCA_mutation_iter, LCA_parent_mutation_count_change_out, parsimony_score_change, LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    , debug, debug_in_iter, debug_begin, end, src_branch, dst_branch, old_score
#endif
                    );
        LCA_mutation_iter++;
        }else if (LCA->children.size()<=1) {
        MAT::Mutation temp(src_count_change.get_position());
        nuc_one_hot par_state=src_count_change.get_par_state();
        temp.set_par_mut(par_state, par_state);
        temp.set_auxillary(0, (~par_state)&0xf, 0);
        LCA_src_match(dst_add_count_iter, dst_add_count_end, src_count_change, major_allele, is_dst_terminal, is_src_terminal, temp, LCA_parent_mutation_count_change_out, parsimony_score_change, LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        , debug, debug_in_iter, debug_begin, end, src_branch, dst_branch, old_score
#endif
        );
        }
        else {
                parsimony_score_change +=is_src_terminal? src_count_change.get_default_change_terminal():src_count_change.get_default_change_internal();
                if (dst_add_count_iter!=dst_add_count_end&&dst_add_count_iter->get_position() ==
                src_count_change.get_position()) {
                parsimony_score_change +=is_dst_terminal?dst_add_count_iter->get_default_change_terminal():
                    dst_add_count_iter->get_default_change_internal();
                dst_add_count_iter++;
                }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            set_LCA_par_score_change(src_count_change, debug, debug_in_iter,debug_begin,
                                     end, parsimony_score_change - old_score);
#endif
        }

    }

    while (dst_add_count_iter != dst_add_count_end) {
        LCA_dst(LCA_mutation_iter, LCA_mutation_end, dst_add_count_iter, LCA,
                src_terminal, dst_terminal,
                LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,
                debug, debug_in_iter,debug_begin, end, src_branch, dst_branch
#endif
                );
        dst_add_count_iter++;
    }

    while (LCA_mutation_iter != LCA_mutation_end) {
        LCA_no_match(LCA, src_terminal, dst_terminal, LCA_mutation_iter,
                     LCA_parent_mutation_count_change_out,
                     parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                     , debug, debug_in_iter, end,
                     src_branch, dst_branch
#endif
                     );
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    while (debug_in_iter < end) {
        debug.emplace_back(debug_in_iter->position,debug_in_iter->par_allele,debug_in_iter->par_allele,0,Mutation_Count_Change());
        debug_in_iter++;
    }
#endif
}
