#include "process_individual_mutation.hpp"
#include "process_each_node.hpp"
#include <cstdio>

static int get_major_allele_three(nuc_one_hot first_child,nuc_one_hot second_child,nuc_one_hot third_child,nuc_one_hot& major_allele_out){
    major_allele_out=first_child&second_child&third_child;
    if (major_allele_out) {
        return 0;
    }
    uint32_t second_third=(second_child<<4)|third_child;
    uint32_t first_second_third=(first_child<<8)|second_third;
    uint32_t second_third_first=(second_third<<4)|first_child;
    uint32_t any_two=first_second_third|second_third_first;
    if (any_two) {
        major_allele_out=0xf&(any_two|(any_two>>4)|(any_two>>8));
        return 1;
    }
    major_allele_out=first_child|second_child|third_child;
    return 2;
}

static void LCA_no_match(
    const MAT::Node *LCA,const MAT::Node* src_branch,bool is_src_terminal,const MAT::Node* dst_branch,bool is_dst_terminal,
    MAT::Mutations_Collection::const_iterator &LCA_mutation_iter,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
    const std::vector<LCA_merged_states>::const_iterator &end
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    int old_score = parsimony_score_change;
#endif
/*if (LCA_mutation_iter->get_position()==241) {
    fputc('a', stderr);
}*/
    nuc_one_hot major_allele=LCA_mutation_iter->get_all_major_allele();
    if (is_dst_terminal) {
        major_allele=LCA_mutation_iter->get_mut_one_hot();
        register_change_from_new_state(LCA_parent_mutation_count_change_out, 0,*LCA_mutation_iter, major_allele);
    } else if (is_src_terminal) {
        major_allele = src_branch==LCA->children[0]?LCA_mutation_iter->get_right_child_state():LCA_mutation_iter->get_left_child_state();
        parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, 0,*LCA_mutation_iter, major_allele);
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    test_allele_count_out_LCA(
        LCA, src_branch, is_src_terminal, dst_branch,
        is_dst_terminal? (uint8_t)LCA_mutation_iter->get_mut_one_hot()
                        : (uint8_t)0,
        *LCA_mutation_iter, major_allele, parsimony_score_change - old_score,
        LCA_parent_mutation_count_change_out, debug, debug_in_iter, end);
#endif
    LCA_mutation_iter++;
}

static void LCA_dst_match(const MAT::Mutation* LCA_mutation_iter,
        const Mutation_Count_Change_Collection::const_iterator &dst_add_count_iter,
        const MAT::Node *LCA,const MAT::Node* src_branch, bool is_src_terminal,const MAT::Node* dst_branch, bool is_dst_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end,int old_score
#endif
){
        nuc_one_hot major_allele;
        if (is_dst_terminal) {
            major_allele = increment_mutation_count(
                LCA_parent_mutation_count_change_out, *LCA_mutation_iter,
                *dst_add_count_iter, parsimony_score_change);
                parsimony_score_change++;
        } else if (is_src_terminal) {
            major_allele=dst_add_count_iter->get_new_state();
            parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, 0,*LCA_mutation_iter, major_allele);
        } else {
            int new_par_score=get_new_major_allele_binary_node(LCA->children[0]==dst_branch?LCA_mutation_iter->get_right_child_state():LCA_mutation_iter->get_left_child_state(), dst_add_count_iter->get_new_state(), major_allele);
            parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, new_par_score, *LCA_mutation_iter, major_allele);
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
        const MAT::Node *LCA, const MAT::Node* src_branch,bool is_src_terminal,const MAT::Node* dst_branch, bool is_dst_terminal,
        Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
        int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end
#endif
) {
    while (LCA_mutation_iter != LCA_mutation_end&&LCA_mutation_iter->get_position() <
               dst_add_count_iter->get_position()
           ) {
        LCA_no_match(LCA,src_branch, is_src_terminal,dst_branch, is_dst_terminal,LCA_mutation_iter,
                     LCA_parent_mutation_count_change_out,
                     parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                     , debug, debug_in_iter, end
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
        LCA_dst_match(&(*LCA_mutation_iter), dst_add_count_iter, LCA, src_branch,is_src_terminal, dst_branch,is_dst_terminal, LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        , debug, debug_in_iter, begin, end, old_score
#endif
        );
        LCA_mutation_iter++;
    }else {
        parsimony_score_change += is_dst_terminal?dst_add_count_iter->get_default_change_terminal():dst_add_count_iter->get_default_change_internal();
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        set_LCA_par_score_change(*dst_add_count_iter, debug, debug_in_iter, begin,end,
                                 parsimony_score_change - old_score);
#endif
    }
}

static void LCA_src_match(Mutation_Count_Change_Collection::const_iterator& dst_add_count_iter,Mutation_Count_Change_Collection::const_iterator& dst_add_count_end,const Mutation_Count_Change&src_count_change, nuc_one_hot& major_allele,bool is_dst_terminal,bool is_src_terminal,const MAT::Mutation& LCA_mutation,Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,int &parsimony_score_change,const MAT::Node* LCA,const MAT::Node* src_branch, const MAT::Node* dst_branch
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        ,
        std::vector<state_change_hist_dbg> &debug,
        std::vector<LCA_merged_states>::const_iterator &debug_in_iter,
        const std::vector<LCA_merged_states>::const_iterator &begin,
        const std::vector<LCA_merged_states>::const_iterator &end,
        int old_score
#endif
) {
    /*if(LCA_mutation.get_position()==28654){
        fputc('a', stderr);
    }*/
    nuc_one_hot dst_added=0;
    if (dst_add_count_iter != dst_add_count_end &&
        dst_add_count_iter->get_position() == src_count_change.get_position()) {
        if (is_src_terminal) {
            major_allele=dst_add_count_iter->get_new_state();
            parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, 0, LCA_mutation, major_allele);
        }
        else if (is_dst_terminal) {
            int new_mut_count=get_major_allele_three(src_branch==LCA->children[0]?LCA_mutation.get_right_child_state():LCA_mutation.get_left_child_state(), src_count_change.get_new_state(), dst_add_count_iter->get_new_state(),major_allele);
            parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, new_mut_count, LCA_mutation, major_allele);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            dst_added = dst_add_count_iter->get_incremented();
#endif
        }else {
            int new_mut_count=get_new_major_allele_binary_node(src_count_change.get_new_state(), dst_add_count_iter->get_new_state(), major_allele);
            parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, new_mut_count, LCA_mutation, major_allele);
        }
        dst_add_count_iter++;
    } else if (is_src_terminal) {
        major_allele=src_branch==LCA->children[0]?LCA_mutation.get_right_child_state():LCA_mutation.get_left_child_state();
        parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, 0, LCA_mutation, major_allele);
    } else if (is_dst_terminal) {
        int new_mut_count=get_major_allele_three(LCA_mutation.get_mut_one_hot(), src_count_change.get_new_state(), src_branch==LCA->children[0]?LCA_mutation.get_right_child_state():LCA_mutation.get_left_child_state(), major_allele);
        parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, new_mut_count, LCA_mutation, major_allele);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        dst_added = LCA_mutation.get_mut_one_hot();
#endif
    } else {
        int new_mut_count=get_new_major_allele_binary_node(src_count_change.get_new_state(), src_branch==LCA->children[0]?LCA_mutation.get_right_child_state():LCA_mutation.get_left_child_state(),major_allele);
        parsimony_score_change+=register_change_from_new_state(LCA_parent_mutation_count_change_out, new_mut_count, LCA_mutation, major_allele);
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    test_allele_count_out_LCA(LCA, src_branch, is_src_terminal, dst_branch,
                              dst_added, LCA_mutation, major_allele,
                              parsimony_score_change - old_score,
                              LCA_parent_mutation_count_change_out,
                              debug, debug_in_iter, end);
#endif
}
void get_LCA_mutation_binary_node(
    const MAT::Node *LCA, const MAT::Node *src_branch,bool is_src_terminal,
    const MAT::Node *dst_branch, bool is_dst_terminal,
    const Mutation_Count_Change_Collection &from_src_remove,
    const Mutation_Count_Change_Collection &from_dst_add,
    Mutation_Count_Change_Collection &LCA_parent_mutation_count_change_out,
    int &parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<LCA_merged_states> &dbg_states_in,
    std::vector<MAT::Node *> node_stack
#endif
) {
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
    int mut_size = from_src_remove.size();
    for (size_t idx = 0; idx < mut_size; idx++) {
        Mutation_Count_Change src_count_change = from_src_remove[idx];
        while ((dst_add_count_iter != dst_add_count_end )&&
               dst_add_count_iter->get_position() < src_count_change.get_position()) {
            LCA_dst(LCA_mutation_iter, LCA_mutation_end, dst_add_count_iter,
                    LCA, src_branch,is_src_terminal,dst_branch,is_dst_terminal,
                    LCA_parent_mutation_count_change_out,
                    parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    , debug, debug_in_iter,debug_begin ,end
#endif
                );
            dst_add_count_iter++;
        }
        while (LCA_mutation_iter != LCA_mutation_end&&LCA_mutation_iter->get_position() <
                   src_count_change.get_position()
                   ) {
            LCA_no_match(
                LCA, src_branch,is_src_terminal,dst_branch,is_dst_terminal, LCA_mutation_iter,
                LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,debug, debug_in_iter, end
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
                    LCA_src_match(dst_add_count_iter, dst_add_count_end, src_count_change, major_allele, is_dst_terminal, is_src_terminal, *LCA_mutation_iter, LCA_parent_mutation_count_change_out, parsimony_score_change, LCA,src_branch,dst_branch
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    , debug, debug_in_iter, debug_begin, end, old_score
#endif
                    );
        LCA_mutation_iter++;
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
                src_branch,is_src_terminal ,dst_branch,is_dst_terminal,
                LCA_parent_mutation_count_change_out, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,
                debug, debug_in_iter,debug_begin, end
#endif
                );
        dst_add_count_iter++;
    }

    while (LCA_mutation_iter != LCA_mutation_end) {
        LCA_no_match(LCA,src_branch, is_src_terminal,dst_branch,is_dst_terminal, LCA_mutation_iter,
                     LCA_parent_mutation_count_change_out,
                     parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                     , debug, debug_in_iter, end
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
