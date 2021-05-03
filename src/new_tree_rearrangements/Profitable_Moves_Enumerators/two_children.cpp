#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
int get_new_major_allele_binary_node(nuc_one_hot left_child,nuc_one_hot right_child,nuc_one_hot& major_allele_out){
    major_allele_out=left_child&right_child;
    if (!major_allele_out) {
        major_allele_out=left_child|right_child;
        return 1;
    }
    return 0;
}

int register_change_from_new_state(Mutation_Count_Change_Collection &out,int new_mut_count,
                            const MAT::Mutation &pos,nuc_one_hot new_state) {
if(new_state!=pos.get_all_major_allele()){
    out.emplace_back(pos);
    nuc_one_hot incremented_allele=new_state&(~pos.get_all_major_allele());
    nuc_one_hot decremented_allele=pos.get_all_major_allele()&(~new_state);
    Mutation_Count_Change &last_added = out.back();
    last_added.set_change(decremented_allele, incremented_allele,new_state);
    }
    int old_mutation_count=!(pos.get_left_child_state()&pos.get_right_child_state());
    return new_mut_count-old_mutation_count;
}

void get_two_child_intermediate_node_mutations(
    const MAT::Node *this_node,
    const MAT::Node *changed_child,
    const Mutation_Count_Change_Collection &this_node_mutation_count_change,
    Mutation_Count_Change_Collection &parent_node_mutation_count_change,
    int &parent_parsimony_score_change
// const New_Tie_Collection_t &parent_ties_in,
// New_Tie_Collection_t &decrease_if_match_parent_state_out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg>& debug,
    const std::vector<MAT::Node *>& node_stack
#endif
){
    bool changed_first_children;
    if (changed_child==this_node->children[0]) {
        changed_first_children=true;
    }else {
        changed_first_children=false;
        assert(changed_child==this_node->children[1]);
    }
    auto mutation_iter=this_node->mutations.begin();
    auto mutation_end=this_node->mutations.end();
    #ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    auto debug_iter = debug.begin();
    auto debug_end=debug.end();
    #endif
    for(const auto& this_mut:this_node_mutation_count_change){
        consume_parent_mutations(
            mutation_iter, mutation_end, this_mut
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_iter,debug_end
#endif
            );
        if (mutation_iter!=mutation_end&&mutation_iter->get_position()==this_mut.get_position()) {
            nuc_one_hot major_allele;
            nuc_one_hot left_child_state;
            nuc_one_hot right_child_state;
            if (changed_first_children) {
                left_child_state=this_mut.get_new_state();
                right_child_state=mutation_iter->get_right_child_state();
            }else {
                left_child_state=mutation_iter->get_left_child_state();
                right_child_state=this_mut.get_new_state();
            }
            int new_count=get_new_major_allele_binary_node(left_child_state, right_child_state, major_allele);
            int mutation_count_change=register_change_from_new_state(parent_node_mutation_count_change, new_count, *mutation_iter, major_allele);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_count_out(this_node, *mutation_iter, major_allele, 0,
                                  parent_node_mutation_count_change,
                                  mutation_count_change,
                                  debug_iter,debug_end, node_stack);
#endif
            parent_parsimony_score_change+=mutation_count_change;
            mutation_iter++;
        }else {
            int mutation_count_change= this_mut.get_default_change_internal();
            parent_parsimony_score_change+=mutation_count_change;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            update_dbg_vector_score_only(this_mut.get_position(), debug_iter,debug_end,
                                         mutation_count_change);
#endif
        }
    }
}

void get_child_removed_binary_node(
    Mutation_Count_Change_Collection &parent_mutation_count_change_out,
    //,New_Tie_Collection_t &decrease_if_match_parent_state,
    const MAT::Node *src, int &parent_parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    std::vector<state_change_hist_dbg> &debug,
    const std::vector<MAT::Node *> &node_stack
#endif
){
    MAT::Node* parent_node=src->parent;
    bool removed_left_child=parent_node->children[0]==src;
    assert(removed_left_child||parent_node->children[1]==src);
    for(const auto& mutation:parent_node->mutations){
        nuc_one_hot new_major_allele=removed_left_child?mutation.get_right_child_state():mutation.get_left_child_state();
        int score_change=register_change_from_new_state(parent_mutation_count_change_out, 0, mutation, new_major_allele);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            test_allele_out_init(
                src->parent, src, mutation,
                score_change, new_major_allele,
                0, parent_mutation_count_change_out, debug);
#endif
        parent_parsimony_score_change+=score_change;
    }
}
