#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <vector>

static int get_parsimony_score_change( std::unordered_map<int, Fitch_Sankoff_Result>* src_tip_fs_result,Dst_Mut& target,MAT::Mutations_Collection& tip_mutation){
    int score_change=0;
    tip_mutation.reserve(target.mut.size());
    for(MAT::Mutation& m: target.mut){
        Fitch_Sankoff_Result& result=(*src_tip_fs_result)[m.position];
        char par_state=one_hot_to_two_bit(m.mut_nuc);
        std::pair<int, char> if_placed=Fitch_Sankoff::get_child_score_on_par_nuc(par_state,result.tip_score);
        result.states.back()=if_placed.second;
        score_change+=(result.original_tip_score-if_placed.first);
        if (if_placed.second!=par_state) {
            MAT::Mutation temp(m);
            temp.par_nuc=temp.mut_nuc;
            temp.mut_nuc=two_bit_to_one_hot(if_placed.second);
            tip_mutation.push_back(m);
        }
    }
    return score_change;
}

static void check_mergable(MAT::Node* parent,MAT::Mutations_Collection& mutations,std::list<Merge_Discriptor*> possible_merges){
    for(MAT::Node* c:parent->children){
        
    }
}

struct check_individual_move_profitable {
    tbb::mutex& mutex;
    std::vector<Dst_Mut>& targets;
    Profitable_Moves* result;
    std::unordered_map<int, Fitch_Sankoff_Result>* src_tip_fs_result;
    void operator()( const tbb::blocked_range<size_t>& r )const{
        for( size_t i=r.begin(); i!=r.end(); ++i ){
            MAT::Mutations_Collection new_tip_mutation;
            Dst_Mut& target=targets[i];
            auto parsimony_score_change=get_parsimony_score_change(src_tip_fs_result,target,new_tip_mutation);
            if (parsimony_score_change>=0||parsimony_score_change>result->score_change) {
                break;
            }
            
        } 
    }
};
Profitable_Moves* Profitable_Moves_Enumerator::operator() (Possible_Moves* in)const{
    if(!in){
        return nullptr;
    }
    
    delete in->src_tip_fs_result;
    delete in;
}