#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "conflicts.hpp"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <list>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_set>
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

struct Merge_Comparator{
    bool operator()(const Merge_Discriptor& first, const Merge_Discriptor& second)const {
        return first.shared_mutations.size()>second.shared_mutations.size();
    }
};
static void check_mergable(MAT::Node* parent,MAT::Mutations_Collection& mutations,std::vector<Merge_Discriptor> possible_merges){
    Merge_Discriptor merge;
    possible_merges.reserve(parent->children.size());
    for(MAT::Node* c:parent->children){
        c->mutations.set_difference(mutations, merge.to_merge_width_unique_mutations, merge.src_unique_mutations, merge.shared_mutations);
        if (merge.shared_mutations.size()==0) {
            possible_merges.push_back(merge);
        }
        merge.src_unique_mutations.clear();
        merge.to_merge_width_unique_mutations.clear();
        merge.shared_mutations.clear();
    }
    possible_merges.shrink_to_fit();
    if(possible_merges.size()>=1){
        std::sort(possible_merges.begin(),possible_merges.end(),Merge_Comparator());
    }
}

size_t remove_conflicts(MAT::Node* src,std::vector<Move> possible_moves,std::vector<bool>& is_conflict){
    size_t non_conflicting=possible_moves.size();
    Move* last_nonconflicting=nullptr;
    for (Move& m:possible_moves){
        if(check_mut_conflict(src, m.states)
        ||check_mut_conflict(m.dst, m.states)
        ||check_loop_conflict(src, m.dst)){
            non_conflicting--;
            is_conflict.push_back(true);
        }else {
            is_conflict.push_back(false);
            last_nonconflicting=&m;
        }
    }
    if (non_conflicting==1) {
        register_loop_conflict(src, last_nonconflicting->dst);
        register_mut_conflict(src, last_nonconflicting->states);
        register_mut_conflict(last_nonconflicting->dst, last_nonconflicting->states);
    }
    return non_conflicting;
}
struct check_individual_move_profitable {
    tbb::mutex& mutex;
    std::vector<Dst_Mut>& targets;
    Profitable_Moves* result;
    std::unordered_map<int, Fitch_Sankoff_Result>* src_tip_fs_result;
    std::vector<size_t>& good_idx;
    void operator()( const tbb::blocked_range<size_t>& r )const{
        for( size_t i=r.begin(); i!=r.end(); ++i ){
            MAT::Mutations_Collection new_tip_mutation;
            Dst_Mut& target=targets[i];
            auto parsimony_score_change=get_parsimony_score_change(src_tip_fs_result,target,new_tip_mutation);
            if (parsimony_score_change>=0||parsimony_score_change>result->score_change) {
                return;
            }
            
            std::vector<Merge_Discriptor> possible_merges;
            check_mergable(target.dst, new_tip_mutation, possible_merges);
            if(parsimony_score_change>result->score_change) return;
            {
                std::lock_guard<tbb::mutex> lock(mutex);
                if(parsimony_score_change>result->score_change) return;
                
                if(parsimony_score_change<result->score_change){
                    good_idx.clear();
                    result->moves.clear();
                    result->score_change=parsimony_score_change;
                } 
                result->moves.emplace_back();
                Move& out=result->moves.back();
                out.dst=target.dst;
                out.new_tip_mutations.swap(new_tip_mutation);
                out.merger.swap(possible_merges);
                good_idx.emplace_back(i);
            }
        } 
    }
};


static void fill_result(Possible_Moves* in,Profitable_Moves* result,std::vector<size_t>& good_idx){
    result->src_tip_fs_result=in->src_tip_fs_result;
    result->shared=in->shared;
    result->src=in->src;
    auto iter=result->moves.begin();
        for(auto idx:good_idx){
            iter->states.reserve(in->dst[idx].mut.size());
            for(auto& m:in->dst[idx].mut){
                iter->states.push_back(&((*in->src_tip_fs_result)[m.position]));
            }
        }
}
Profitable_Moves* Profitable_Moves_Enumerator::operator() (Possible_Moves* in)const{
    if(!in){
        return nullptr;
    }
    tbb::mutex mutex;
    Profitable_Moves* result=new Profitable_Moves;
    std::vector<size_t> good_idx;
    result->score_change=0;
    good_idx.reserve(2);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,in->dst.size()),check_individual_move_profitable({mutex,in->dst,result,in->src_tip_fs_result,good_idx}));
    if(result->score_change<0){
        fill_result(in,result,good_idx);
        size_t non_conflicting;
        std::vector<bool> is_conflict;
        is_conflict.reserve(result->moves.size());
        {
            std::lock_guard<tbb::mutex> lock(conflict_guard);
            non_conflicting=remove_conflicts(in->src, result->moves, is_conflict);
        }
        if (non_conflicting) {
            if(non_conflicting!=result->moves.size()){
                std::vector<Move> filtered;
                for (size_t i=0; i<is_conflict.size(); i++) {
                    if(!is_conflict[i]){
                        filtered.push_back(result->moves[i]);
                    }
                }
                result->moves.swap(filtered);
            }
            delete in;
            return result;
        }
        postponed.push_back(in->src);
    }
    
    delete in;
    delete in->src_tip_fs_result;
    delete result;
    return nullptr;
}