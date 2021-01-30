#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

static void distribute(MAT::Mutations_Collection& src, std::vector<std::vector<char>>& dst){
    assert(src.size()==dst.size());
    auto iter=src.begin();
    for(auto& a:dst){
        a.push_back(iter->mut_nuc);
        iter++;
    }
} 


static std::vector<std::vector<char>>  get_original_states(std::vector<MAT::Node*> dfs_ordered_nodes, const std::pair<size_t, size_t> &range,MAT::Mutations_Collection& target){
    MAT::Node* start_node=dfs_ordered_nodes[range.first];
    std::vector<std::vector<char>> result(target.size());
    
    for(std::vector<char>& a:result){
        a.reserve(range.second-range.first);    
    }
    
    for(Mutation_Annotated_Tree::Mutation& m : target){
        m.mut_nuc=get_genotype(start_node, m);
    }
    distribute(target, result);

    for (size_t idx=range.first+1; idx<range.second; idx++) {
        size_t loci_idx=0;
        auto par_idx=dfs_ordered_nodes[idx]->parent->index-range.first;
        for(Mutation_Annotated_Tree::Mutation& m : target){
            m.mut_nuc=result[loci_idx][par_idx];
            loci_idx++;
        }
        dfs_ordered_nodes[idx]->mutations.batch_find(target);
        #ifndef NDEBUG
        for(Mutation_Annotated_Tree::Mutation& m : target){
            assert(m.mut_nuc==get_genotype(dfs_ordered_nodes[idx], m));
        }
        #endif
        distribute(target, result);
    }
    return result;
}

Possible_Moves* Parsimony_Score_Calculator::operator()(Possible_Moves* in)const{
    if(!in){
        return nullptr;
    }
    
    in->shared.range=Fitch_Sankoff::dfs_range(in->src, dfs_ordered_nodes);
    std::vector<std::vector<char>> original_states=get_original_states(dfs_ordered_nodes, in->shared.range, *(in->to_search));
    std::unordered_map<int, Fitch_Sankoff_Result>* src_tip_fs_result=new std::unordered_map<int, Fitch_Sankoff_Result>;
    std::vector<Fitch_Sankoff_Result*> results;
    results.reserve(in->to_search->size());
    src_tip_fs_result->reserve(2*in->to_search->size());
    
    MAT::Mutations_Collection& mutations=*(in->to_search);
    for(MAT::Mutation& m: mutations){
        auto iter=src_tip_fs_result->emplace(m.position,Fitch_Sankoff_Result());
        results.push_back(&(iter.first->second));
    }
    tbb::parallel_for(tbb::blocked_range<size_t>(0,in->to_search->size()),
    [&](tbb::blocked_range<size_t> r) {
        for(size_t j=r.begin(); j<r.end(); ++j){
            auto this_result=results[j];
            this_result->original_state=original_states[j];
            Fitch_Sankoff::Scores_Type scores;
            Fitch_Sankoff::States_Type& states=this_result->states;
            this_result->original_tip_score=Fitch_Sankoff::sankoff_backward_pass(in->shared.range,mutations[j],dfs_ordered_nodes,scores,states,original_states[j]);
            this_result->tip_score=scores.back();
        }
    });
    return in;
}