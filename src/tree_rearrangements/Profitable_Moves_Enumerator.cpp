#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <iterator>
#include <list>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_set>
#include <utility>
#include <vector>

static void distribute(MAT::Mutations_Collection& src, std::vector<std::vector<char>>& dst){
    assert(src.size()==dst.size());
    auto iter=src.begin();
    for(auto& a:dst){
        a.push_back(iter->mut_nuc);
        iter++;
    }
} 


static std::vector<std::vector<char>>  get_original_states(std::vector<MAT::Node*> dfs_ordered_nodes, const std::pair<size_t, size_t> &range,MAT::Mutations_Collection target){
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
MAT::Node* get_mutation_path(MAT::Mutations_Collection& mutations,MAT::Node* src, MAT::Node* dst){
    std::unordered_set<MAT::Node*> src_to_root;
    std::vector<MAT::Node*> src_to_root_path;
    assert(src!=dst);
    MAT::Node* LCA=src->parent;
    while (LCA) {
        src_to_root.insert(LCA);
        src_to_root_path.push_back(LCA);
        LCA=LCA->parent;
    }

    LCA=dst;
    std::vector<MAT::Node*> dst_to_root_path;
    while (!src_to_root.count(LCA)) {
        dst_to_root_path.push_back(LCA);
        LCA=LCA->parent;
    }
    assert(dst_to_root_path.back()->parent==LCA);
    
    mutations=src->mutations;
    for(MAT::Node* n:src_to_root_path){
        mutations.merge(n->mutations, MAT::Mutations_Collection::INVERT_MERGE);
        if (n==LCA) {
            break;
        }
    }

    for (auto iter=dst_to_root_path.rbegin(); iter<dst_to_root_path.rend(); iter++) {
        mutations.merge((*iter)->mutations, MAT::Mutations_Collection::MERGE);
    }

    return LCA;
}

static void patch_sankoff_result(size_t start_idx, Fitch_Sankoff_Result* out, MAT::Node* src, MAT::Node* dst){
    Fitch_Sankoff::set_internal_score(*dst, out->scores, start_idx, out->states, src);
    MAT::Node* changing_node=dst->parent;
    while (changing_node->index>=start_idx) {
        Fitch_Sankoff::set_internal_score(*changing_node, out->scores, start_idx, out->states);
    }
}
static int calculate_parsimony_score_change(std::pair<size_t, size_t>& range, Fitch_Sankoff_Result* out, MAT::Node* src, MAT::Node* dst,char par_nuc){
    int start_idx=range.second-1;
    int original_parsimony_score=out->scores.back()[par_nuc];
    patch_sankoff_result(start_idx, out, src, src->parent);
    patch_sankoff_result(start_idx, out, src, dst);
    return out->scores.back()[par_nuc]-original_parsimony_score;
}

void Profitable_Moves_Enumerator::operator() (Possible_Move* in)const{
    MAT::Mutations_Collection mutations;
    MAT::Node* LCA=get_mutation_path(mutations, in->src, in->dst);
    std::pair<size_t, size_t> range=Fitch_Sankoff::dfs_range(LCA,dfs_ordered_nodes);
    std::vector<std::vector<char>> original_states=get_original_states(dfs_ordered_nodes, range, mutations);
    std::atomic_int score_changes(0);
    std::vector<Fitch_Sankoff_Result*> unchanged_states(mutations.size());
    std::vector<Fitch_Sankoff_Result*> moved_states(mutations.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,mutations.size()),[&](tbb::blocked_range<size_t> r){
        for (size_t i=r.begin();i<r.end() ; i++) {
            Fitch_Sankoff::sankoff_backward_pass(range, mutations[i], dfs_ordered_nodes, unchanged_states[i]->scores, unchanged_states[i]->states, original_states[i]);
            Fitch_Sankoff_Result* patched=new Fitch_Sankoff_Result(*unchanged_states[i]);
            int this_change=calculate_parsimony_score_change(range,patched,in->src,in->dst,get_genotype(LCA->parent, mutations[i]));
            moved_states[i]=patched;
            score_changes.fetch_add(this_change,std::memory_order_relaxed);
        }
    });
    if (score_changes<0) {
        Move* out=new Move;
        out->score_change=score_changes;
        out->src=in->src;
        out->dst=in->dst;
        for (size_t i=0; i<mutations.size(); i++) {
            moved_states[i]->mutation=mutations[i];
            moved_states[i]->original_state.swap(original_states[i]);
        }
        out->states.swap(moved_states);
        profitable_moves.push_back(out);
    }else{
        for(auto c:moved_states){
            delete c;
        }
    }
    for(auto c:unchanged_states){
        delete c;
    }
    delete in;
}