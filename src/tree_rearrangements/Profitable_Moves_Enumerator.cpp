#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <ios>
#include <iterator>
#include <list>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static void distribute(MAT::Mutations_Collection& src, std::vector<Fitch_Sankoff::States_Type>& dst,MAT::Node* this_node){
    assert(src.size()==dst.size());
    for(size_t idx=0;idx<src.size();idx++){
        #ifndef NDEBUG
        dst[idx].emplace_back(src[idx].mut_nuc,this_node);
        #else
        dst[idx].push_back(src[idx].mut_nuc);
        #endif
        assert(src[idx].mut_nuc==get_genotype(this_node,src[idx]));
    }
} 

static void remove_child(MAT::Node* child_to_remove){
        auto parent=child_to_remove->parent;
        auto iter=std::find(parent->children.begin(),parent->children.end(),child_to_remove);
        assert(iter!=parent->children.end());
        parent->children.erase(iter);
}


static std::vector<char>  get_original_states(std::vector<MAT::Node*> dfs_ordered_nodes, const std::pair<size_t, size_t> &range,MAT::Mutations_Collection target,std::vector<Fitch_Sankoff::States_Type>& result){
    std::vector<char> start_node_parent_states(target.size());
    MAT::Node* start_node=dfs_ordered_nodes[range.first];
    result.reserve(target.size());
    for(size_t i=0;i<target.size();i++){
        Mutation_Annotated_Tree::Mutation &m = target[i];
        char parent_state=(start_node->parent)?get_genotype(start_node->parent, m):m.ref_nuc;
        m.mut_nuc=parent_state;
        start_node_parent_states[i]=parent_state;
        auto mut_iter=start_node->mutations.find(m.position);
        if(mut_iter!=start_node->mutations.end()){
            m.mut_nuc=mut_iter->mut_nuc;
        }
        result.emplace_back();
        result.back().reserve(range.second-range.first);    
    }
    distribute(target, result,start_node);

    for (size_t idx=range.first+1; idx<range.second; idx++) {
        size_t loci_idx=0;
        auto par_idx=dfs_ordered_nodes[idx]->parent->index-range.first;
        for(Mutation_Annotated_Tree::Mutation& m : target){
            m.mut_nuc=result[loci_idx][par_idx];
            loci_idx++;
        }
        dfs_ordered_nodes[idx]->mutations.batch_find(target);
        distribute(target, result,dfs_ordered_nodes[range.first+result[0].size()]);
    }
    return start_node_parent_states;
}

MAT::Node* get_mutation_path(MAT::Mutations_Collection& mutations,MAT::Node* src, MAT::Node* dst,std::vector<MAT::Node*>& path){
    std::unordered_set<MAT::Node*> src_to_root;
    std::vector<MAT::Node*> src_to_root_path;
    assert(src!=dst);
    MAT::Node* LCA=src;
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
    assert(dst_to_root_path.empty()||dst_to_root_path.back()->parent==LCA);

    mutations.clear();
    for(MAT::Node* n:src_to_root_path){
        if (n==LCA) {
            break;
        }
        path.push_back(n);
        mutations.merge(n->mutations, MAT::Mutations_Collection::KEEP_SELF);
    }
    path.push_back(LCA);
    for (auto iter=dst_to_root_path.rbegin(); iter<dst_to_root_path.rend(); iter++) {
        path.push_back(*iter);
        mutations.merge((*iter)->mutations, MAT::Mutations_Collection::KEEP_SELF);
    }

    return LCA;
}

static void patch_sankoff_result(size_t offset, Fitch_Sankoff_Result* out, MAT::Node* src, MAT::Node* dst,size_t end_idx){
    if(dst->is_leaf()){
        assert(out->scores[offset-dst->index].node==dst);
        Fitch_Sankoff::Score_Type leaf_score(out->scores[offset-dst->index]);
        Fitch_Sankoff::set_internal_score(*dst, out->scores, offset, src,&leaf_score);
    }else{
        Fitch_Sankoff::set_internal_score(*dst, out->scores, offset, src);
    }
    MAT::Node* changing_node=dst->parent;
    while (changing_node&&changing_node->index>end_idx) {
        Fitch_Sankoff::set_internal_score(*changing_node, out->scores, offset);
        changing_node=changing_node->parent;
    }
}
static int calculate_parsimony_score_change(std::pair<size_t, size_t>& range, Fitch_Sankoff_Result* out, MAT::Node* src, MAT::Node* dst,char par_nuc, MAT::Node* LCA){
    assert(LCA->index==range.first);
    int offset=range.second-1;
    char par_nuc_idx=one_hot_to_two_bit(par_nuc);
    int original_parsimony_score=out->scores.back()[par_nuc_idx];
    size_t end_idx=range.first;
    patch_sankoff_result(offset, out, src, src->parent,end_idx);
    patch_sankoff_result(offset, out, src, dst,end_idx);
    //Deal with LCA separately, as it is influenced by both src and dst, especially the case when src is child of LCA, then src need to be processed, or dst is LCA itself (move to parent of parent), dst is to be processed
    assert(!LCA->is_leaf());
    if(dst==LCA||src->parent==LCA){
        Fitch_Sankoff::set_internal_score(*LCA, out->scores, offset, src);
    }else{
        Fitch_Sankoff::set_internal_score(*LCA, out->scores, offset);
    }
    return out->scores.back()[par_nuc_idx]-original_parsimony_score;
}

void Profitable_Moves_Enumerator::operator() (Possible_Move* in)const{
    MAT::Mutations_Collection mutations;
    Move *out = new Move;
    MAT::Node* LCA=get_mutation_path(mutations, in->src, in->dst,out->path);
    assert(LCA!=in->src);
    if(mutations.empty()){
        delete in;
        delete out;
        return;
    }
    std::pair<size_t, size_t> range=Fitch_Sankoff::dfs_range(LCA,dfs_ordered_nodes);
    std::vector<Fitch_Sankoff::States_Type> original_states;
    std::vector<char> LCA_parent_states=get_original_states(dfs_ordered_nodes, range, mutations,original_states);
    std::atomic_int score_changes(0);
    std::vector<Fitch_Sankoff_Result*> unchanged_states(mutations.size());
    std::vector<Fitch_Sankoff_Result*> moved_states(mutations.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,mutations.size()),[&](tbb::blocked_range<size_t> r){
        for (size_t i=r.begin();i<r.end() ; i++) {
            unchanged_states[i]=new Fitch_Sankoff_Result;
            Fitch_Sankoff::sankoff_backward_pass(range, dfs_ordered_nodes, unchanged_states[i]->scores,original_states[i]);
            Fitch_Sankoff_Result* patched=new Fitch_Sankoff_Result(*unchanged_states[i]);
            int this_change=calculate_parsimony_score_change(range,patched,in->src,in->dst,LCA_parent_states[i],LCA);
            moved_states[i]=patched;
            score_changes.fetch_add(this_change,std::memory_order_relaxed);
        }
    });
/*
#ifndef NDEBUG
    MAT::Tree new_tree(in->src->tree);
    MAT::Node* new_src=new_tree.get_node(in->src->identifier);
    MAT::Node* new_LCA=new_tree.get_node(LCA->identifier);
    MAT::Node* old_parent=new_src->parent;
    remove_child(new_src);
    MAT::Node* new_dst=new_tree.get_node(in->dst->identifier);
    MAT::Node* added=new_dst->add_child(new_src);
    std::vector<MAT::Node*> new_dfs_ordered_nodes=new_tree.depth_first_expansion();
    std::pair<size_t, size_t> new_range=Fitch_Sankoff::dfs_range(new_LCA,new_dfs_ordered_nodes);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,mutations.size()),[&](tbb::blocked_range<size_t> r){
        for (size_t i=r.begin();i<r.end() ; i++) {
            std::unordered_map<std::string, std::array<int,4>> original_scores;
            for(auto&& score:moved_states[i]->scores){
                if(added&&score.node==in->dst){
                    original_scores.emplace(new_dst->identifier,score.score);
                }else{
                    original_scores.emplace(score.node->identifier,score.score);
                }
            }
            Fitch_Sankoff::Scores_Type new_score;
            Fitch_Sankoff::States_Type new_original_states(new_range.second-new_range.first);
            for(auto&& state:original_states[i]){
                auto new_corresponding_node=new_tree.get_node(state.node->identifier);
                new_original_states[new_corresponding_node->index-new_range.first].state=state.state;
                new_original_states[new_corresponding_node->index-new_range.first].node=new_corresponding_node;
            }
            if(added){
                new_original_states[new_dst->index - new_range.first].state =original_states[i][in->dst->index - range.first].state;
                new_original_states[new_dst->index - new_range.first].node =added;
            }
            Fitch_Sankoff::sankoff_backward_pass(new_range, new_dfs_ordered_nodes, new_score, new_original_states);
            for(auto&& score: new_score){
                if(added==score.node){
                    continue;
                }
                if(original_scores[score.node->identifier]!=score.score){
                    fprintf(stderr,"score mismatch at original tree node index %zu \n",in->src->tree->get_node(score.node->identifier)->index);
                    assert(false);
                }
            }
        }
    });
#endif
*/
    if (score_changes<0) {
        out->LCA=LCA;
        out->score_change=score_changes;
        out->src=in->src;
        out->dst=in->dst;
        for (size_t i=0; i<mutations.size(); i++) {
            moved_states[i]->range=range;
            moved_states[i]->mutation=mutations[i];
            moved_states[i]->LCA_parent_state=LCA_parent_states[i];
            moved_states[i]->original_state.swap(original_states[i]);
        }
        out->states.swap(moved_states);
        profitable_moves.push_back(out);
    }else{
        for(auto c:moved_states){
            delete c;
        }
        delete out;
    }
    for(auto c:unchanged_states){
        delete c;
    }
    delete in;
}