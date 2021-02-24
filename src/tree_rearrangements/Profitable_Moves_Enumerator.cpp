#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
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
#include <tbb/queuing_rw_mutex.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static void remove_child(MAT::Node* child_to_remove){
        auto parent=child_to_remove->parent;
        auto iter=std::find(parent->children.begin(),parent->children.end(),child_to_remove);
        assert(iter!=parent->children.end());
        parent->children.erase(iter);
}

static void patch_sankoff_result(size_t offset, Fitch_Sankoff::Scores_Type& out, MAT::Node* src, MAT::Node* dst,size_t end_idx){
    if(dst->is_leaf()){
        assert(out[offset-dst->index].node==dst);
        Fitch_Sankoff::Score_Type leaf_score(out[offset-dst->index]);
        Fitch_Sankoff::set_internal_score(*dst, out, offset, src,&leaf_score);
    }else{
        Fitch_Sankoff::set_internal_score(*dst, out, offset, src);
    }
    MAT::Node* changing_node=dst->parent;
    while (changing_node&&changing_node->index>end_idx) {
        Fitch_Sankoff::set_internal_score(*changing_node, out, offset);
        changing_node=changing_node->parent;
    }
}

static int calculate_parsimony_score_change(std::pair<size_t, size_t>& range, Fitch_Sankoff::Scores_Type& out, MAT::Node* src, MAT::Node* dst,char par_nuc, MAT::Node* LCA){
    assert(LCA->index==range.first);
    int offset=range.second-1;
    char par_nuc_idx=one_hot_to_two_bit(par_nuc);
    int original_parsimony_score=out.back()[par_nuc_idx];
    size_t end_idx=range.first;
    patch_sankoff_result(offset, out, src, src->parent,end_idx);
    patch_sankoff_result(offset, out, src, dst,end_idx);
    //Deal with LCA separately, as it is influenced by both src and dst, especially the case when src is child of LCA, then src need to be processed, or dst is LCA itself (move to parent of parent), dst is to be processed
    assert(!LCA->is_leaf());
    if(dst==LCA||src->parent==LCA){
        Fitch_Sankoff::set_internal_score(*LCA, out, offset, src);
    }else{
        Fitch_Sankoff::set_internal_score(*LCA, out, offset);
    }
    return out.back()[par_nuc_idx]-original_parsimony_score;
}

static void copy_scores(const Fitch_Sankoff_Result &ori,
                       Fitch_Sankoff::Scores_Type &dst,
                       std::pair<size_t, size_t> &new_range) {
    assert(new_range.first >= ori.range.first);
    assert(new_range.second <= ori.range.second);
    std::copy(ori.scores.begin() + ori.range.second - new_range.second,
              ori.scores.end() - new_range.first + ori.range.first,
              std::back_inserter(dst));
    assert(dst.front().node->index==new_range.second-1);
    assert(dst.back().node->index==new_range.first);
}
/*
static void copy_states(const Fitch_Sankoff_Result &ori,
                       Fitch_Sankoff::States_Type &dst,
                       std::pair<size_t, size_t> &new_range) {
    assert(new_range.first >= ori.range.first);
    assert(new_range.second <= ori.range.second);
    std::copy(ori.scores.begin() - ori.range.first - new_range.first,
              ori.scores.end() + new_range.second - ori.range.second,
              std::back_inserter(dst));
    assert(dst.front().node->index==new_range.first);
    assert(dst.back().node->index==new_range.second -1);
}
*/
/*
static char extract_LCA_parent_state(const Fitch_Sankoff_Result &ori,
                                     MAT::Node *LCA) {
    if (LCA == ori.original_state->front().node) {
        return ori.LCA_parent_state;
    } else {
        size_t LCA_parent_idx = LCA->parent->index;
        assert(ori.range.first <= LCA_parent_idx &&
               ori.range.second >= LCA_parent_idx);
        const Fitch_Sankoff::States_Type &temp = *(ori.original_state);
        const Fitch_Sankoff::State_Type &parent_state =
            temp[LCA_parent_idx - ori.range.first];
        assert(parent_state.node == LCA->parent);
        return parent_state.state;
    }
}
*/
struct Test_Move_Profitable {
    std::vector<MAT::Node *> &dfs_ordered_nodes;
    tbb::concurrent_vector<Profitable_Move *> &profitable_moves;
    Candidate_Moves *all_moves;
    tbb::queuing_rw_mutex& mutex;
    const Original_State_t& original_states;
    void operator()(tbb::blocked_range<size_t> &range) const {
        for (size_t move_idx = range.begin(); move_idx < range.end();
             move_idx++) {
            Move_info *in = &(all_moves->moves[move_idx]);
            MAT::Node *LCA = in->LCA;
            std::pair<size_t, size_t> range =
                Fitch_Sankoff::dfs_range(LCA, dfs_ordered_nodes);
            std::atomic_int score_changes(0);
            std::vector<Fitch_Sankoff_Result_Final> moved_states(
                in->FS_results.size());
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, in->FS_results.size()),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t mut_idx = r.begin(); mut_idx < r.end();
                         mut_idx++) {
                        Fitch_Sankoff::Scores_Type &patched =
                            moved_states[mut_idx].scores;
                        patched.clear();
                        Fitch_Sankoff_Result *parent_fs_result =
                            in->FS_results[mut_idx];
                        bool clear_after=parent_fs_result->scores.empty();
                        if(clear_after){
                            Fitch_Sankoff::sankoff_backward_pass(parent_fs_result->range, dfs_ordered_nodes,parent_fs_result->scores,original_states,parent_fs_result->mutation,parent_fs_result->LCA_parent_state);
                            
                        }
                        copy_scores(*parent_fs_result, patched, range);
                        char LCA_parent_state =(LCA==dfs_ordered_nodes[parent_fs_result->range.first])?parent_fs_result->LCA_parent_state:get_genotype(LCA->parent, parent_fs_result->mutation);
                        moved_states[mut_idx].LCA_parent_state=LCA_parent_state;
                        int this_change = calculate_parsimony_score_change(
                            range, patched, all_moves->src, in->dst,
                            LCA_parent_state, LCA);
                        if(clear_after){
                            parent_fs_result->scores.clear();
                            parent_fs_result->scores.shrink_to_fit();
                        }
                        score_changes.fetch_add(this_change,
                                                std::memory_order_relaxed);
                    }
                });
            /*
            #ifndef NDEBUG
                MAT::Tree new_tree(all_moves->src->tree);
                MAT::Node*
            new_src=new_tree.get_node(all_moves->src->identifier); MAT::Node*
            new_LCA=new_tree.get_node(LCA->identifier); MAT::Node*
            old_parent=new_src->parent; remove_child(new_src); MAT::Node*
            new_dst=new_tree.get_node(in->dst->identifier); MAT::Node*
            added=new_dst->add_child(new_src); std::vector<MAT::Node*>
            new_dfs_ordered_nodes=new_tree.depth_first_expansion();
                std::pair<size_t, size_t>
            new_range=Fitch_Sankoff::dfs_range(new_LCA,new_dfs_ordered_nodes);
                tbb::parallel_for(tbb::blocked_range<size_t>(0,mutations.size()),[&](tbb::blocked_range<size_t>
            r){ for (size_t i=r.begin();i<r.end() ; i++) {
                        std::unordered_map<std::string, std::array<int,4>>
            original_scores; for(auto&& score:moved_states[i]->scores){
                            if(added&&score.node==in->dst){
                                original_scores.emplace(new_dst->identifier,score.score);
                            }else{
                                original_scores.emplace(score.node->identifier,score.score);
                            }
                        }
                        Fitch_Sankoff::Scores_Type new_score;
                        Fitch_Sankoff::States_Type
            new_original_states(new_range.second-new_range.first); for(auto&&
            state:original_states[i]){ auto
            new_corresponding_node=new_tree.get_node(state.node->identifier);
                            new_original_states[new_corresponding_node->index-new_range.first].state=state.state;
                            new_original_states[new_corresponding_node->index-new_range.first].node=new_corresponding_node;
                        }
                        if(added){
                            new_original_states[new_dst->index -
            new_range.first].state =original_states[i][in->dst->index -
            range.first].state; new_original_states[new_dst->index -
            new_range.first].node =added;
                        }
                        Fitch_Sankoff::sankoff_backward_pass(new_range,
            new_dfs_ordered_nodes, new_score, new_original_states); for(auto&&
            score: new_score){ if(added==score.node){ continue;
                            }
                            if(original_scores[score.node->identifier]!=score.score){
                                fprintf(stderr,"score mismatch at original tree
            node index %zu
            \n",all_moves->src->tree->get_node(score.node->identifier)->index);
                                assert(false);
                            }
                        }
                    }
                });
            #endif
            */
            if (score_changes < 0) {
                Profitable_Move *out = new Profitable_Move;
                out->LCA = LCA;
                out->score_change = score_changes;
                out->src = all_moves->src;
                out->dst = in->dst;
                out->path.swap(in->path);
                out->range = range;
                for (size_t i = 0; i < in->FS_results.size(); i++) {
                    moved_states[i].mutation = in->FS_results[i]->mutation;
                }
                out->states.swap(moved_states);
                {
                    tbb::queuing_rw_mutex::scoped_lock lock(mutex,false);
                    profitable_moves.push_back(out);
                }
            }
        }
    }
};

void Profitable_Moves_Enumerator::operator()(Candidate_Moves * in)const{
    tbb::parallel_for(tbb::blocked_range<size_t>(0,in->moves.size()),Test_Move_Profitable{dfs_ordered_nodes,profitable_moves,in,mutex,original_states});
    delete in;
}