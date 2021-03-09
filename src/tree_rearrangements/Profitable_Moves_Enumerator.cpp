#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/Fitch_Sankoff.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <ios>
#include <iterator>
#include <list>
#include <memory>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/queuing_rw_mutex.h>
#include <tbb/spin_mutex.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
int count_mutation(MAT::Node* node,int pos);
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

static std::pair<int,int> calculate_parsimony_score_change(std::pair<size_t, size_t>& range, Fitch_Sankoff::Scores_Type& out, MAT::Node* src, MAT::Node* dst,char par_nuc, MAT::Node* LCA){
    assert(LCA->index==range.first);
    int offset=range.second-1;
    char par_nuc_idx=one_hot_to_two_bit(par_nuc);
    int original_parsimony_score=Fitch_Sankoff::get_child_score_on_par_nuc(par_nuc_idx, out.back()).first;
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
    return std::make_pair(Fitch_Sankoff::get_child_score_on_par_nuc(par_nuc_idx, out.back()).first,original_parsimony_score);
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
    std::vector<Profitable_Move*> &tie;
    int& best_score_change;
    Candidate_Moves *all_moves;
    std::mutex& mutex;
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
            
            #ifdef DETAIL_DEBUG_PARSIMONY_SCORE_MATCH
            MAT::Tree new_tree;
            MAT::Node *new_LCA = new Mutation_Annotated_Tree::Node(*LCA,nullptr,&new_tree);
            new_tree.root=new_LCA;
            MAT::Node *new_src = new_tree.get_node(all_moves->src->identifier);
            //MAT::Node *new_LCA = new_tree.get_node(LCA->identifier);
            MAT::Node *original_src_parent_in_new_tree = new_src->parent;
            MAT::Node *new_dst = new_tree.get_node(in->dst->identifier);
            remove_child(new_src);
            MAT::Node *added = new_dst->add_child(new_src);
            while (original_src_parent_in_new_tree->is_leaf()) {
                MAT::Node* parent_of_parent=original_src_parent_in_new_tree->parent;
                remove_child(original_src_parent_in_new_tree);
                original_src_parent_in_new_tree=parent_of_parent;
            }
            
            std::vector<MAT::Node *> new_dfs_ordered_nodes =
                new_tree.depth_first_expansion();
            std::pair<size_t, size_t> new_range =
                Fitch_Sankoff::dfs_range(new_LCA, new_dfs_ordered_nodes);
            #endif
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, in->FS_results.size()),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t mut_idx = r.begin(); mut_idx < r.end();
                         mut_idx++) {
                        Fitch_Sankoff::Scores_Type &patched =
                            moved_states[mut_idx].scores;
                        patched.clear();
                        std::shared_ptr<Fitch_Sankoff_Result>& parent_fs_result_managing_ptr=
                            in->FS_results[mut_idx];
                        Fitch_Sankoff_Result *parent_fs_result =parent_fs_result_managing_ptr.get();
                        moved_states[mut_idx].mutation=parent_fs_result->mutation;
                        bool clear_after=parent_fs_result->scores.empty();
                        if(clear_after){
                            Fitch_Sankoff::sankoff_backward_pass(parent_fs_result->range, dfs_ordered_nodes,parent_fs_result->scores,original_states,parent_fs_result->mutation);
                        }
                        copy_scores(*parent_fs_result, patched, range);
                        char LCA_parent_state =(LCA==dfs_ordered_nodes[parent_fs_result->range.first])?parent_fs_result->LCA_parent_state:get_genotype(LCA->parent, parent_fs_result->mutation);
                        moved_states[mut_idx].LCA_parent_state=LCA_parent_state;
                        moved_states[mut_idx].LCA_parent_state_original_tree=LCA_parent_state;
                        auto new_old_score = calculate_parsimony_score_change(
                            range, patched, all_moves->src, in->dst,
                            LCA_parent_state, LCA);
                        int this_change =
                            new_old_score.first - new_old_score.second;
                        moved_states[mut_idx].optimized_score=new_old_score.first;
                        #ifndef NDEBUG
                        moved_states[mut_idx].original_topology_score=new_old_score.second;
                        #endif
#ifdef DETAIL_DEBUG_PARSIMONY_SCORE_MATCH
                        int parsimony_score_in_original_tree = 0;
                        for (size_t idx = range.first; idx < range.second;
                             idx++) {
                            parsimony_score_in_original_tree += count_mutation(
                                dfs_ordered_nodes[idx],
                                parent_fs_result->mutation.position);
                        }
                        assert(parsimony_score_in_original_tree>=new_old_score.second);
                        std::unordered_map<std::string, std::array<int, 4>>
                            scores_to_check;
                        for (auto &&score : patched) {
                            if (added && score.node == in->dst) {
                                scores_to_check.emplace(new_dst->identifier,
                                                        score.score);
                            } else {
                                scores_to_check.emplace(score.node->identifier,
                                                        score.score);
                            }
                        }
                        Fitch_Sankoff::Scores_Type new_score;

                        Fitch_Sankoff::sankoff_backward_pass(
                            new_range, new_dfs_ordered_nodes, new_score,
                            original_states, parent_fs_result->mutation);
                        for (auto &&score : new_score) {
                            if (added == score.node) {
                                continue;
                            }
                            if(scores_to_check.count(score.node->identifier)==0){
                                fprintf(stderr, "%s not found\n",score.node->identifier.c_str());
                                assert(false);
                            }
                            auto & patched_score=scores_to_check[score.node->identifier];
                            if (patched_score !=score.score) {
                                fprintf(stderr,
                                        "score mismatch at original tree node "
                                        "index %zu\n",
                                        all_moves->src->tree
                                            ->get_node(score.node->identifier)
                                            ->index);
                                assert(false);
                            }
                        }
#endif
                        parent_fs_result_managing_ptr.reset();
                        score_changes.fetch_add(this_change,
                                                std::memory_order_relaxed);
                    }
                });
#ifdef DETAIL_DEBUG_PARSIMONY_SCORE_MATCH
            new_tree.delete_nodes();
#endif
            if (score_changes<=best_score_change) {
                Profitable_Move *out = new Profitable_Move;
                out->LCA = LCA;
                out->score_change = score_changes;
                out->src = all_moves->src;
                out->dst = in->dst;
                out->path.swap(in->path);
                out->range = range;
                out->states.swap(moved_states);
                std::vector<Profitable_Move*> temp({out});
                temp.reserve(all_moves->moves.size());
                {
		    std::lock_guard<std::mutex> lock(mutex);
                    if (score_changes==best_score_change) {
                        tie.push_back(out);
                        continue;
                    }else if(score_changes<best_score_change){
			best_score_change=score_changes;
                        tie.swap(temp);
                    }
                }
                for(Profitable_Move* move:temp){
                    delete move;
                }
            }
        }
    }
};

Profitable_Moves_From_One_Source* Profitable_Moves_Enumerator::operator()(Candidate_Moves * in)const{
    Profitable_Moves_From_One_Source* out=new Profitable_Moves_From_One_Source;
    out->src=in->src;
    std::mutex mutex;
    int best_score_change=-1;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,in->moves.size()),Test_Move_Profitable{dfs_ordered_nodes,out->profitable_moves,best_score_change,in,mutex,original_states});
    delete in;
    return out;
}
