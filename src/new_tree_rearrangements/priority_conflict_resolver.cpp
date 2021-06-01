#include "priority_conflict_resolver.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <mutex>
#include <utility>
#include <vector>

bool Conflict_Resolver::check_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move)const{
    int best_score=0;
    candidate_move->apply_nodes([&best_score,this](MAT::Node* node){
        best_score=std::min(best_score,potential_crosses[node->bfs_index].parsimony_score_change.load(std::memory_order_acquire));
    });
    if (candidate_move->score_change<best_score) {
        return true;
    }
    return false;
}

static void remove_move(Cross_t &potential_crosses, const Profitable_Moves_ptr_t& other_move,
                        MAT::Node *exclude) {

    other_move->apply_nodes([&](MAT::Node *other_nodes_in_path) {
        if (other_nodes_in_path != exclude) {
            auto& other_node_conflict = potential_crosses[other_nodes_in_path->bfs_index];
            assert(other_node_conflict.moves.size() == 1 &&
               other_node_conflict.moves[0] == other_move);
            other_node_conflict.parsimony_score_change.store(0,std::memory_order_release);
            other_node_conflict.moves.clear();
        }
    });
}

bool Conflict_Resolver::register_single_move_no_conflict(
    Profitable_Moves_ptr_t& candidate_move)  {
    if (candidate_move->src->bfs_index==5319&&candidate_move->dst_to_LCA[0]->bfs_index==1285) {
        fputc('a',stderr);
    }
    if (candidate_move->src->bfs_index==19290&&candidate_move->dst_to_LCA[0]->bfs_index==9626) {
        fputc('a',stderr);
    }
    candidate_move->apply_nodes([&](MAT::Node* node) {
        Conflict_Set& iter = potential_crosses[node->bfs_index];
        iter.parsimony_score_change.store(candidate_move->score_change,std::memory_order_release);
            for (Profitable_Moves_ptr_t& other_move : iter.moves) {
                assert(other_move->score_change >=
                       candidate_move->score_change);
                remove_move(potential_crosses, other_move, node);
                assert(candidate_move!=other_move);
                #ifndef NDEBUG
                //nodes_inside--;
                #endif
                #ifdef CHECK_LEAK
                other_move->destructed=true;
                #else
                delete other_move;
                #endif
            }
            iter.moves.clear();
            iter.moves.push_back(candidate_move);
            assert(candidate_move->score_change==iter.parsimony_score_change);
    });
    assert(potential_crosses[candidate_move->get_src()->bfs_index].moves[0]=candidate_move);
    assert(potential_crosses[candidate_move->get_dst()->bfs_index].moves[0]=candidate_move);
    return true;
}

char Conflict_Resolver::operator()(std::vector<Profitable_Moves_ptr_t>& candidate_move){
    char ret=0;
    Profitable_Moves_ptr_t selected_move=nullptr;
    for (Profitable_Moves_ptr_t& move : candidate_move) {
#ifdef CONFLICT_RESOLVER_DEBUG
        fprintf(log, "Trying %s to %s\n",move->src->identifier.c_str(),move->get_dst()->identifier.c_str());
#endif
        //fflush(log);
        if (check_single_move_no_conflict(move)) {
            std::lock_guard<std::mutex> lk(register_lock);
            if (check_single_move_no_conflict(move)){
            register_single_move_no_conflict(move);
            ret =1;
            selected_move = move;
            break;}
        }
    }

    if(!selected_move&&(!candidate_move.empty())){
        deferred_moves.reserve(deferred_moves.size()+candidate_move.size());
        for (Profitable_Moves_ptr_t move : candidate_move) {
            deferred_moves.emplace_back(move->src->identifier,move->get_dst()->identifier);
        }
    }
    for (Profitable_Moves_ptr_t move : candidate_move) {
        if (move!=selected_move) {
            #ifdef CHECK_LEAK
            move->destructed=true;
            #else
            delete move;
            #endif
        }
    }
    return ret;
}

void Conflict_Resolver::schedule_moves( std::vector<Profitable_Moves_ptr_t>& out){
    #ifdef CONFLICT_RESOLVER_DEBUG
    fputc('\n', log);
    fflush(log);
    std::unordered_map<size_t,std::pair<size_t,Profitable_Moves_ptr_t>> pushed_moves;
    #endif
    for (int idx=potential_crosses.size()-1; idx>=0; idx--) {
        auto& moves_on_this_node=potential_crosses[idx];
        if (!moves_on_this_node.moves.empty()) {
            const Profitable_Moves_ptr_t this_move=moves_on_this_node.moves[0];
            if ((int)this_move->get_src()->bfs_index>=idx&&(int)this_move->get_dst()->bfs_index>=idx) {
                assert(this_move->get_src()->bfs_index==idx||this_move->get_dst()->bfs_index==idx);
    #ifdef CONFLICT_RESOLVER_DEBUG
        auto src_res=pushed_moves.emplace(this_move->src->bfs_index,std::make_pair(idx, this_move));
        if (!src_res.second) {
            auto other_move=*(src_res.first);
            fprintf(stderr, "another move pushed at %zu \n",other_move.first);
            assert(false);
        }
    #endif
                out.push_back(this_move);
                remove_move(potential_crosses, this_move, 0);
            }
        }
    }
    #ifndef NDEBUG
    for(const auto& a:potential_crosses){
        assert(a.moves.empty());
        assert(a.parsimony_score_change==0);
    }
    #endif
}
