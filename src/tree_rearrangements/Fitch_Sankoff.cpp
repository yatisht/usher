#include <algorithm>
#include <array>
#include <cstdio>
#include <limits.h>
#include <string>
#include <utility>
#include <vector>
#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"

typedef Fitch_Sankoff::Score_Type Score_Type;
typedef Fitch_Sankoff::Scores_Type Scores_Type;
typedef Fitch_Sankoff::State_Type State_Type;
typedef Fitch_Sankoff::States_Type States_Type;

static void check_onehot(char in){
  assert(__builtin_popcount(in)==1);
  assert((in&0xf0)==0);  
}

char get_genotype( MAT::Node* node, const Mutation_Annotated_Tree::Mutation& m){
    while (node)
    {
        auto iter=node->mutations.find(m);
        if(iter!=node->mutations.end()){
            return iter->mut_nuc;
        }
        node=node->parent;
    }
    return m.ref_nuc;

}
#ifndef NDEBUG
static char get_genotype_with_regrafting(MAT::Node* node, const Mutation_Annotated_Tree::Mutation& m, MAT::Node* src, MAT::Node* dst){
    while (node)
    {
        auto iter=node->mutations.find(m);
        if(iter!=node->mutations.end()){
            return iter->mut_nuc;
        }
        node=(node==src)?dst:node->parent;
    }
    return m.ref_nuc;

}
#endif
static void fill_nuc(char nuc, Score_Type &out) {
    check_onehot(nuc);
    for (auto i = 0; i < 4; i++) {
        if ((1 << i) & nuc) {
            out[i] = 0;
        } else {
            out[i] =
                2; // 2 is sufficient for not being chosen, and avoid overflow.
        }
    }
}
static void set_leaf_score(MAT::Node &this_node,
                           Score_Type &out,char ori_state) {
    assert(out.node==&this_node);
    //assert(ori_state==get_genotype(&this_node, pos));
    fill_nuc(ori_state, out);
}

std::pair<int, char>
Fitch_Sankoff::get_child_score_on_par_nuc(char par_nuc,
                           const Score_Type &child_scores) {
    assert(par_nuc<4);
    int min_score = INT_MAX;
    char min_state = -1;
    for (char child_nuc = 0; child_nuc < 4; child_nuc++) {
        auto child_score = child_scores[child_nuc];
        auto this_score =
            (child_nuc == par_nuc) ? child_score : child_score + 1;
        if (min_score > this_score) {
            min_score = this_score;
            min_state = child_nuc;
        }
    }
    assert(min_state != -1);
    return std::make_pair(min_score, min_state);
}

static void patch_children(std::vector<MAT::Node*>& in, MAT::Node* changed_child){
    std::vector<MAT::Node*>::iterator iter=std::find(in.begin(),in.end(),changed_child);
    if (iter==in.end()) {
        in.push_back(changed_child);
    }else {
        in.erase(iter);
    }
}

void Fitch_Sankoff::set_internal_score(const MAT::Node &this_node, Scores_Type &out,
                        const int start_idx,MAT::Node* changed_child,Score_Type* leaf_score) {
    // make sure no off by one
    //assert(out.size() - 1 == start_idx - this_node.index);
    for (char par_nuc = 0; par_nuc < 4; par_nuc++) {
        int score = 0;
        // Children are visited after parent in pre-order DFS
        std::vector<MAT::Node*>::const_iterator iter;
        std::vector<MAT::Node*>::const_iterator end;
        std::vector<MAT::Node*> patched_children;
        if(changed_child){
             patched_children=this_node.children;
             if(patched_children.empty()){
                patched_children.push_back(changed_child);
             }else{
                patch_children(patched_children, changed_child);
             }
             iter=patched_children.begin();
             end=patched_children.end();
        }else{
            iter=this_node.children.begin();
            end=this_node.children.end();
        }
        for (;iter<end;iter++) {
            auto child=*iter;
            size_t child_idx = start_idx - child->index;
            // make sure indexes match
            assert(child_idx <= out.size() - 2);
            assert(out[child_idx].node == child);
            // Make sure child values are computed
            assert(out[child_idx][0] != -1);

            auto temp =
                Fitch_Sankoff::get_child_score_on_par_nuc(par_nuc, out[child_idx]);
            score += temp.first;
        }
        if(leaf_score) score+=Fitch_Sankoff::get_child_score_on_par_nuc(par_nuc, *leaf_score).first;
        /* In the whole subtree backward pass, the node whose score are being calculated are always at the end, but when patching the result, it is not necessarily the case.
        assert(out.back().node == &this_node);
        out.back()[par_nuc] = score;
        */
        size_t this_node_idx = start_idx - this_node.index;
        assert(out[this_node_idx].node==&this_node);
        out[this_node_idx][par_nuc]=score;
    }
}

std::pair<size_t, size_t> Fitch_Sankoff::dfs_range(const MAT::Node *start,const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    size_t start_idx = start->index;
    size_t stop_idx = INT_MAX;
    if (!start->parent) {
        return std::make_pair(start_idx, dfs_ordered_nodes.size());
    }
    for (auto sibling : start->parent->children) {
        auto sibling_idx = sibling->index;
        if (sibling_idx > start_idx && sibling_idx < stop_idx) {
            stop_idx = sibling_idx;
        }
    }
    if (stop_idx == INT_MAX) {
        stop_idx = Fitch_Sankoff::dfs_range(start->parent,dfs_ordered_nodes).second;
    }
    assert (stop_idx != INT_MAX);
    assert(check_grand_parent(dfs_ordered_nodes[stop_idx-1], start));
    assert(stop_idx==dfs_ordered_nodes.size()||!check_grand_parent(dfs_ordered_nodes[stop_idx], start));
    assert(start_idx < stop_idx);
    return std::make_pair(start_idx, stop_idx);
}

char get_original_state(MAT::Node* node, const Original_State_t& original_states, const MAT::Mutation& mutation){
    const Mutation_Set& this_node_muts=original_states.at(node->identifier); 
    auto mut_iter=this_node_muts.find(mutation);
    return (mut_iter==this_node_muts.end())?mutation.ref_nuc:mut_iter->mut_nuc;
}
int Fitch_Sankoff::sankoff_backward_pass(const std::pair<size_t, size_t> &range,
                           const std::vector<MAT::Node *> &dfs_ordered_nodes,
                           Scores_Type &scores,const Original_State_t& original_state,const MAT::Mutation& mutation,char starting_node_parent_state) {
    //Going from the end to start
    size_t start_idx = range.second-1;
    assert(scores.empty());
    scores.reserve(range.second - range.first);
    for (auto iter = dfs_ordered_nodes.begin() + range.second-1;
         iter >= dfs_ordered_nodes.begin() + range.first; iter--) {
#ifndef NDEBUG
        scores.emplace_back(*iter);
#else
        scores.emplace_back();
#endif
        Score_Type& score_array = scores.back();
        if ((*iter)->is_leaf()) {
            
            set_leaf_score(**iter, score_array,get_original_state(*iter,original_state,mutation));
        } else {
            set_internal_score(**iter, scores, start_idx);
        }
    }
    assert(scores.size()== (range.second - range.first));
    return get_child_score_on_par_nuc(one_hot_to_two_bit(starting_node_parent_state), scores.back()).first;
}

static void set_mutation(MAT::Node *node, char state, char par_state, const MAT::Mutation& mutationOri) {    
    check_onehot(par_state);
    check_onehot(state);
    if (par_state == state) {
        node->mutations.remove(mutationOri.position);
        assert(node->mutations.find(mutationOri.position)==node->mutations.end());
    } else {
        MAT::Mutation mutation(mutationOri);
        mutation.par_nuc = par_state;
        mutation.mut_nuc = state;
        node->mutations.insert(mutation, 0);
    }
}

#define get_score_idx(dfs_idx) (score_offset-(dfs_idx))
#define get_state_idx(dfs_idx) ((dfs_idx)-state_offset)

static void fill_dst_states(MAT::Node* to_fill,const size_t score_offset,const size_t state_offset,const MAT::Mutation& mutation,Scores_Type& scores,States_Type& states,MAT::Node* to_move,MAT::Node* dst){
    assert(to_fill==dst||(!to_fill->is_leaf()));

    size_t this_state_idx = get_state_idx(to_fill->index);
    assert(states[this_state_idx].node==to_fill);

    if(states[this_state_idx]>=4){
        MAT::Node *parent = to_fill->parent;
        size_t parent_state_idx = get_state_idx(parent->index);
        assert(states[parent_state_idx].node == parent);

        fill_dst_states(parent, score_offset, state_offset, mutation, scores,
                        states, to_move, dst);
        assert((1 << states[parent_state_idx]) ==
               get_genotype_with_regrafting(parent, mutation, to_move, dst));

        size_t this_score_idx = get_score_idx(to_fill->index);
        assert(scores[this_score_idx].node == to_fill);

        states[this_state_idx] =
            Fitch_Sankoff::get_child_score_on_par_nuc(states[parent_state_idx],
                                                      scores[this_score_idx])
                .second;
        set_mutation(to_fill, 1 << states[this_state_idx],
                     1 << states[parent_state_idx], mutation);
    }
}

void Fitch_Sankoff::sankoff_forward_pass(const std::pair<size_t, size_t> &range,
                          std::vector<MAT::Node *> &dfs_ordered_nodes,const MAT::Mutation &mutation,const Original_State_t& original_states,
                          Scores_Type& scores, char starting_node_parent_state,MAT::Node* to_move,MAT::Node* dst,MAT::Node* new_leaf
) {
    // Score vector runs backward
    //first node->last element
    //index=size-1-(dfs_index-first_element_index)
    size_t score_offset=-1+range.second;

    //States in 2bit format , original_states vector in one-hot format goes forward
    #ifndef NDEBUG
    States_Type states;
    states.reserve(range.second-range.first);
    for (size_t node_idx=range.first; node_idx<range.second; node_idx++) {
        states.push_back(State_Type(32,dfs_ordered_nodes[node_idx]));
    }
    #else
    States_Type states(range.second-range.first,32);
    #endif
    size_t state_offset=range.first;

    //starting node have states passed in as argument, and I don't want to have different index for original_states and state, so here goes the special treatment for the starting node.
    assert(get_genotype(dfs_ordered_nodes[range.first]->parent, mutation)==starting_node_parent_state);

    //in 2 bit format needed by get_score_by_par_nuc
    starting_node_parent_state=one_hot_to_two_bit(starting_node_parent_state);

    assert(states[0].node==dfs_ordered_nodes[range.first]);
    Score_Type& last_score=scores[range.second-range.first-1];
    states[0]=get_child_score_on_par_nuc(starting_node_parent_state, last_score).second;
    assert(last_score.node==dfs_ordered_nodes[range.first]);
    set_mutation(dfs_ordered_nodes[range.first], 1<<states[0], 1<<starting_node_parent_state, mutation);

    for (size_t dfs_idx=range.first+1; dfs_idx<range.second; dfs_idx++) {
        //Index conversions
        auto this_node=dfs_ordered_nodes[dfs_idx];
        size_t score_idx=get_score_idx(dfs_idx);
        assert(this_node==scores[score_idx].node);

        MAT::Node* parent_node;

        if(this_node==to_move){
            fill_dst_states(dst, score_offset, state_offset, mutation, scores, states, to_move, dst);
            parent_node=dst;
        }else{
            parent_node=this_node->parent;
        }

        size_t parent_state_idx=get_state_idx(parent_node->index);
        assert(states[parent_state_idx].node==parent_node);
        char parent_state=states[parent_state_idx];
        assert(get_genotype_with_regrafting(parent_node, mutation,to_move,dst)==(1<<states[parent_state_idx]));

        size_t this_state_idx=get_state_idx(dfs_idx);
        assert(states[this_state_idx].node==this_node);

        //If to_move is moved to dst which is a leaf node, it will become an internal node, so treat it as an internal node
        if(this_node->is_leaf()&&this_node!=dst) {
            char original_state=get_original_state(this_node, original_states, mutation);
            set_mutation(this_node, original_state, 1<<parent_state, mutation);
            assert(get_genotype_with_regrafting(this_node, mutation,to_move,dst)==original_state);
            continue;
        }


        char this_state=get_child_score_on_par_nuc(parent_state, scores[score_idx]).second;

        states[this_state_idx]=this_state;
        set_mutation(this_node, 1<<this_state, 1<<parent_state, mutation);
        assert(get_genotype_with_regrafting(this_node, mutation,to_move,dst)==(1<<states[this_state_idx]));

        //deal with leaf node turned into internal node
        if(this_node->is_leaf()&&this_node==dst) {
            char original_state=get_original_state(this_node, original_states, mutation);
            set_mutation(new_leaf, original_state, 1<<this_state, mutation);
            assert(get_genotype_with_regrafting(new_leaf, mutation,to_move,dst)==original_state);
        }
    }
}
