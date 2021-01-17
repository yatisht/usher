#include "mutation_annotated_tree.hpp"
#include <array>
#include <limits.h>
#include <utility>
#include <vector>
#include "Fitch_Sankoff.hpp"

typedef Fitch_Sankoff::Score_Type Score_Type;
typedef Fitch_Sankoff::Scores_Type Scores_Type;
typedef Fitch_Sankoff::State_Type State_Type;
typedef Fitch_Sankoff::States_Type States_Type;
static char get_genotype(MAT::Node* node, const Mutation_Annotated_Tree::Mutation& m){
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
static void fill_nuc(char nuc, Score_Type &out,State_Type& state) {
    for (auto i = 0; i < 4; i++) {
        if ((1 << i) & nuc) {
            out[i] = 0;
            state=i;
        } else {
            out[i] =
                2; // 2 is sufficient for not being chosen, and avoid overflow.
        }
    }
}
static void set_leaf_score(MAT::Node &this_node, const MAT::Mutation &pos,
                           Score_Type &out,State_Type& state) {
    assert(out.node==&this_node);
    fill_nuc(get_genotype(&this_node, pos), out,state);
    
}

std::pair<int, char>
Fitch_Sankoff::get_child_score_on_par_nuc(char par_nuc,
                           Score_Type &child_scores) {
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

static void set_internal_score(const MAT::Node &this_node, Scores_Type &out,
                               const int start_idx, States_Type &states) {
    // make sure no off by one
    assert(out.size() - 1 == start_idx - this_node.index);
    for (char par_nuc = 0; par_nuc < 4; par_nuc++) {
        int score = 0;
        // Children are visited after parent in pre-order DFS
        for (auto child : this_node.children) {
            size_t child_idx = start_idx - child->index;
            // make sure indexes match
            assert(child_idx <= out.size() - 2);
            assert(out[child_idx].node == child);
            assert(states[child_idx].node == child);
            // Make sure child values are computed
            assert(out[child_idx][0] != -1);

            auto temp =
                get_child_score_on_par_nuc(par_nuc, out[child_idx]);
            score += temp.first;
            assert((states[child_idx] &(3<<(2* par_nuc)))==0||(states[child_idx] &(3<<(2* par_nuc)))==temp.second);
            states[child_idx] |= (temp.second << (2 * par_nuc));
        }
        assert(out.back().node == &this_node);
        out.back()[par_nuc] = score;
    }
}

std::pair<size_t, size_t> Fitch_Sankoff::dfs_range(const MAT::Node *start) {
    size_t start_idx = start->index;
    size_t stop_idx = INT_MAX;
    if (!start->parent) {
        return std::make_pair(start_idx, stop_idx);
    }
    for (auto sibling : start->parent->children) {
        auto sibling_idx = sibling->index;
        if (sibling_idx > start_idx && sibling_idx < stop_idx) {
            stop_idx = sibling_idx;
        }
    }
    if (stop_idx == INT_MAX) {
        stop_idx = Fitch_Sankoff::dfs_range(start->parent).second;
    }
    assert(stop_idx !=
           INT_MAX); // shouldn't be running on last node (which is a leaf node)
    assert(start_idx < stop_idx);
    return std::make_pair(start_idx, stop_idx);
}

void Fitch_Sankoff::sankoff_backward_pass(const std::pair<size_t, size_t> &range,
                           const MAT::Mutation &mutation,
                           const std::vector<MAT::Node *> &dfs_ordered_nodes,
                           Scores_Type &scores, States_Type &states) {
    //Going from the end to start
    size_t start_idx = range.second-1;
    assert(scores.empty());
    assert(states.empty());
    for (auto iter = dfs_ordered_nodes.begin() + range.second-1;
         iter >= dfs_ordered_nodes.begin() + range.first; iter--) {
#ifndef NDEBUG
        scores.emplace_back(*iter);
#else
        scores.emplace_back();
#endif

#ifndef NDEBUG
        states.emplace_back(*iter);
#else
        states.emplace_back(0);
#endif

        Score_Type& score_array = scores.back();
        if (!(*iter)->not_sample()) {
            set_leaf_score(**iter, mutation, score_array,states.back());
        } else {
            set_internal_score(**iter, scores, start_idx, states);
        }
    }
    assert(scores.size()== (range.second - range.first));
    assert(scores.size() == states.size());
}

static void set_mutation(MAT::Node *node, char state, char par_state, const MAT::Mutation& mutationOri) {
    assert(par_state==get_genotype(node->parent,mutationOri));
    if(par_state == state){
        node->mutations.remove(mutationOri.position);
        return;
    }
    MAT::Mutation mutation(mutationOri);
    mutation.par_nuc=par_state;
    mutation.mut_nuc=state;
    node->mutations.insert(mutation,0);
}



void Fitch_Sankoff::sankoff_forward_pass(const std::pair<size_t, size_t> &range,
                          States_Type &states,
                          std::vector<MAT::Node *> &dfs_ordered_nodes,const MAT::Mutation &mutation,char ancestor_state 
#ifndef NDEBUG
,Scores_Type &child_scores 
#endif
) {
    //first node->last element
    //index=size-1-(dfs_index-first_element_index)
    size_t offset=states.size()-1+range.first;
    //set the state of first node
    assert(dfs_ordered_nodes[range.first]==states.back().node);
    set_mutation(dfs_ordered_nodes[range.first],1<<states.back(),ancestor_state,mutation);
    //for the rest
    for (size_t dfs_idx=range.first+1; dfs_idx<range.second; dfs_idx++) {
        auto this_node=dfs_ordered_nodes[dfs_idx];
        size_t state_idx=offset-dfs_idx;
        assert(this_node==states[state_idx].node);

        if(!this_node->not_sample()) continue;

        size_t parent_state_idx=offset-this_node->parent->index;
        assert(states[parent_state_idx].node==this_node->parent);
        char parent_state=states[parent_state_idx];
        assert(parent_state<4);
        
        char this_state=3&(states[state_idx]>>(2*parent_state));
        assert(this_state==get_child_score_on_par_nuc(parent_state, child_scores[state_idx]).second);
        
        states[state_idx]=this_state;
        set_mutation(this_node, 1<<this_state, 1<<parent_state, mutation);
    }
}

