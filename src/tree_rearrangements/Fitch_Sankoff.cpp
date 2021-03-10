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
MAT::Node* get_parent(MAT::Node* this_node,const all_moves_bare_type& edits){
    auto iter=edits.find(this_node);
    if (iter==edits.end()) {
        return this_node->parent;
    }else {
        return iter->second;
    }
}
#ifdef DETAIL_DEBUG_FS_CORRECT
char get_genotype_with_regrafting(MAT::Node* node, const Mutation_Annotated_Tree::Mutation& m, const all_moves_bare_type& edits,MAT::Node* LCA_parent_node,char LCA_parent_state){
    while (node!=LCA_parent_node)
    {
        auto iter=node->mutations.find(m);
        if(iter!=node->mutations.end()){
            return iter->mut_nuc;
        }
        node=get_parent(node, edits);
    }
    return LCA_parent_state;

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
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(out.node==&this_node);
    #endif
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
#ifdef DETAIL_DEBUG_FS_CORRECT
static void set_internal_score_original(const MAT::Node &this_node, Scores_Type &out,
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
        #ifdef DETAIL_DEBUG_INDEX_MATCH
            assert(out[child_idx].node == child);
        #endif
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
        #ifdef DETAIL_DEBUG_INDEX_MATCH
        assert(out[this_node_idx].node==&this_node);
        #endif
        out[this_node_idx][par_nuc]=score;
    }
}
#endif
static void compare_set_arg(int op1, int op2, char flag1, char flag2,
                            int &min_out, char &arg_min_out) {
    min_out = op1;
    arg_min_out = flag1;
    if (op1 > op2) {
        min_out = op2;
        arg_min_out = flag2;
    } else if (op1 == op2) {
        arg_min_out = flag1 | flag2;
    }
}

static void find_min_4(const Fitch_Sankoff::Score_Type &in, int &min_out,
                       char &arg_min_out) {
    int first_min;
    char first_min_flag;
    compare_set_arg(in[0], in[1], 1, 2, first_min, first_min_flag);
    int second_min;
    char second_min_flag;
    compare_set_arg(in[2], in[3], 4, 8, second_min, second_min_flag);
    compare_set_arg(first_min, second_min, first_min_flag, second_min_flag,
                    min_out, arg_min_out);
    /*
    assert(min_out==*std::min_element(in.score.begin(),in.score.end()));
    for (char i=0; i<3; i++) {
        if(arg_min_out&(1<<i)){
            assert(in[i]==min_out);
        }
    }*/
}

static void proc_score(const Fitch_Sankoff::Score_Type &score,std::vector<char> &min_score_possible_nuc,std::vector<int> &min_scores){
    int min_score;
    char arg_min_one_hot;
    find_min_4(score, min_score, arg_min_one_hot);
    min_score_possible_nuc.push_back(arg_min_one_hot);
    min_scores.push_back(min_score);
}
static void proc_child(const MAT::Node *child, size_t start_idx,
                       const Fitch_Sankoff::Scores_Type &scores,
                       std::vector<char> &min_score_possible_nuc,
                       std::vector<int> &min_scores) {
    size_t child_idx = start_idx - child->index;
    /* make sure indexes match*/
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(child_idx <= scores.size() - 2);
    #endif
    const Fitch_Sankoff::Score_Type &this_score=scores[child_idx];
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(this_score.node == child);
    #endif
    /* Make sure child values are computed */
    assert(this_score[0] != -1);
    proc_score(this_score, min_score_possible_nuc, min_scores);
    
}

static void set_score_from_collated_child_score(size_t final_child_size,const std::vector<char>& min_score_possible_nuc,const std::vector<int>& min_scores,Fitch_Sankoff::Score_Type &this_node_score){
    assert(min_score_possible_nuc.size() == final_child_size);
    assert(min_scores.size() == final_child_size);
    for (char par_nuc = 0; par_nuc < 4; par_nuc++) {
        int new_score = 0;
        char par_nuc_one_hot = 1 << par_nuc;
        for (size_t child_idx = 0; child_idx < final_child_size; child_idx++) {
            new_score += min_scores[child_idx];
            if (!(par_nuc_one_hot & min_score_possible_nuc[child_idx])) {
                new_score++;
            }
        }
        this_node_score[par_nuc]=new_score;
    }
}
/*
bool patch_internal_score(const MAT::Node &this_node,
                          Fitch_Sankoff::Scores_Type &out, const int start_idx,
                          const ConfirmedMove &changed_child) {

    // make sure no off by one
    // assert(out.size() - 1 == start_idx - this_node.index);
    std::vector<char> min_score_possible_nuc;
    std::vector<int> min_scores;
    size_t final_child_size = this_node.children.size() +
                           changed_child.added.size() -
                           changed_child.removed.size();

    min_score_possible_nuc.reserve(final_child_size);
    min_scores.reserve(final_child_size);
    const auto &removed_children = changed_child.removed;
    for (MAT::Node *child : this_node.children) {
        if (std::find(removed_children.begin(), removed_children.end(),
                      child) != removed_children.end()) {
            continue;
        }
        proc_child(child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    for (MAT::Node *child : changed_child.added) {
        proc_child(child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    size_t this_node_idx = start_idx - this_node.index;
    Fitch_Sankoff::Score_Type &this_node_score = out[this_node_idx];
    assert(this_node_score.node == &this_node);
    return score_changed;
}
*/
#define set_internal_score_premeable \
std::vector<char> min_score_possible_nuc;\
std::vector<int> min_scores;\
min_score_possible_nuc.reserve(final_child_size);\
min_scores.reserve(final_child_size);\

static void set_internal_score(const MAT::Node *this_node, Scores_Type &out,const int start_idx){
    size_t final_child_size = this_node->children.size();
    set_internal_score_premeable

    #ifdef DETAIL_DEBUG_FS_CORRECT
    set_internal_score_original(*this_node, out, start_idx, nullptr, nullptr);
    Score_Type old_score=out.back();
    #endif

    for (MAT::Node *child : this_node->children) {
        proc_child(child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    // In the whole subtree backward pass, the node whose score are being calculated are always at the end
    
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(out.back().node == this_node);
    #endif
    
    set_score_from_collated_child_score(final_child_size,min_score_possible_nuc,min_scores,out.back());

    #ifdef DETAIL_DEBUG_FS_CORRECT
    assert(old_score==out.back());
    #endif
}

void Fitch_Sankoff::set_internal_score_patched(const MAT::Node &this_node, Scores_Type &out,
        const int start_idx,MAT::Node* changed_child,Score_Type* leaf_score){
    assert(!this_node.children.empty()||changed_child);
    size_t final_child_size = this_node.children.size()+1+(leaf_score!=nullptr);
    set_internal_score_premeable
    size_t this_node_idx = start_idx - this_node.index;
    Score_Type& this_score=out[this_node_idx];
    
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(this_score.node==&this_node);
    #endif
    
    #ifdef DETAIL_DEBUG_FS_CORRECT
    set_internal_score_original(this_node, out, start_idx, changed_child, leaf_score);
    Score_Type old_score=out[this_node_idx];
    #endif

    for (MAT::Node *child : this_node.children) {
        if(child==changed_child){
            assert(!leaf_score);
            final_child_size-=2;
            changed_child=nullptr;
            continue;
        }
        proc_child(child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    if(changed_child){
        proc_child(changed_child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    if(leaf_score){
        proc_score(*leaf_score, min_score_possible_nuc, min_scores);
    }
    // In the whole subtree backward pass, the node whose score are being calculated are always at the end
    set_score_from_collated_child_score(final_child_size,min_score_possible_nuc,min_scores,this_score);
    
    #ifdef DETAIL_DEBUG_FS_CORRECT
    assert(old_score==out[this_node_idx]);
    #endif
}

void Fitch_Sankoff::set_internal_score_patched(const MAT::Node &this_node, Scores_Type &out,
        const int start_idx){
    assert(!this_node.children.empty());
    size_t final_child_size = this_node.children.size();
    set_internal_score_premeable
    size_t this_node_idx = start_idx - this_node.index;
    Score_Type& this_score=out[this_node_idx];
    
    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(this_score.node==&this_node);
    #endif
    
    #ifdef DETAIL_DEBUG_FS_CORRECT
    set_internal_score_original(this_node, out, start_idx, nullptr, nullptr);
    Score_Type old_score=out[this_node_idx];
    #endif

    for (MAT::Node *child : this_node.children) {
        proc_child(child, start_idx, out, min_score_possible_nuc, min_scores);
    }
    // In the whole subtree backward pass, the node whose score are being calculated are always at the end
    set_score_from_collated_child_score(final_child_size,min_score_possible_nuc,min_scores,this_score);

    #ifdef DETAIL_DEBUG_FS_CORRECT
    assert(old_score==out[this_node_idx]);
    #endif

}
/*

*/
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
    #ifdef DETAIL_DEBUG_FS_CORRECT
    assert(check_grand_parent(dfs_ordered_nodes[stop_idx-1], start));
    assert(stop_idx==dfs_ordered_nodes.size()||!check_grand_parent(dfs_ordered_nodes[stop_idx], start));
    #endif
    assert(start_idx < stop_idx);
    return std::make_pair(start_idx, stop_idx);
}

char get_original_state(MAT::Node* node, const Original_State_t& original_states, const MAT::Mutation& mutation){
    const Mutation_Set& this_node_muts=original_states.at(node->identifier); 
    auto mut_iter=this_node_muts.find(mutation);
    return (mut_iter==this_node_muts.end())?mutation.ref_nuc:mut_iter->mut_nuc;
}
void Fitch_Sankoff::sankoff_backward_pass(const std::pair<size_t, size_t> &range,
                           const std::vector<MAT::Node *> &dfs_ordered_nodes,
                           Scores_Type &scores,const Original_State_t& original_state,const MAT::Mutation& mutation) {
    //Going from the end to start
    size_t start_idx = range.second-1;
    assert(scores.empty());
    scores.reserve(range.second - range.first);
    for (auto iter = dfs_ordered_nodes.begin() + range.second-1;
         iter >= dfs_ordered_nodes.begin() + range.first; iter--) {
#ifdef DETAIL_DEBUG_INDEX_MATCH
        scores.emplace_back(*iter);
#else
        scores.emplace_back();
#endif
        Score_Type& score_array = scores.back();
        if ((*iter)->is_leaf()) {
            
            set_leaf_score(**iter, score_array,get_original_state(*iter,original_state,mutation));
        } else {
            set_internal_score(*iter, scores, start_idx);
        }
    }
    assert(scores.size()== (range.second - range.first));
    //return get_child_score_on_par_nuc(one_hot_to_two_bit(starting_node_parent_state), scores.back()).first;
}

static void set_mutation(MAT::Node *node, char state, char par_state, const MAT::Mutation& mutationOri,States_To_Set& states_to_update) {    
    check_onehot(par_state);
    check_onehot(state);
    if (par_state == state) {
        node->mutations.remove(mutationOri.position);
    #ifdef DETAIL_DEBUG_FS_CORRECT
        assert(node->mutations.find(mutationOri.position)==node->mutations.end());
    #endif
    } else {
        MAT::Mutation mutation(mutationOri);
        mutation.par_nuc = par_state;
        mutation.mut_nuc = state;
        node->mutations.insert(mutation, 0);
    }
    auto state_iter=states_to_update.find(node);
    if (state_iter!=states_to_update.end()) {
        assert(!node->is_leaf());
        *state_iter->second=state;
    }
}

#define get_score_idx(dfs_idx) (score_offset-(dfs_idx))
#define get_state_idx(dfs_idx) ((dfs_idx)-state_offset)

static void fill_dst_states(MAT::Node *to_fill, const size_t score_offset,
                            const size_t state_offset,
                            const MAT::Mutation &mutation, Scores_Type &scores,
                            States_Type &states, MAT::Node *dst,
                            const all_moves_bare_type &edits,
                            States_To_Set &states_to_update,const Original_State_t& original_states,MAT::Node* starting_node_parent_node,char starting_node_parent_state_one_hot) {
    // assert(to_fill==dst||(!to_fill->is_leaf()));

    size_t this_state_idx = get_state_idx(to_fill->index);
    auto &to_fill_state = states[this_state_idx];
#ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(to_fill_state.node == to_fill);
#endif

    if (to_fill_state >= 4) {
        MAT::Node *parent = get_parent(to_fill, edits);
        size_t parent_state_idx = get_state_idx(parent->index);
        const auto &parent_state = states[parent_state_idx];
#ifdef DETAIL_DEBUG_INDEX_MATCH
        assert(parent_state.node == parent);
#endif

        fill_dst_states(parent, score_offset, state_offset, mutation, scores,
                        states, dst, edits, states_to_update,original_states,starting_node_parent_node,starting_node_parent_state_one_hot);
#ifdef DETAIL_DEBUG_FS_CORRECT
        assert((1 << parent_state) ==
               get_genotype_with_regrafting(parent, mutation, edits,starting_node_parent_node,starting_node_parent_state_one_hot));
#endif
        if (to_fill->is_leaf() && (to_fill != dst)) {
            to_fill_state = one_hot_to_two_bit(
                get_original_state(to_fill, original_states, mutation));
        } else {
            size_t this_score_idx = get_score_idx(to_fill->index);
            const auto &this_score = scores[this_score_idx];
#ifdef DETAIL_DEBUG_INDEX_MATCH
            assert(this_score.node == to_fill);
#endif
            to_fill_state = Fitch_Sankoff::get_child_score_on_par_nuc(parent_state, this_score).second;
        }
        set_mutation(to_fill, 1 << to_fill_state, 1 << parent_state, mutation,
                     states_to_update);
    }
}

void Fitch_Sankoff::sankoff_forward_pass(const std::pair<size_t, size_t> &range,
                          std::vector<MAT::Node *> &dfs_ordered_nodes,const MAT::Mutation &mutation,const Original_State_t& original_states,
                          Scores_Type& scores, char starting_node_parent_state_one_hot,MAT::Node* to_move,MAT::Node* dst,MAT::Node* new_leaf,States_To_Set& states_to_update,const all_moves_bare_type& edits
) {
    #ifndef NDEBUG
    if (!states_to_update.empty()) {
        for(const auto& temp:states_to_update){
            assert(temp.first->index>=range.first&&temp.first->index<range.second);
        }
    }
    #endif
    // Score vector runs backward
    //first node->last element
    //index=size-1-(dfs_index-first_element_index)
    size_t score_offset=-1+range.second;
//#ifdef DETAIL_DEBUG_FS_CORRECT
    MAT::Node* starting_node_parent_node=get_parent(dfs_ordered_nodes[range.first],edits);
//#endif
    //States in 2bit format , original_states vector in one-hot format goes forward
    #ifdef DETAIL_DEBUG_INDEX_MATCH
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
    //in 2 bit format needed by get_score_by_par_nuc
    char starting_node_parent_state_2_bit=one_hot_to_two_bit(starting_node_parent_state_one_hot);

    #ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(states[0].node==dfs_ordered_nodes[range.first]);
    #endif

    Score_Type& last_score=scores[range.second-range.first-1];
    states[0]=get_child_score_on_par_nuc(starting_node_parent_state_2_bit, last_score).second;
#ifdef DETAIL_DEBUG_INDEX_MATCH
    assert(last_score.node==dfs_ordered_nodes[range.first]);
#endif

    set_mutation(dfs_ordered_nodes[range.first], 1<<states[0], 1<<starting_node_parent_state_2_bit, mutation,states_to_update);

    for (size_t dfs_idx=range.first+1; dfs_idx<range.second; dfs_idx++) {
        //Index conversions
        auto this_node=dfs_ordered_nodes[dfs_idx];
        size_t score_idx=get_score_idx(dfs_idx);

#ifdef DETAIL_DEBUG_INDEX_MATCH
        assert(this_node==scores[score_idx].node);
#endif

        MAT::Node* parent_node=get_parent(this_node, edits);
        size_t parent_state_idx=get_state_idx(parent_node->index);
#ifdef DETAIL_DEBUG_INDEX_MATCH
        assert(states[parent_state_idx].node==parent_node);
#endif

        char parent_state;
        
        fill_dst_states(parent_node, score_offset, state_offset, mutation, scores, states, dst,edits,states_to_update,original_states,starting_node_parent_node,starting_node_parent_state_one_hot);
#ifdef DETAIL_DEBUG_FS_CORRECT
        assert(get_genotype_with_regrafting(parent_node, mutation,edits,starting_node_parent_node,starting_node_parent_state_one_hot)==(1<<states[parent_state_idx]));
#endif
         parent_state=states[parent_state_idx];

        size_t this_state_idx=get_state_idx(dfs_idx);

#ifdef DETAIL_DEBUG_INDEX_MATCH
        assert(states[this_state_idx].node==this_node);
#endif

        //If to_move is moved to dst which is a leaf node, it will become an internal node, so treat it as an internal node
        if(this_node->is_leaf()&&this_node!=dst) {
            char original_state=get_original_state(this_node, original_states, mutation);
            set_mutation(this_node, original_state, 1<<parent_state, mutation,states_to_update);
    #ifdef DETAIL_DEBUG_FS_CORRECT
            assert(get_genotype_with_regrafting(this_node, mutation,edits,starting_node_parent_node,starting_node_parent_state_one_hot)==original_state);
            #endif
            continue;
        }


        char this_state=get_child_score_on_par_nuc(parent_state, scores[score_idx]).second;

        states[this_state_idx]=this_state;
        set_mutation(this_node, 1<<this_state, 1<<parent_state, mutation,states_to_update);
#ifdef DETAIL_DEBUG_FS_CORRECT
        assert(get_genotype_with_regrafting(this_node, mutation,edits,starting_node_parent_node,starting_node_parent_state_one_hot)==(1<<states[this_state_idx]));
#endif
        //deal with leaf node turned into internal node
        if(this_node->is_leaf()&&this_node==dst) {
            char original_state=get_original_state(this_node, original_states, mutation);
            set_mutation(new_leaf, original_state, 1<<this_state, mutation,states_to_update);
    #ifdef DETAIL_DEBUG_FS_CORRECT
            assert(get_genotype_with_regrafting(new_leaf, mutation,edits,starting_node_parent_node,starting_node_parent_state_one_hot)==original_state);
            #endif
        }
    }
}
