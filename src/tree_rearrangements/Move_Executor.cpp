#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/Fitch_Sankoff.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>
MAT::Node* Move_Executor::get_parent(MAT::Node* node) const{
        auto iter=new_parents_map.find(node);
        if(iter==new_parents_map.end()){
            return node->parent;
        }else{
            return (MAT::Node*)iter->second;
        }
}
static size_t get_parsimony_score(MAT::Node* root){
    int par_score=root->mutations.size();
    for(auto child:root->children){
        par_score+=get_parsimony_score(child);
    }
    return par_score;
}
int count_mutation(MAT::Node* node,int pos){
    unsigned int count=0;
    if (node) {
        count=std::count_if(node->mutations.begin(), node->mutations.end(),[pos](const MAT::Mutation& m){return m.position==pos;} );
    }
    assert(count<=1);
    return count;
}
/*
#ifndef NDEBUG
void check_original_states(Fitch_Sankoff::States_Type& original_states, Mutation_Annotated_Tree::Mutation& mutation,const Original_State_t& mut){
    for(Fitch_Sankoff::State_Type& state:original_states){
        auto iter=mut.find(state.node->identifier);
        if (iter!=mut.end()) {
            assert(state.node->is_leaf());
            auto mut_iter=iter->second.find(mutation);
            if(mut_iter!=iter->second.end()){
                assert(state.state==mut_iter->mut_nuc);
            }else {
                assert(state.state==mutation.ref_nuc);
            }
        }else{
            assert(!state.node->is_leaf());
        }
    }
}
#endif
*/
struct Node_To_Edit {
    MAT::Node *this_Node;
    ConfirmedMove edits;
    const all_moves_bare_type* all_other_moves;
    void operator++() {
        edits.added.clear();
        edits.removed.clear();
        this_Node = get_parent(this_Node,*all_other_moves);
    }
    operator bool() { return this_Node; }
    // Max Heap, want the child come out before parent
    bool operator<(const Node_To_Edit &other) const {
        return this_Node->index < other.this_Node->index;
    }

};

static void moves_bare_type_to_node_to_edit(const all_moves_bare_type &in,
                                            std::vector<Node_To_Edit> &out,MAT::Node* this_src, MAT::Node* this_dst) {
    Pending_Moves_t temp;
    for (const auto &move : in) {
        if (move.first->is_leaf() || move.second->is_leaf()) {
            // State of leaf node won't change, and it was assumed to stay
            // constant in profitable move enumerator, not a concern.
            continue;
        }
        ConfirmedMove ttt;
        auto op_node = temp.insert(std::make_pair(move.first->parent, ttt));
        op_node.first->second.removed.push_back(move.first);
        assert(op_node.first->first == move.first->parent);
        op_node = temp.insert(std::make_pair(move.second, ttt));
        op_node.first->second.added.push_back(move.first);
    }
    for (const auto &node_with_children_changed : temp) {
        out.push_back(Node_To_Edit{node_with_children_changed.first,
                       node_with_children_changed.second,&in});
    }
    std::make_heap(out.begin(), out.end());
}
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
}
static void proc_child(const MAT::Node *child, size_t start_idx,
                       const Fitch_Sankoff::Scores_Type &out,
                       std::vector<char> &min_score_possible_nuc,
                       std::vector<int> &min_scores) {
    size_t child_idx = start_idx - child->index;
    /* make sure indexes match*/
    assert(child_idx <= out.size() - 2);
    assert(out[child_idx].node == child);
    /* Make sure child values are computed */
    assert(out[child_idx][0] != -1);
    int min_score;
    char arg_min_one_hot;
    find_min_4(out[child_idx], min_score, arg_min_one_hot);
    min_score_possible_nuc.push_back(arg_min_one_hot);
    min_scores.push_back(min_score);
}

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
    assert(min_score_possible_nuc.size() == final_child_size);
    assert(min_scores.size() == final_child_size);
    bool score_changed = false;
    size_t this_node_idx = start_idx - this_node.index;
    Fitch_Sankoff::Score_Type &this_node_score = out[this_node_idx];
    assert(this_node_score.node == &this_node);
    for (char par_nuc = 0; par_nuc < 4; par_nuc++) {
        int new_score = 0;
        char par_nuc_one_hot = 1 << par_nuc;
        for (size_t child_idx = 0; child_idx < final_child_size; child_idx++) {
            new_score += min_scores[child_idx];
            if (par_nuc_one_hot & min_score_possible_nuc[child_idx]) {
                new_score++;
            }
        }
        if (this_node_score[par_nuc] != new_score) {
            score_changed = true;
            this_node_score[par_nuc] = new_score;
        }
    }
    return score_changed;
}

static void patch_score(std::pair<size_t, size_t> &range,
                        Fitch_Sankoff_Result_Final &to_patch,
                        const all_moves_bare_type &in,MAT::Node* this_src, MAT::Node* this_dst) {
    std::vector<Node_To_Edit> edits_in_subtree;
    moves_bare_type_to_node_to_edit(in, edits_in_subtree,this_src,this_dst);
    int offset = range.second - 1;
    while (!edits_in_subtree.empty()) {
        std::pop_heap(edits_in_subtree.begin(), edits_in_subtree.end());
        while (edits_in_subtree.front().this_Node==edits_in_subtree.back().this_Node) {
            auto& to_remove_content=edits_in_subtree.back().edits;
            if (!to_remove_content.empty()) {
                assert(edits_in_subtree.front().edits.empty());
                edits_in_subtree.front().edits.swap(to_remove_content);
            }
            edits_in_subtree.pop_back();
        }
        auto &this_changed_node = edits_in_subtree.back();
        if (!patch_internal_score(*this_changed_node.this_Node, to_patch.scores,
                                  offset, this_changed_node.edits)) {
            edits_in_subtree.pop_back();
        } else {
            ++this_changed_node;
            if ((!this_changed_node)||this_changed_node.this_Node->index<range.first) {
                edits_in_subtree.pop_back();
            } else {
                std::push_heap(edits_in_subtree.begin(),
                               edits_in_subtree.end());
            }
        }
    }
}
static void remove_child(MAT::Node *child_to_remove) {
    auto parent = child_to_remove->parent;
    auto iter = std::find(parent->children.begin(), parent->children.end(),
                          child_to_remove);
    assert(iter != parent->children.end());
    parent->children.erase(iter);
}

#ifdef DETAIL_DEBUG_MOVE_EXE_SCORE_RECALC
struct New_Temp_Tree {
    MAT::Tree new_tree;
    MAT::Tree *old_tree;
    std::pair<size_t, size_t> new_range;
    std::vector<MAT::Node *> new_dfs_ordered_nodes;
    const Original_State_t &original_states;
    MAT::Node *added;
    MAT::Node *new_dst;
    MAT::Node* old_dst;
    New_Temp_Tree(MAT::Node *LCA, const all_moves_bare_type &other_moves_to_apply,
                  const Original_State_t &original_states,MAT::Node* this_src,MAT::Node* this_dst)
        : original_states(original_states) {
            old_dst=this_dst;
        MAT::Node *new_LCA =
            new Mutation_Annotated_Tree::Node(*LCA, nullptr, &new_tree);
        old_tree = LCA->tree;
        new_tree.root = new_LCA;
        for (const auto &move : other_moves_to_apply) {
            if (move.first==this_src) {
                assert(move.second==this_dst);
                continue;
            }
            MAT::Node *new_src = new_tree.get_node(move.first->identifier);
            // MAT::Node *new_LCA = new_tree.get_node(LCA->identifier);
            MAT::Node *new_dst = new_tree.get_node(move.second->identifier);
            remove_child(new_src);
            assert(new_dst->add_child(new_src) == nullptr);
        }
        MAT::Node *new_src = new_tree.get_node(this_src->identifier);
        // MAT::Node *new_LCA = new_tree.get_node(LCA->identifier);
        MAT::Node *original_src_parent_in_new_tree = new_src->parent;
        new_dst = new_tree.get_node(this_dst->identifier);
        remove_child(new_src);
        added = new_dst->add_child(new_src);
        while (original_src_parent_in_new_tree->is_leaf()) {
            MAT::Node *parent_of_parent =
                original_src_parent_in_new_tree->parent;
            remove_child(original_src_parent_in_new_tree);
            original_src_parent_in_new_tree = parent_of_parent;
        }
        new_dfs_ordered_nodes = new_tree.depth_first_expansion();
        new_range = Fitch_Sankoff::dfs_range(new_LCA, new_dfs_ordered_nodes);
    }
    void check_score(const Fitch_Sankoff::Scores_Type &patched,
                     const MAT::Mutation &mutation) {
                        std::unordered_map<std::string, std::array<int, 4>>
                            scores_to_check;
                        for (auto &&score : patched) {
                            if (added && score.node == old_dst) {
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
                            original_states, mutation);
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
                                        old_tree
                                            ->get_node(score.node->identifier)
                                            ->index);
                                assert(false);
                            }
                        }
    }
};
#endif
void Move_Executor::operator()(tbb::blocked_range<size_t> &r) const {
    for (auto i = r.begin(); i < r.end(); i++) {
        Profitable_Move *this_move = moves[i];
        all_moves_bare_type moves_in_subtree;
        for (const auto &move : this_move->other_moves_in_subtree) {
            if (!(move.first->is_leaf() || move.second->is_leaf())) {
                // State of leaf node won't change, and it was assumed to stay
                // constant in profitable move enumerator, not a concern.
                moves_in_subtree.insert(move);
            }
        }
        MAT::Node* new_leaf=nullptr;
        moves_in_subtree.emplace(this_move->src, this_move->dst);
        /*if (!moves_in_subtree.empty()) {
            #ifdef DETAIL_DEBUG_MOVE_EXE_SCORE_RECALC
            New_Temp_Tree moved_tree(this_move->LCA,moves_in_subtree,ori,this_move->src,this_move->dst);
            #endif
            
            for (Fitch_Sankoff_Result_Final &this_mutation_states :this_move->states) {
                patch_score(this_move->range, this_mutation_states, this_move->other_moves_in_subtree,this_move->src,this_move->dst);
                moved_tree.check_score(this_mutation_states.scores, this_mutation_states.mutation);
            }
        }
#ifndef NDEBUG
        for (Fitch_Sankoff_Result_Final &this_mutation_states :
             this_move->states) {
            
                auto parent = get_parent(this_move->LCA);
                while (parent) {
                    auto mutation_iter = parent->mutations.find(
                        this_mutation_states.mutation.position);
                    if (mutation_iter != parent->mutations.end()) {
                        assert(this_mutation_states.LCA_parent_state ==
                               mutation_iter->mut_nuc);
                        goto END;
                    }
                    parent=get_parent(parent);
                }
                assert(this_mutation_states.LCA_parent_state ==
                       this_mutation_states.mutation.ref_nuc);
            END:;
            
        }
#endif
            /* Apply move no matter whether it reduce parsimony score as a simplification
            if (this_mutation_states.LCA_parent_state_original_tree !=
                this_mutation_states.LCA_parent_state) {
                int new_score = Fitch_Sankoff::get_child_score_on_par_nuc(
                                    one_hot_to_two_bit(this_mutation_states.LCA_parent_state),
                                    this_mutation_states.scores.back())
                                    .first;
                assert(new_score <= this_mutation_states.optimized_score + 1);
                this_move->score_change +=
                    (new_score - this_mutation_states.optimized_score);
            }
        }
        if(this_move->score_change>=0){
            continue;
        }* /
#ifndef NDEBUG
        for(auto m:tree_edits){
            const auto& other_removed=m.second.removed;
            assert(std::find(other_removed.begin(),other_removed.end(),this_move->src)==other_removed.end());
        }
#endif
*/
        // Register Move
        ConfirmedMove temp;
        auto op_node = tree_edits.insert(std::make_pair(this_move->src->parent, temp));
        op_node.first->second.removed.push_back(this_move->src);
        assert(op_node.first->first==this_move->src->parent);
        op_node = tree_edits.insert(std::make_pair(this_move->dst, temp));
        std::vector<MAT::Node*>& to_add=op_node.first->second.added;
        if(this_move->dst->is_leaf()){
            if(to_add.empty()){
                new_leaf=new MAT::Node();
                new_leaf->parent=this_move->dst;
                new_leaf->identifier=this_move->dst->identifier;
                to_add.push_back(new_leaf);
            }else{
                new_leaf=to_add.front();
                assert(new_leaf->identifier==this_move->dst->identifier);
            }
        }
        to_add.push_back(this_move->src);

        assert(new_parents_map.insert(std::make_pair(this_move->src, this_move->dst)).second);
        for(const auto& move:this_move->other_moves_in_subtree){
            if (!(move.first->is_leaf() || move.second->is_leaf()) ){
                // State of leaf node won't change, and it was assumed to stay
                // constant in profitable move enumerator, not a concern.
                moves_in_subtree.insert(move);
            }
        }
        //moves_in_subtree.emplace(this_move->src,this_move->dst);
        // conflicts checked at the last step
        for (size_t state_idx=0;state_idx<this_move->states.size();state_idx++) {
            Fitch_Sankoff_Result_Final& m=this_move->states[state_idx];
            Fitch_Sankoff::sankoff_forward_pass(
                this_move->range, dfs_ordered_nodes, m.mutation, ori,
                m.scores, m.LCA_parent_state, this_move->src, this_move->dst,new_leaf,m.other_move_LCA_parent_states_to_update,moves_in_subtree);
            /*
            #ifndef NDEBUG
            int this_mut_score=count_mutation(new_leaf, m.mutation.position);
            for (size_t node_idx=this_move->range.first; node_idx<this_move->range.second; node_idx++) {
                this_mut_score+=count_mutation(dfs_ordered_nodes[node_idx], m.mutation.position);
            }
            assert(this_mut_score==m.optimized_score);
            #endif
            */
        }
/*
#ifndef NDEBUG
        Original_State_t copy(ori);
        Mutation_Set parental;
        MAT::Node* ancestor=get_parent(this_move->LCA);
        while(ancestor){
            for(const auto& m:ancestor->mutations){
                    parental.insert(m);
            }
            ancestor=get_parent(ancestor);
        }
        for (auto iter = parental.begin(); iter != parental.end();) {
            if (iter->mut_nuc == iter->ref_nuc) {
                auto prev_iter = iter;
                iter++;
                parental.erase(prev_iter);
            } else {
                iter++;
            }
        }
        check_samples_worker_with_pending_moves(this_move->LCA, parental, copy,tree_edits);
#endif
*/
#ifdef NDEBUG
        delete this_move;
#endif
    }
}
