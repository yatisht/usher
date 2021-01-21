#include "Fitch_Sankoff.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement.hpp"
#include "usher_graph.hpp"
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "check_samples.hpp"
/* This is for non-adjacent moves
class Merged_Mutation_Iterator {
    class range {
        std::vector<MAT::Mutation>::const_iterator start;
        std::vector<MAT::Mutation>::const_iterator end;

      public:
        range(std::vector<MAT::Mutation> in)
            : start(in.begin()), end(in.end()) {}
        const std::vector<MAT::Mutation>::const_iterator &operator->() const {
            return start;
        }
        const std::vector<MAT::Mutation>::const_iterator &operator*() const {
            return start;
        }
        void operator++(int) { start++; }
        operator bool() { return start != end; }
    };
    struct comparator {
        bool operator()(const range &lhs, const range rhs) {
            return lhs->position > rhs->position;
        }
    };
    std::vector<range> heap;
    static std::vector<range>
    create_from_node_helper(const MAT::Node *parent, const MAT::Node *exclude) {
        std::vector<range> iterators;

        return iterators;
    }

  public:
    Merged_Mutation_Iterator(const MAT::Node *parent,
                             const MAT::Node *exclude) {
        heap.reserve(parent->children.size());
        for (MAT::Node *child : parent->children) {
            if (child == exclude)
                continue;
            if (child->mutations.size()) {
                heap.emplace_back(child->mutations);
            }
        }
        std::make_heap(heap.begin(), heap.end(), comparator());
    }
    void operator++(int) {
        std::pop_heap(heap.begin(), heap.end(), comparator());
        range &incremented_iter = heap.back();
        if (incremented_iter) {
            std::push_heap(heap.begin(), heap.end(), comparator());
        } else {
            heap.pop_back();
        }
    }
    const std::vector<MAT::Mutation>::const_iterator &operator->() const {
        return *(heap.front());
    }
    operator bool() { return !heap.empty(); }
};
*/
typedef std::vector<std::vector<char>> Original_States;
struct change_bundle{
    MAT::Node* new_parent;
    std::vector<Fitch_Sankoff::States_Type> states;
    Original_States original_states;
#ifndef NDEBUG
    std::vector<Fitch_Sankoff::Scores_Type> scores;
#endif

    void swap(change_bundle& other){
        new_parent=other.new_parent;
        states.swap(other.states);
        original_states.swap(other.original_states);
#ifndef NDEBUG
        scores.swap(other.scores);
#endif        
    }
};

static MAT::Mutation reverse_mutation(MAT::Mutation to_reverse) {
    auto old_par_nuc = to_reverse.par_nuc;
    to_reverse.par_nuc = to_reverse.mut_nuc;
    to_reverse.mut_nuc = old_par_nuc;
    return to_reverse;
}

static void distribute(MAT::Mutations_Collection& src, std::vector<std::vector<char>>& dst){
    assert(src.size()==dst.size());
    auto iter=src.begin();
    for(auto& a:dst){
        a.push_back(iter->mut_nuc);
        iter++;
    }
} 


Original_States get_original_states(std::vector<MAT::Node*> dfs_ordered_nodes, const std::pair<size_t, size_t> &range,MAT::Mutations_Collection target){
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
static void insert_node(MAT::Node *parent, MAT::Node *to_insert,
                        MAT::Tree &tree) {
    MAT::Mutations_Collection best_common;
    MAT::Mutations_Collection best_sibling_unique;
    MAT::Mutations_Collection best_insert_unique;
    MAT::Node *best_sibling;
    for (auto child : parent->children) {
        MAT::Mutations_Collection this_common;
        MAT::Mutations_Collection this_sibling_unique;
        MAT::Mutations_Collection this_insert_unique;
        child->mutations.set_difference(to_insert->mutations,
                                        this_sibling_unique, this_insert_unique,
                                        this_common);
        if (this_common.size() > best_common.size()) {
            best_sibling = child;
            best_common.swap(this_common);
            best_sibling_unique.swap(this_sibling_unique);
            best_insert_unique.swap(this_insert_unique);
        }
    }

    auto &parent_children = parent->children;
    if (best_common.size() == 0) {
        parent_children.push_back(to_insert);
    } else {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              best_sibling);
        parent_children.erase(iter);

        auto new_node = tree.create_node(
            std::to_string(++tree.curr_internal_node), parent->identifier);
        new_node->mutations.swap(best_common);

        to_insert->mutations.swap(best_insert_unique);
        new_node->children.push_back(to_insert);

        best_sibling->mutations.swap(best_sibling_unique);
        new_node->children.push_back(best_sibling);

        parent_children.push_back(new_node);
    }
}

template <bool move_to_parent>
static void
apply_move(MAT::Node *this_node, const std::pair<size_t, size_t> &range,
            std::vector<MAT::Node *> &dfs_ordered_nodes,MAT::Tree &tree,change_bundle& bundle
           ) {
    MAT::Node *new_parent=bundle.new_parent;
    std::vector<Fitch_Sankoff::States_Type>& states_all_pos=bundle.states;
    #ifndef NDEBUG
    std::vector<Fitch_Sankoff::Scores_Type>& scores_all_pos=bundle.scores;
    #endif
    Mutation_Annotated_Tree::Mutations_Collection mutations(move_to_parent?this_node->parent->mutations:new_parent->mutations);
    auto mutations_begin=mutations.begin();
    #ifndef NDEBUG
    auto mutations_end=mutations.end();
    Sample_Mut_Type ori;
    auto node_check=this_node->parent->parent;
    check_samples(node_check,ori);
    #endif
    // start moving
    MAT::Node* parent = this_node->parent;
    std::vector<MAT::Node*>& parent_children = parent->children;
    auto iter =
        std::find(parent_children.begin(), parent_children.end(), this_node);
    parent_children.erase(iter);
    if (parent_children.size() <= 1 && parent->not_sample()) {
        auto& grandpa_child=parent->parent->children;
        if(parent_children.size()==1){
            Mutation_Annotated_Tree::Mutations_Collection temp;
            auto child=parent_children[0];
            parent->mutations.merge_out(child->mutations, temp, 0);
            child->mutations.swap(temp);
            grandpa_child.push_back(child);
            child->parent=parent->parent;
        }
        grandpa_child.erase(std::find(grandpa_child.begin(),grandpa_child.end(),parent));

    }
    // Try merging with sibling in the new location
    insert_node(new_parent,this_node,tree);
    this_node->parent=new_parent;
    // apply adjustment to subtree
    #ifndef NDEBUG
    auto score_iter=scores_all_pos.begin();
    #endif
    assert(states_all_pos.size()==mutations_end-mutations_begin);
    std::unordered_map<MAT::Node*, MAT::Node*> ori_child;
    auto ori_state_iter=bundle.original_states.begin();
    for (auto states : states_all_pos) {
        assert(states.back().node==this_node);
        Fitch_Sankoff::sankoff_forward_pass(
            range, states, dfs_ordered_nodes, *mutations_begin,
            move_to_parent ? mutations_begin->par_nuc : mutations_begin->mut_nuc,*ori_state_iter,ori_child,tree
#ifndef NDEBUG
            ,*score_iter
#endif
        );
        mutations_begin++;
        ori_state_iter++;
    #ifndef NDEBUG
    score_iter++;
    #endif
    }
    assert(mutations_begin == mutations_end);
    assert(ori_state_iter==bundle.original_states.end());
    #ifndef NDEBUG
    check_samples(node_check,ori);
    #endif
}

template <bool move_to_parent>
static int check_move_profitable(
    MAT::Node *this_node, const std::pair<size_t, size_t> &range,
    MAT::Node *new_parent, const std::vector<MAT::Node *> &dfs_ordered_nodes,
    change_bundle& out
    ) {
    std::vector<Fitch_Sankoff::States_Type>& states_all_pos=out.states;
    #ifndef NDEBUG
    std::vector<Fitch_Sankoff::Scores_Type>& scores_all_pos=out.scores;
    #endif
    out.new_parent=new_parent;
    auto& mutations=move_to_parent?this_node->parent->mutations:new_parent->mutations;
    states_all_pos.reserve(mutations.size());
    int score_change = 0; // change in parsimony score if moved to new location
    auto original_genotype=get_original_states(dfs_ordered_nodes, range, mutations);
    auto genotype_iter=original_genotype.begin();
    for (auto&this_mutation: mutations) {
        Fitch_Sankoff::Scores_Type scores;
        Fitch_Sankoff::States_Type states;
        Fitch_Sankoff::sankoff_backward_pass(range, this_mutation,
                                             dfs_ordered_nodes, scores, states,*genotype_iter);

        char cur_parent_nuc = 31-__builtin_clz((unsigned int)(move_to_parent ? this_mutation.mut_nuc
                                             : this_mutation.par_nuc));
        char new_parent_nuc = 31-__builtin_clz((unsigned int)(move_to_parent ? this_mutation.par_nuc
                                             : this_mutation.mut_nuc));

        // score for current parent
        auto curr_result = Fitch_Sankoff::get_child_score_on_par_nuc(
            cur_parent_nuc, scores.back());
        auto new_result = Fitch_Sankoff::get_child_score_on_par_nuc(
            new_parent_nuc, scores.back());
        // always prepare states for new parent, as nothing will be done if
        // there is no improvement and this_node stay with its original parent.
        states.back() = new_result.second;
        score_change += (new_result.first - curr_result.first);
        states_all_pos.push_back(states);
        #ifndef NDEBUG
        scores_all_pos.push_back(scores);
        #endif
        genotype_iter++;
    }
    assert(states_all_pos.size()==mutations.size());
    if (score_change >= 0) {
        states_all_pos.clear();
        #ifndef NDEBUG
        scores_all_pos.clear();
        #endif
    }else{
        out.original_states.swap(original_genotype);
    }
    return score_change;
}

/**
 * @brief See whether there will be any improvement by moving current node as a
 * child of parent of parent of a child of its sibling. NNI can be achieved with
 * such moves:
 *       parent of parent
 *          /       \
 *         A        parent
 *                   /  \
 *                  B   this_node
 * Move this_node to parent of parent (parent and B fused):
 *      parent of parent
 *          /    |      \
 *         A     B      this_node
 * Make a new parent for A and B if they share mutations:
 *      parent of parent
 *         /      \
 *        *      this_node
 *       / \
 *      B   A
 * Finish a NNI of A and this_node
 * Moving under siblings unlikely be NNI
 *          parent
 *      /      |    \
 * this_node   B     sibling
 *                      \
 *                        A
 * Move this under sibling
 *    parent
 *  /       \
 *  B     sibling
 *        /    \
 *       A   this_node
 * @param this_node
 * @return bool Whether there is an improvement from moving this_node
 */

bool Tree_Rearrangement::move_nearest(
    MAT::Node *this_node, std::vector<MAT::Node *> &dfs_ordered_nodes,
    MAT::Tree &tree) {
    if (!this_node->parent || !this_node->parent->parent) {
        return false;
    }
    auto range = Fitch_Sankoff::dfs_range(this_node,dfs_ordered_nodes);
    // moving to parent
    change_bundle best;

    auto best_score = check_move_profitable<true>(
        this_node, range, this_node->parent->parent, dfs_ordered_nodes, best);

    for (auto sibling : this_node->parent->children) {
        if (sibling == this_node)
            continue;
        change_bundle this_bundle;
        auto this_score=check_move_profitable<false>(this_node, range, sibling,
                dfs_ordered_nodes,this_bundle) ;
        if (this_score< best_score) {
            best.swap(this_bundle);
            best_score=this_score;
        }
    }
    if (best_score >= 0) {
        return false;
    }
    if (best.new_parent == this_node->parent->parent) {
        apply_move<true>(this_node, range, dfs_ordered_nodes, tree, best);
    } else {
        apply_move<false>(this_node, range, dfs_ordered_nodes, tree, best);
    }
    return true;
}