#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
#include <algorithm>
#include <string>
#ifndef NDEBUG
#include "../mutation_annotated_tree.hpp"
namespace MAT = Mutation_Annotated_Tree;
#include "../Fitch_Sankoff.hpp"
#include "../check_samples.hpp"
#include <unordered_set>

extern std::unordered_map<MAT::Mutation,
       std::unordered_map<std::string, nuc_one_hot> *,
       Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
       mutated_positions;

struct Mutation_Count {
    int position;
    int count;
    std::unordered_map<std::string, std::vector<MAT::Node *>> mut_count;
    Mutation_Count() {}
    Mutation_Count(int position) : position(position), count(0) {}
    // Mutation_Count(const MAT::Mutation& mut,int
    // count):position(mut.get_position()),count(count){}
};

static void count_mutations(MAT::Node *start,
                            std::vector<Mutation_Count> &mutations_count) {
    size_t mut_idx = 0;
    for (const MAT::Mutation &m : start->mutations) {
        while (mut_idx != mutations_count.size() &&
                mutations_count[mut_idx].position < m.get_position()) {
            mut_idx++;
        }
        if (mut_idx != mutations_count.size() &&
                mutations_count[mut_idx].position == m.get_position() &&
                m.is_valid()) {
            mutations_count[mut_idx].count++;
            auto res = mutations_count[mut_idx].mut_count.emplace(
                           start->parent->identifier, 0);
            res.first->second.push_back(start);
        }
    }
    for (auto child : start->children) {
        count_mutations(child, mutations_count);
    }
}

static void remove_child(MAT::Node *child_to_remove) {
    auto parent = child_to_remove->parent;
    auto iter = std::find(parent->children.begin(), parent->children.end(),
                          child_to_remove);
    assert(iter != parent->children.end());
    parent->children.erase(iter);
}
MAT::Node *get_mutation_path(MAT::Mutations_Collection &mutations,
                             MAT::Node *src, MAT::Node *dst,
                             MAT::Node * &src_branch_node) {
    std::unordered_set<MAT::Node *> src_to_root;
    std::vector<MAT::Node *> src_to_root_path;
    assert(src != dst);
    mutations.reserve(src->mutations.size());
    for (const auto &mut : src->mutations) {
        if (mut.get_par_one_hot() != mut.get_all_major_allele()) {
            mutations.push_back(mut);
        }
    }
    MAT::Node *LCA = src->parent;
    while (LCA) {
        src_to_root.insert(LCA);
        src_to_root_path.push_back(LCA);
        LCA = LCA->parent;
    }

    LCA = dst;
    std::vector<MAT::Node *> dst_to_root_path;
    while (!src_to_root.count(LCA)) {
        dst_to_root_path.push_back(LCA);
        LCA = LCA->parent;
    }
    assert(dst_to_root_path.empty() || dst_to_root_path.back()->parent == LCA);

    mutations.clear();
    for (MAT::Node *n : src_to_root_path) {
        if (n == LCA) {
            break;
        }
        src_branch_node=n;
        mutations.merge(n->mutations, MAT::Mutations_Collection::KEEP_SELF);
    }
    for (auto iter = dst_to_root_path.rbegin(); iter < dst_to_root_path.rend();
            iter++) {
        mutations.merge((*iter)->mutations,
                        MAT::Mutations_Collection::KEEP_OTHER);
    }

    return LCA;
}
int get_parsimmony_score_only(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,const MAT::Tree* ori_tree) {
    MAT::Tree new_tree;
    MAT::Mutations_Collection mutations;
    MAT::Node * src_branch_node;
    MAT::Node *actual_LCA = get_mutation_path(mutations, src, dst, src_branch_node);
    mutations.merge(actual_LCA->mutations,
                    MAT::Mutations_Collection::KEEP_SELF);
    while (actual_LCA != LCA) {
        actual_LCA = actual_LCA->parent;
    }
    MAT::Node *new_LCA =
        new Mutation_Annotated_Tree::Node(*LCA, nullptr, &new_tree, false);
    new_tree.curr_internal_node = ori_tree->curr_internal_node;
    new_tree.root = new_LCA;
    MAT::Node *new_src = new_tree.get_node(src->identifier);
    std::vector<Mutation_Count> original_mutation_count;
    original_mutation_count.reserve(mutations.size());
    for (const auto &mut : mutations) {
        original_mutation_count.emplace_back(mut.get_position());
    }
    count_mutations(LCA, original_mutation_count);
    // MAT::Node *new_LCA = new_tree.get_node(LCA->identifier);
    MAT::Node *original_src_parent_in_new_tree = new_src->parent;
    MAT::Node *new_dst = new_tree.get_node(dst->identifier);
    remove_child(new_src);
    if (new_dst->identifier!=actual_LCA->identifier) {
        auto &dst_parent_children = new_dst->parent->children;
        dst_parent_children.erase(std::find(
                                      dst_parent_children.begin(), dst_parent_children.end(), new_dst));
        auto split_node = new_tree.create_node(
                              std::to_string(++new_tree.curr_internal_node), new_dst->parent);
        split_node->children.push_back(new_dst);
        new_dst->parent = split_node;
        split_node->children.push_back(new_src);
        new_src->parent = split_node;
    } else {
        auto& LCA_children=new_dst->children;
        auto new_src_branch_node=new_tree.get_node(src_branch_node->identifier);
        LCA_children.erase(std::find(LCA_children.begin(),LCA_children.end(),new_src_branch_node));
        auto split_node = new_tree.create_node(
                              std::to_string(++new_tree.curr_internal_node), new_dst);
        split_node->children.push_back(new_src_branch_node);
        new_src_branch_node->parent=split_node;
        split_node->children.push_back(new_src);
        new_src->parent = split_node;
    }
    std::vector<MAT::Node *> new_bfs_ordered_nodes =
        new_tree.breadth_first_expansion();
    int parsimony_score_change = 0;
    for (size_t idx = 0; idx < mutations.size(); idx++) {
        int position = mutations[idx].get_position();
        nuc_one_hot LCA_parent_state = get_parent_state(LCA, position);
        std::vector<uint8_t> boundary1_major_allele(
            new_bfs_ordered_nodes.size() + 8);
        MAT::Mutation mut(position);
        const auto& non_ref_muts = mutated_positions[mut];
        FS_backward_pass(new_bfs_ordered_nodes, boundary1_major_allele,
                         *non_ref_muts, MAT::Mutation::refs[position]);
        std::vector<uint8_t> states_out(new_bfs_ordered_nodes.size());
        std::vector<std::vector<MAT::Node *>> mutation_count(
                                               new_bfs_ordered_nodes.size());
        int new_parsimony_score = FS_forward_assign_states_only(
                                      new_bfs_ordered_nodes, boundary1_major_allele, LCA_parent_state,
                                      states_out, mutation_count);
        int this_change =
            new_parsimony_score - original_mutation_count[idx].count;
        if (this_change) {
            for (auto &temp : original_mutation_count[idx].mut_count) {
                MAT::Node *corresponding_node_in_new_tree =
                    new_tree.get_node(temp.first);
            }
        }
        parsimony_score_change += this_change;
    }
    new_tree.delete_nodes();
    return parsimony_score_change;
}
#endif
#endif