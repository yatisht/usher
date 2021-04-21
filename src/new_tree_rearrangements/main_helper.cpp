#include "mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "Fitch_Sankoff.hpp"
#include "src/new_tree_rearrangements/Twice_Bloom_Filter.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
namespace MAT=Mutation_Annotated_Tree;
extern std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
static void add_root(MAT::Tree *tree) {
    MAT::Node *old_root = tree->root;
    MAT::Node *new_root = new MAT::Node();
    new_root->tree=tree;
    new_root->identifier = std::to_string(++tree->curr_internal_node);
    std::string &node_name = new_root->identifier;
    new_root->children.push_back(old_root);
    tree->all_nodes.emplace(node_name, new_root);
    old_root->parent=new_root;
    tree->root=new_root;
}
Mutation_Annotated_Tree::Tree load_tree(char* path,Original_State_t& origin_states){
    Mutation_Annotated_Tree::Tree t =
        Mutation_Annotated_Tree::load_mutation_annotated_tree(path);
    fprintf(stderr, "original score:%zu\n", t.get_parsimony_score());

    auto bfs_ordered_nodes = t.breadth_first_expansion();

    for (MAT::Node *node : bfs_ordered_nodes) {
        for (const MAT::Mutation &m : node->mutations) {
            mutated_positions.emplace(
                m, new std::unordered_map<std::string, nuc_one_hot>);
        }
        node->tree = &t;
    }

    
    check_samples(t.root, origin_states, &t);

    if (t.root->children.size()>1) {
        add_root(&t);
        bfs_ordered_nodes = t.breadth_first_expansion();
    }
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&origin_states, &bfs_ordered_nodes, &output](
            const std::pair<MAT::Mutation,
                            std::unordered_map<std::string, nuc_one_hot> *>
                &pos) {
            std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
            for (auto &sample : origin_states) {
                auto iter = sample.second.find(pos.first);
                if (iter != sample.second.end()) {
                    mutated->emplace(sample.first, iter->get_mut_one_hot());
                }
            }
            Fitch_Sankoff_Whole_Tree(bfs_ordered_nodes, pos.first, *mutated,
                                     output);
        });
    tbb::affinity_partitioner ap;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto &to_refill = output[i];
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size());
            }
        },
        ap);
    return t;
}
void fix_condensed_nodes(MAT::Tree *tree) {
    std::vector<MAT::Node *> nodes_to_fix;
    for (auto iter : tree->all_nodes) {
        if (tree->condensed_nodes.count(iter.first) &&
            (!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for (auto node : nodes_to_fix) {
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier,
                          std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier, node);
    }
}


static void find_nodes_with_recurrent_mutations(
    const std::vector<MAT::Node *> &all_nodes,
    std::vector<MAT::Node *> &output) {
    Twice_Bloom_Filter filter;
    for (MAT::Node *n : all_nodes) {
        for (const MAT::Mutation &m : n->mutations) {
            filter.insert(m.get_position());
        }
    }
    for (MAT::Node *n : all_nodes) {
        for (const MAT::Mutation &m : n->mutations) {
            if (filter.query(m.get_position())) {
                output.push_back(n);
                break;
            }
        }
    }
}
 void
find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                   tbb::concurrent_vector<MAT::Node*> &output) {
    std::vector<MAT::Node *> nodes_with_recurrent_mutations;
    find_nodes_with_recurrent_mutations(bfs_ordered_nodes,
                                        nodes_with_recurrent_mutations);

    std::unordered_set<size_t> pushed_nodes;
    for (MAT::Node *this_node : nodes_with_recurrent_mutations) {
        while (this_node->parent) {
            auto result = pushed_nodes.insert(this_node->bfs_index);
            if (result.second) {
                output.push_back(this_node);
                this_node = this_node->parent;
            } else {
                break;
            }
        }
    }
}
