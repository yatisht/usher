#include "priority_conflict_resolver.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tbb/parallel_for_each.h"
#include "tree_rearrangement_internal.hpp"
#include <cstdio>
#include <string>
#include <utility>
#include <vector>
MAT::Node *get_LCA(MAT::Node *src, MAT::Node *dst);
bool check_not_ancestor(MAT::Node *dst, MAT::Node *src);
int main(int argc, char **argv) {
    MAT::Tree t;
    t.load_detatiled_mutations(
        "/scratch/home/cheng/usher/testout/2021-01-20/2021-01-20-msa-ambiguious-tree-patch-itermediate0.pb");
    auto bfs_ordered_nodes = t.breadth_first_expansion();
    for (auto node : bfs_ordered_nodes) {
        node->tree = &t;
    }
    Original_State_t origin_states;
    check_samples(t.root, origin_states, &t);
    // save_final_tree(t,origin_states,"tttttt");
    for (MAT::Node *node : bfs_ordered_nodes) {
        for (const MAT::Mutation &m : node->mutations) {
            mutated_positions.emplace(
                m, new std::unordered_map<std::string, nuc_one_hot>);
        }
        node->tree = &t;
    }
    /*tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&origin_states](
            const std::pair<MAT::Mutation,
                            std::unordered_map<std::string, nuc_one_hot> *>
                &pos) {
            std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
            for (auto &sample : origin_states) {
                auto iter = sample.second.find(pos.first);
                if (iter != sample.second.end()) {
                    mutated->emplace(sample.first,
                                     iter->get_all_major_allele());
                }
            }
        });*/

    FILE *log = fopen("testout/try_move_debug", "w");
    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    Deferred_Move_t deferred_moves;
    Conflict_Resolver resolver(bfs_ordered_nodes.size(), deferred_moves
#ifdef CONFLICT_RESOLVER_DEBUG
                               ,
                               log
#endif
    );
    std::vector<Profitable_Moves_ptr_t> all_moves{};
    fprintf(stderr, "%zu \n", t.get_parsimony_score());
    output_t out;
    tbb::concurrent_vector<MAT::Node *> nodes_to_search;
    FILE *moves = fopen(argv[2], "r");
    char src[BUFSIZ];
    char dst[BUFSIZ];
    while (fscanf(moves, "Trying %s to %s\n", src, dst) != EOF) {
        MAT::Node *src_node = t.get_node(src);
        MAT::Node *dst_node = t.get_node(dst);
        if (!(src_node && dst_node)) {
            continue;
        }
        output_t out;
        individual_move(src_node, dst_node, get_LCA(src_node, dst_node), out);
        if (!out.moves.empty()) {
            resolver(out.moves);
        }
        // deferred_moves.push_back(std::make_pair(src,
        // std::vector<std::string>{dst}));
    }
    resolver.schedule_moves(all_moves);
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
                ,
                origin_states
#endif
    );
    all_moves.clear();

    while (!deferred_moves.empty()) {
        bfs_ordered_nodes=t.breadth_first_expansion();
        Deferred_Move_t deferred_moves_next;
        Conflict_Resolver resolver(bfs_ordered_nodes.size(), deferred_moves_next
#ifdef CONFLICT_RESOLVER_DEBUG
                                   ,
                                   log
#endif
        );
        tbb::parallel_for(tbb::blocked_range<size_t>(0,deferred_moves.size()),[&deferred_moves,&resolver,&t](const
         tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i <r.end(); i++) {
            MAT::Node *src = t.get_node(deferred_moves[i].first);
            if (src) {
                output_t out;
                for (auto dst_id : deferred_moves[i].second) {
                    auto dst = t.get_node(dst_id);
                    if (dst && check_not_ancestor(dst, src)) {
                        individual_move(src, dst, get_LCA(src, dst), out);
                    }
                }
                if (!out.moves.empty()) {
                    resolver(out.moves);
                }
            }
        }
        });
        resolver.schedule_moves(all_moves);
        apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
                    ,
                    origin_states
#endif
        );
        all_moves.clear();
        deferred_moves = std::move(deferred_moves_next);
    }
    check_samples(t.root, origin_states, &t);
    return 0;
}
/*
change of -1 @ 24389
change of -1 @ 24390
change of 1 @ 28881
change of 1 @ 28882
change of 1 @ 28883
*/