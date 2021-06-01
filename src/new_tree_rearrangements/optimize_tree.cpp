#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "priority_conflict_resolver.hpp"
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst){
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src=src->parent;
        }
        else if (src->dfs_index<dst->dfs_index) {
            dst=dst->parent;
        }
    }
    return src;
}

bool check_not_ancestor(MAT::Node* dst,MAT::Node* src){
    while (dst) {
        if (dst==src) {
            return false;
        }
        dst=dst->parent;
    }
    return true;
}

size_t optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
              tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
              MAT::Tree &t,int radius
              #ifndef NDEBUG
              , Original_State_t origin_states
              #ifdef CONFLICT_RESOLVER_DEBUG
              ,FILE* log
            #endif
            #endif
              ) {
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", bfs_ordered_nodes.size());
    fprintf(stderr, "Internal node size %zu\n", t.curr_internal_node);
    for(auto node:bfs_ordered_nodes){
        node->changed=false;
    }
    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    Deferred_Move_t deferred_moves;
    Conflict_Resolver resolver(bfs_ordered_nodes.size(),deferred_moves
#ifdef CONFLICT_RESOLVER_DEBUG
    ,log
#endif
    );
    output_t out;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes_to_search.size()),
                      [&nodes_to_search, &resolver,
                       &deferred_nodes,radius](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              output_t out;
                              find_profitable_moves(nodes_to_search[i], out, radius);
                              if (!out.moves.empty()) {
                                  deferred_nodes.push_back(
                                      out.moves[0]->get_src());
                                  resolver(out.moves);
                              }
                          }
                      });
    std::vector<Profitable_Moves_ptr_t> all_moves;
    resolver.schedule_moves(all_moves);
    apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes
#ifdef CHECK_STATE_REASSIGN
    ,origin_states
#endif
    );
    while (!deferred_moves.empty()) {
        bfs_ordered_nodes=t.breadth_first_expansion();
        {Deferred_Move_t deferred_moves_next;
        Conflict_Resolver resolver(bfs_ordered_nodes.size(),deferred_moves_next
#ifdef CONFLICT_RESOLVER_DEBUG
    ,log
#endif
    );
        tbb::parallel_for(tbb::blocked_range<size_t>(0,deferred_moves.size()),[&deferred_moves,&resolver,&t](const tbb::blocked_range<size_t>& r){
            for (size_t i=r.begin(); i<r.end(); i++) {
                MAT::Node* src=t.get_node(deferred_moves[i].first);
                MAT::Node* dst=t.get_node(deferred_moves[i].second);
                if (src&&dst) {
                    output_t out;
                        if (check_not_ancestor(dst, src)) {
                            individual_move(src,dst,get_LCA(src, dst),out);
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
    ,origin_states
#endif
    );
        all_moves.clear();
        deferred_moves=std::move(deferred_moves_next);}
    
    }
    #ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
    #endif
    nodes_to_search = std::move(deferred_nodes);
    return t.get_parsimony_score();
}
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
void save_final_tree(MAT::Tree &t, Original_State_t& origin_states,
                            const std::string &output_path) {
#ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
#endif
    std::vector<MAT::Node *> dfs = t.depth_first_expansion();
#ifndef NDEBUG
    check_samples(t.root, origin_states, &t);
#endif
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
                      [&dfs](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                                  dfs[i]->mutations.remove_invalid();
                          }
                      });
    fix_condensed_nodes(&t);
    fprintf(stderr, "%zu condensed_nodes",t.condensed_nodes.size());
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, output_path);
}