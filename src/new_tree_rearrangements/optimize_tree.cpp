#include "tree_rearrangement_internal.hpp"
#include "priority_conflict_resolver.hpp"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
size_t optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
              tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
              MAT::Tree &t, Original_State_t origin_states,int radius
              #ifdef CONFLICT_RESOLVER_DEBUG
              ,FILE* log
            #endif
              ) {
    fprintf(stderr, "%zu nodes to search \n", nodes_to_search.size());
    fprintf(stderr, "Node size: %zu\n", bfs_ordered_nodes.size());
    fprintf(stderr, "Internal node size %zu\n", t.curr_internal_node);

    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    Conflict_Resolver resolver(bfs_ordered_nodes.size()
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
    check_samples(t.root, origin_states, &t);
    nodes_to_search = std::move(deferred_nodes);
    return t.get_parsimony_score();
}
std::unordered_map<MAT::Mutation,
                   std::unordered_map<std::string, nuc_one_hot> *,
                   Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
void save_final_tree(MAT::Tree &t, Original_State_t origin_states,
                            const std::string &output_path) {
    std::vector<MAT::Node *> dfs = t.depth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs.size()),
                      [&dfs](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              dfs[i]->mutations.remove_invalid();
                          }
                      });
    fix_condensed_nodes(&t);
    fprintf(stderr, "%d condensed_nodes",t.condensed_leaves.size());
    check_samples(t.root, origin_states, &t);
    Mutation_Annotated_Tree::save_mutation_annotated_tree(t, output_path);
}