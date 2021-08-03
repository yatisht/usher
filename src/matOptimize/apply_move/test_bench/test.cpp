#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "src/matOptimize/priority_conflict_resolver.hpp"
#include <cstddef>
#include <cstdio>
#include "tbb/parallel_for_each.h"
#include <random>
#include <chrono>
#include <vector>
std::unordered_map<MAT::Mutation,
    std::unordered_map<std::string, nuc_one_hot> *,
    Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator>
    mutated_positions;
static Profitable_Moves_ptr_t make_move(MAT::Node* src,MAT::Node* dst) {
    std::vector<MAT::Node*> src_to_LCA;
    std::vector<MAT::Node*> dst_to_LCA;
    while (src!=dst) {
        if (src->dfs_index>dst->dfs_index) {
            src_to_LCA.push_back(src);
            src=src->parent;
        } else if (src->dfs_index<dst->dfs_index) {
            dst_to_LCA.push_back(dst);
            dst=dst->parent;
        }
    }
    dst_to_LCA.push_back(src);
    Profitable_Moves_ptr_t ret=new Profitable_Moves;
    ret->dst_to_LCA=std::move(dst_to_LCA);
    ret->src_to_LCA=std::move(src_to_LCA);
    ret->LCA=src;
    ret->src=ret->src_to_LCA[0];
    return ret;
}
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_PRIMARY_MOVE
                 ,
                 const Original_State_t& original_state
#endif
                ) ;
int main(int argc, char** argv) {
    Original_State_t origin_states;
    MAT::Tree t;
    t.load_detatiled_mutations(argv[1]);
    fprintf(stderr, "%zu \n",t.get_parsimony_score());
    auto bfs_ordered_nodes = t.breadth_first_expansion();
    check_samples(t.root, origin_states, &t);
    /*for (MAT::Node *node : bfs_ordered_nodes) {
        for (const MAT::Mutation &m : node->mutations) {
            mutated_positions.emplace(
                m, new std::unordered_map<std::string, nuc_one_hot>);
        }
        node->tree = &t;
    }
        tbb::parallel_for_each(
        mutated_positions.begin(), mutated_positions.end(),
        [&origin_states](
            const std::pair<MAT::Mutation,
                            std::unordered_map<std::string, nuc_one_hot> *>
                &pos) {
            std::unordered_map<std::string, nuc_one_hot> *mutated = pos.second;
            for (auto &sample : origin_states) {
                auto iter = sample.second.find(pos.first);
                if (iter != sample.second.end()) {
                    mutated->emplace(sample.first, iter->get_all_major_allele());
                }
            }
        });*/
    std::vector<Profitable_Moves_ptr_t> all_moves{};

    tbb::concurrent_vector<MAT::Node *> deferred_nodes;
    if (argc==3) {
        FILE* moves=fopen(argv[2], "r");
        char src[BUFSIZ];
        char dst[BUFSIZ];
        Conflict_Resolver resolver(bfs_ordered_nodes.size(),fopen("/dev/null", "w"));
        while (fscanf(moves, "Trying %s to %s\n",src,dst)!=EOF) {
            if (std::string(src)=="43419") {
                fputc('a',stderr);
            }
            MAT::Node* src_node=t.get_node(src);
            MAT::Node* dst_node=t.get_node(dst);
            if(!(src_node&&dst_node)) {
                continue;
            }
            Profitable_Moves_ptr_t move=make_move(src_node,dst_node);
            //all_moves.push_back(move);
            move->score_change=-1;
            std::vector<Profitable_Moves_ptr_t>temp{move};
            resolver(temp);
        }
        resolver.schedule_moves(all_moves);
        apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes,origin_states);
        fprintf(stderr, "%zu \n",t.get_parsimony_score());
        return 0;
    }

    /*
        //From https://www.cplusplus.com/reference/random/mersenne_twister_engine/mersenne_twister_engine/
        //unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
        unsigned seed1 = 234;
        std::mt19937 g1 (seed1);
        while(true){
        all_moves.clear();
        std::vector<MAT::Node*> dfs_ordered_nodes=t.depth_first_expansion();
        t.breadth_first_expansion();
        std::uniform_int_distribution<int> src_idx_dist(1,dfs_ordered_nodes.size()-1);
        std::uniform_int_distribution<int> score_change_idx_dist(-5,-1);
        Conflict_Resolver resolver(dfs_ordered_nodes.size());
        for(int i=0;i<1;i++){
            size_t src_idx=src_idx_dist(g1);
            MAT::Node* src_node=dfs_ordered_nodes[src_idx];
            size_t dst_range=dfs_ordered_nodes.size()-src_node->dfs_end_index+src_idx;
            size_t dst_idx_raw=g1()%dst_range;
            if (dst_idx_raw>=src_idx) {
                dst_idx_raw+=(src_node->dfs_end_index-src_idx);
            }
            if(dst_idx_raw==src_node->parent->dfs_index){
                i--;
                continue;
            }
            Profitable_Moves_ptr_t move=make_move(src_node, dfs_ordered_nodes[dst_idx_raw]);
            move->score_change=score_change_idx_dist(g1);
            std::vector<Profitable_Moves_ptr_t>temp{move};
            resolver(temp);
        }
        resolver.schedule_moves(all_moves);
        apply_moves(all_moves, t, bfs_ordered_nodes, deferred_nodes,origin_states);
        Mutation_Annotated_Tree::save_mutation_annotated_tree(t, "last_tree.pb");}
        */
}