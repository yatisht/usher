#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "process_each_node.hpp"
#include <cstdio>
#include <unordered_set>
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
static void check_LCA(MAT::Node* src, MAT::Node* dst,MAT::Node* LCA_to_check){
    std::unordered_set<int> path;
    while (src!=LCA_to_check) {
        path.insert(src->bfs_index);
        src=src->parent;
    }
    while (dst!=LCA_to_check) {
        assert(!path.count(dst->bfs_index));
        dst=dst->parent;
    }
    assert(dst==LCA_to_check);
}
#endif

bool dst_branch(const MAT::Node *LCA,
           const Mutation_Count_Change_Collection &mutations,
           int &parsimony_score_change,
           std::vector<MAT::Node *> &node_stack_from_dst, MAT::Node *this_node,
           Mutation_Count_Change_Collection &parent_added,bool& early_stop,int src_side_max_improvement
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
           ,std::vector<Mutation_Count_Change_Collection> &debug_from_dst
#endif
           );

void output_result(MAT::Node *&src, MAT::Node *&dst, MAT::Node *&LCA,
               int &parsimony_score_change, output_t &output,
               const std::vector<MAT::Node *> &node_stack_from_src,
               std::vector<MAT::Node *> &node_stack_from_dst,
               std::vector<MAT::Node *> &node_stack_above_LCA) {
    if (parsimony_score_change <= output.score_change) {
        if (parsimony_score_change < output.score_change) {
            output.score_change = parsimony_score_change;
            for (auto t : output.moves) {
                delete t;
            }
            output.moves.clear();
        }
        Profitable_Moves_ptr_t new_move = new Profitable_Moves;
        new_move->score_change = parsimony_score_change;
        new_move->src_to_LCA = std::move(node_stack_from_src);
        new_move->dst_to_LCA = std::move(node_stack_from_dst);
        new_move->src = src;
        new_move->LCA = LCA;
        assert(dst == LCA || new_move->dst_to_LCA.back()->parent == LCA);
        assert(src->parent == LCA ||
               new_move->src_to_LCA.back()->parent == LCA);
        MAT::Node* last_node=nullptr;
        for (const auto node : node_stack_above_LCA) {
            if (node!=last_node) {
                new_move->dst_to_LCA.push_back(node);
                last_node=node;
            }
        }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        std::unordered_set<int> node_idx_set;
        node_idx_set.insert(src->bfs_index);
        for (auto node : new_move->src_to_LCA) {
            assert(node_idx_set.insert(node->bfs_index).second);
        }
        for (auto node : new_move->dst_to_LCA) {
            assert(node_idx_set.insert(node->bfs_index).second);
        }
#endif
        output.moves.push_back(new_move);
    }
}
int get_parsimmony_score_only(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,const MAT::Tree* ori_tree);
MAT::Node *check_move_profitable_LCA(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int& parsimony_score_change,
    const std::vector<MAT::Node *> &node_stack_from_dst,
    Mutation_Count_Change_Collection &parent_added,
    const std::vector<MAT::Node *>& node_stack_from_src,
    std::vector<MAT::Node *>& node_stack_above_LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,
    std::vector<Mutation_Count_Change_Collection>& debug_above_LCA
#endif
) ;
extern unsigned int early_stop_saving;
/**
 * @brief Calculate Parsimony score change of this individual move without applying it
 * @param src
 * @param dst
 * @param LCA
 * @param mutations mutations needed above src node for moving src under parent of dst
 * @param root_mutations_altered fitch set change of src_branch_node from removing src
 * @param parsimony_score_change
 * @param node_stack_from_src nodes from src to LCA
 * @param debug_above_LCA
 * @param output Most profitable moves from this src node
 * @return The highest node had their Fitch set changed
 */
int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    int parsimony_score_change, output_t &output,
    const std::vector<MAT::Node *> node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<Mutation_Count_Change_Collection> debug_from_src,
    const MAT::Tree* tree
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    check_LCA(src, dst, LCA);
    assert(src->parent!=dst);
    assert(src!=dst);
    std::vector<Mutation_Count_Change_Collection> debug_from_dst;
#endif
    std::vector<MAT::Node *> node_stack_from_dst({});

    assert(dst);
    bool early_stop=false;
    Mutation_Count_Change_Collection dst_added;
    //Going up from dst node to LCA node to adjust state assignment
    if (LCA != dst) {
        if(!dst_branch(LCA, mutations, parsimony_score_change,
                  node_stack_from_dst, dst, dst_added,early_stop,src->parent==LCA?src->mutations.size():0
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
, debug_from_dst
#endif
                  )){
                      return 1;
        }
    }

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<Mutation_Count_Change_Collection> debug_above_LCA;
#endif
    std::vector<MAT::Node *> node_stack_above_LCA;
    //Adjust LCA node and above
    MAT::Node* ancestor=check_move_profitable_LCA(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change, node_stack_from_dst, dst_added, node_stack_from_src, node_stack_above_LCA
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
     ,debug_above_LCA
#endif
    );
    if(!ancestor){
        return 1;
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    //fprintf(stderr, "LCA idx: %zu",LCA_a->bfs_index);
    auto ref_score=get_parsimmony_score_only(src,dst,LCA,tree);
    assert(parsimony_score_change == ref_score);
#endif
    assert((dst==LCA&&node_stack_above_LCA[0]==LCA)||node_stack_from_dst[0]==dst);
    if (early_stop) {
        if (parsimony_score_change<0) {
            fprintf(stderr, "early stop will miss %s to %s of change %d\n",src->identifier.c_str(),dst->identifier.c_str(),parsimony_score_change);
        }else {
            early_stop_saving++;
        }
    }
    output_result(src, dst, LCA, parsimony_score_change, output,
              node_stack_from_src, node_stack_from_dst, node_stack_above_LCA);
    return parsimony_score_change;
}