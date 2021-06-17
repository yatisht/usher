#include "mutation_annotated_tree.hpp"
#include <vector>
#include "check_samples.hpp"
#pragma once
namespace MAT = Mutation_Annotated_Tree;
extern std::unordered_map<MAT::Mutation, std::unordered_map<std::string, nuc_one_hot>*,Mutation_Pos_Only_Hash,
                       Mutation_Pos_Only_Comparator>
        mutated_positions;
struct Profitable_Moves{
    int score_change;
    MAT::Node* src;
    std::vector<MAT::Node*> src_to_LCA;
    std::vector<MAT::Node*> dst_to_LCA;
    MAT::Node* LCA;
    bool new_node;
    Profitable_Moves():score_change(0){}
    Profitable_Moves(MAT::Node* src, MAT::Node* dst){
        score_change=-1;
        this->src=src;
        dst_to_LCA.push_back(dst);
    }
    MAT::Node* get_src()const{
        return src;
    }
    MAT::Node* get_dst()const{
        return dst_to_LCA.front();
    }
    template<typename F>
    void apply_nodes(F f){
        f(src);
        for(auto node: src_to_LCA){
            f(node);
        }
        for(auto node: dst_to_LCA){
            f(node);
        }
    } 
};
typedef Profitable_Moves* Profitable_Moves_ptr_t;
struct output_t{
    int score_change;
    std::vector<Profitable_Moves_ptr_t> moves;
    output_t():score_change(-1){}
};
void find_profitable_moves(Mutation_Annotated_Tree::Node *src, output_t &out,int radius
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,MAT::Tree* tree
#endif
);
int individual_move(Mutation_Annotated_Tree::Node* src,Mutation_Annotated_Tree::Node* dst,Mutation_Annotated_Tree::Node* LCA,output_t& out
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,MAT::Tree* tree
#endif
);
Mutation_Annotated_Tree::Tree load_tree(const std::string& path,Original_State_t& origin_states);
Mutation_Annotated_Tree::Tree load_vcf_nh_directly(const std::string& nh_path,const std::string& vcf_path,Original_State_t& origin_states);
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_STATE_REASSIGN
                 ,
                 const Original_State_t& original_state
#endif
);
void fix_condensed_nodes(MAT::Tree *tree) ;
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                   tbb::concurrent_vector<MAT::Node*> &output,bool is_first,int radius) ;
                   void add_root(MAT::Tree *tree) ;
void VCF_input(const char * name,MAT::Tree& tree);

size_t optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
              tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
              MAT::Tree &t,int radius,FILE* log
              #ifndef NDEBUG
              , Original_State_t origin_states
            #endif
              );
void save_final_tree(MAT::Tree &t, Original_State_t& origin_states,const std::string &output_path);
//For removing nodes with no valid mutations between rounds
void clean_tree(MAT::Tree& t);
void populate_mutated_pos(const Original_State_t& origin_state);