#include "../mutation_annotated_tree.hpp"
#include <array>
#include <unordered_map>
#include <utility>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
char get_genotype(MAT::Node* node, const Mutation_Annotated_Tree::Mutation& m);
namespace Fitch_Sankoff {
/*#ifndef NDEBUG
struct State_Type {
    char state;
    const MAT::Node *const node;
    State_Type(const MAT::Node *n) : node(n) {
        state=0;
    }
    State_Type() : State_Type(nullptr) {}
    operator char &() { return state; }
    char operator=(char a) { return (state = a); }
};
struct Score_Type {
    std::array<int, 4> score;
    const MAT::Node *const node;
    Score_Type(const MAT::Node *n) : node(n) {}
    Score_Type() : Score_Type(nullptr) { score[0] = -1; }
    int &operator[](size_t a) { return score[a]; }
};
#else*/
typedef char State_Type;
typedef std::array<int, 4> Score_Type;
//#endif
typedef std::vector<State_Type> States_Type;
typedef std::vector<Score_Type> Scores_Type;

std::pair<size_t, size_t> dfs_range(const MAT::Node *start,std::vector<MAT::Node *> &dfs_ordered_nodes);

int sankoff_backward_pass(const std::pair<size_t, size_t> &range,
                           const MAT::Mutation &mutation,
                           const std::vector<MAT::Node *> &dfs_ordered_nodes,
                           Scores_Type &scores, States_Type &states,std::vector<char> original_state);
void sankoff_forward_pass(const std::pair<size_t, size_t> &range,
                          States_Type &states,
                          std::vector<MAT::Node *> &dfs_ordered_nodes,const MAT::Mutation &mutation,std::vector<char> original_state
//#ifndef NDEBUG
//,Scores_Type &scores 
//#endif
);
std::pair<int, char>
get_child_score_on_par_nuc(char par_nuc,
                           Score_Type &child_scores);
} // namespace Fitch_Sankoff