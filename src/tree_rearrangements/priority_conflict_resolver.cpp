#include "priority_conflict_resolver.hpp"
#include "src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstdio>
#include <utility>
#include <vector>

bool Conflict_Resolver::check_single_move_no_conflict(Profitable_Moves_ptr_t& candidate_move)const{
    int best_score=0;
    for(auto node:candidate_move->path){
        auto iter=potential_crosses.find(node);
        if (iter!=potential_crosses.end()) {
            best_score=std::min(best_score,iter->second.parsimony_score);
        }
    }
    if (candidate_move->score_change<best_score) {
        return true;
    }
    return false;
}

static void remove_move(Cross_t &potential_crosses, Profitable_Moves_ptr_t& other_move,
                        MAT::Node *exclude) {
    for (MAT::Node *other_nodes_in_path : other_move->path) {
        if (other_nodes_in_path != exclude) {
            auto other_node_iter = potential_crosses.find(other_nodes_in_path);
            assert(other_node_iter != potential_crosses.end());
            assert(other_node_iter->second.moves.size() == 1 &&
               other_node_iter->second.moves[0] == other_move);
            potential_crosses.erase(other_node_iter);
        }
    }
}

void Conflict_Resolver::register_single_move_no_conflict(
    Profitable_Moves_ptr_t& candidate_move) const {
    for (auto node : candidate_move->path) {
        auto iter = potential_crosses.find(node);
        if (iter != potential_crosses.end()) {
            for (Profitable_Moves_ptr_t& other_move : iter->second.moves) {
                assert(other_move->score_change >=
                       candidate_move->score_change);
                remove_move(potential_crosses, other_move, node);
                #ifndef NDEBUG
                nodes_inside--;
                #endif
                #ifdef CHECK_LEAK
                other_move->destructed=true;
                #else
                delete other_move;
                #endif
            }
            iter->second.parsimony_score = candidate_move->score_change;
            iter->second.moves.clear();

            iter->second.moves.push_back(candidate_move);
        } else {
            potential_crosses.insert(std::make_pair(
                node,
                Conflict_Set{candidate_move->score_change,
                             std::vector<Profitable_Moves_ptr_t>{candidate_move}}));
        }
    }
}

char Conflict_Resolver::operator()(Profitable_Moves_From_One_Source* candidate_move)const{
    if (candidate_move->profitable_moves.empty()) {
        delete candidate_move;
        return 0;
    }
    char ret=MOVE_FOUND_MASK;
    Profitable_Moves_ptr_t selected_move=nullptr;
    for (Profitable_Moves_ptr_t& move : candidate_move->profitable_moves) {
        if (check_single_move_no_conflict(move)) {
            register_single_move_no_conflict(move);
            ret |= NONE_CONFLICT_MOVE_MASK;
            selected_move = move;
            break;
        }
    }
    if (!selected_move) {
        assert(ret==MOVE_FOUND_MASK);
        deferred_nodes.push_back(candidate_move->src);
    }
    #ifndef NDEBUG
    if(selected_move){
        nodes_inside++;
    }
    #endif
    for (Profitable_Moves_ptr_t move : candidate_move->profitable_moves) {
        if (move!=selected_move) {
            #ifdef CHECK_LEAK
            move->destructed=true;
            #else
            delete move;
            #endif
        }
    }

    delete candidate_move;
    return ret;
}

struct Move_Hash{
    size_t operator()(const Profitable_Moves_ptr_t& move) const{
        size_t src_idx=move->src->index;
        size_t dst_idx= move->dst->index;
        return src_idx ^( dst_idx + 0x9e3779b9 + (src_idx << 6) + (src_idx >> 2));
    }
};
struct Move_Eq{
    bool operator()(const Profitable_Moves_ptr_t& lhs,const Profitable_Moves_ptr_t& rhs)const{
        return (lhs->src==rhs->src)&&(lhs->dst==rhs->dst);
    }
};
struct Dependency_Node;
struct Dependency_Hash{
    size_t operator()(const Dependency_Node* in)const;
};
struct Dependency_Eq{
    bool operator()(const Dependency_Node* lhs,const Dependency_Node* rhs)const;
};

struct Dependency_Node {
    mutable int incoming_count;
    Profitable_Moves_ptr_t move;
    std::unordered_set<Dependency_Node *, Dependency_Hash, Dependency_Eq>
        children;
    Dependency_Node(const Profitable_Moves_ptr_t &move)
        : incoming_count(0), move(move) {}
    void add_reverse_dependency(Dependency_Node *dependency,
                                Fitch_Sankoff_Result_Final &this_state) {
        if (children.insert(dependency).second) {
            dependency->incoming_count++;
        }
        MAT::Node *lca_parent = dependency->move->LCA->parent;
        if (lca_parent) {
            assert(lca_parent->index >= move->range.first &&
                   lca_parent->index < move->range.second);
            for (Fitch_Sankoff_Result_Final &dependent_state :
                 dependency->move->states) {
                if (dependent_state.mutation.position ==
                    this_state.mutation.position) {
                    this_state.other_move_LCA_parent_states_to_update.insert(
                        std::make_pair(lca_parent,
                                       &dependent_state.LCA_parent_state));
                    return;
                }
            }
            assert(false);
        }
    }
};
size_t Dependency_Hash::operator()(const Dependency_Node* in)const{
    return Move_Hash()(in->move);
}
bool Dependency_Eq::operator()(const Dependency_Node* lhs,const Dependency_Node* rhs)const{
    return Move_Eq()(lhs->move,rhs->move);
}

typedef std::unordered_map<int,std::vector<Dependency_Node*>> Mut_Move_Map_t;
typedef std::unordered_map<MAT::Node*, Mut_Move_Map_t,Node_Idx_Hash,Node_Idx_Eq> Mut_t;

//typedef std::unordered_map<Profitable_Move*, Dependency_Node,Move_Hash,Move_Eq> Dependency_t;
typedef std::vector<Dependency_Node*> Dependency_t;
static void register_dependency(Profitable_Moves_ptr_t& move, Mut_t&repeated_mutations, Dependency_t& dependencies){
    // assert(repeated_mutations.find(move->LCA)==repeated_mutations.end());
    Dependency_Node *dep = new Dependency_Node(move);
    dependencies.push_back(dep);
    for (Fitch_Sankoff_Result_Final &e : move->states) {
        // for LCA only to check what it needs to update
        auto LCA_iter =
            repeated_mutations.insert({move->LCA, Mut_Move_Map_t()});
        Mut_Move_Map_t &mutation_to_insert = LCA_iter.first->second;
        auto LCA_mutation_insert_iter =
            mutation_to_insert.find(e.mutation.position);
        if (LCA_mutation_insert_iter != mutation_to_insert.end()) {
            for (Dependency_Node *reverse_dependency :
                 LCA_mutation_insert_iter->second) {
                dep->add_reverse_dependency(reverse_dependency, e);
            }
        }
        for(const auto& move_in_subtree:mutation_to_insert ){
            for (Dependency_Node* individual_move:move_in_subtree.second){
                move->other_moves_in_subtree.insert(std::make_pair(individual_move->move->src, individual_move->move->dst));
            }
        }
        auto LCA_ancestor_node = move->LCA->parent;
        while (LCA_ancestor_node) {
            auto LCA_ancestor_node_iter =
                repeated_mutations.insert({LCA_ancestor_node, Mut_Move_Map_t()})
                    .first;
            auto LCA_ancestor_node_mutation_result =
                LCA_ancestor_node_iter->second.insert(
                    std::make_pair(e.mutation.position,
                                   std::vector<Dependency_Node *>({dep})));
            if (!LCA_ancestor_node_mutation_result.second) {
                LCA_ancestor_node_mutation_result.first->second.push_back(dep);
            }
            LCA_ancestor_node = LCA_ancestor_node->parent;
        }
    }
}

struct Move_LCA_Comparator{
    bool operator()(const Profitable_Moves_ptr_t& lhs,const Profitable_Moves_ptr_t& rhs)const {
        return lhs->LCA->index>rhs->LCA->index;
    }
};

struct Move_LCA_Equal{
    bool operator()(const Profitable_Moves_ptr_t& lhs,const Profitable_Moves_ptr_t& rhs)const {
        if(lhs->LCA==rhs->LCA){
            assert(lhs==rhs);
            return true;
        }
        return false;
    }
};

static void order_moves(Cross_t& found_moves,std::vector<Profitable_Moves_ptr_t>& out){
    for(const auto& found_move:found_moves){
        out.push_back(found_move.second.moves.front());
    }
    std::sort(out.begin(),out.end(),Move_LCA_Comparator());

    auto end=std::unique(out.begin(),out.end(),Move_LCA_Equal());
    out.erase(end,out.end());
};
#ifndef NDEBUG
static void check_all_dep(Profitable_Moves_ptr_t& this_move,std::vector<Profitable_Moves_ptr_t>& other_moves,Mut_t& repeated_mutations){
    const std::pair<size_t, size_t> this_range=this_move->range;
    for(Profitable_Moves_ptr_t other_move:other_moves ){
        auto other_LCA_parent=other_move->LCA->parent;
        if(!other_LCA_parent) continue;
        size_t other_LCA_parent_idx=other_LCA_parent->index;
        if (other_LCA_parent_idx>=this_range.first&&other_LCA_parent_idx<this_range.second) {
            for(const auto& this_mutation:this_move->states){
                for(const auto& other_mutation:other_move->states){
                    if (this_mutation.mutation.position==other_mutation.mutation.position) {
                        auto state_to_set_iter=this_mutation.other_move_LCA_parent_states_to_update.find(other_LCA_parent);
                        if(state_to_set_iter==this_mutation.other_move_LCA_parent_states_to_update.end()){
                            fprintf(stderr, "mutation at %d of LCA parent with index %zu not counted, LCA idx %zu, this move LCA idx %zu \n",this_mutation.mutation.position,other_LCA_parent_idx,other_move->LCA->index,this_move->LCA->index);
                            auto other_LCA_parent_rep_node_iter=repeated_mutations.find(other_LCA_parent);
                            if (other_LCA_parent_rep_node_iter==repeated_mutations.end()) {
                                fprintf(stderr, "node not found in repeated mutation map\n");
                                assert(false);
                            }else {
                                auto other_LCA_parent_mut_iter=other_LCA_parent_rep_node_iter->second.find(this_mutation.mutation.position);
                                if(other_LCA_parent_mut_iter==other_LCA_parent_rep_node_iter->second.end()){
                                fprintf(stderr, "node found but mutation not found in repeated mutation map\n");
                                }
                                auto dep=other_LCA_parent_mut_iter->second;
                                for(auto& tmp:dep){
                                    if (tmp->move==other_move) {
                                        fprintf(stderr, "move found\n");
                                assert(false);
                                    }
                                }
                                assert(false);
                            }
                        }
                        assert(state_to_set_iter->second==&other_mutation.LCA_parent_state);
                    }
                }
            }
        }
    }
}
#endif
static void get_all_deps(Cross_t& found_moves,Dependency_t& dependencies){
    Mut_t repeated_mutations;
    std::vector<Profitable_Moves_ptr_t> ordered_moves;
    order_moves(found_moves,ordered_moves);
    for(auto& found_move:ordered_moves){
        register_dependency(found_move, repeated_mutations, dependencies);
        #ifndef NDEBUG
        check_all_dep(found_move, ordered_moves,repeated_mutations);
        #endif
    }
    assert(dependencies.size()==ordered_moves.size());
}

static void order_moves(Dependency_t& dependencies,std::vector<Profitable_Moves_ptr_t>& out){
    Dependency_t next_round;
    #ifndef NDEBUG
    bool have_progress=false;
    #endif
    while (!dependencies.empty()) {
        next_round.clear();
        for(auto move:dependencies){
            if (move->incoming_count) {
                next_round.push_back(move);
            }else {
                #ifndef NDEBUG
                have_progress=true;
                #endif
                for(auto reverse_dependency:move->children){
                    reverse_dependency->incoming_count--;
                }
                out.push_back(move->move);
                delete move;
            }
        }
        assert(have_progress);
        dependencies.swap(next_round);
    }
}

void schedule_moves(Cross_t& found_moves, std::vector<Profitable_Moves_ptr_t>& out,int& nodes_inside){
    Dependency_t dependencies;
    get_all_deps(found_moves, dependencies);
    assert(dependencies.size()==nodes_inside);
    nodes_inside=0;
    size_t dependency_size=dependencies.size();
    order_moves(dependencies, out);
    assert(dependency_size==out.size());
}