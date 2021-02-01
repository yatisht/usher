#include "tree_rearrangement_internal.hpp"
bool check_mut_conflict(MAT::Node* node, const std::vector<Fitch_Sankoff_Result*>& states);
bool check_loop_conflict(MAT::Node* src, MAT::Node* dst);
void register_mut_conflict(MAT::Node* node, const std::vector<Fitch_Sankoff_Result*>& states);
void register_loop_conflict(MAT::Node* src, MAT::Node* dst);

struct multiple_moves_resolver{
    tbb::concurrent_vector<Profitable_Moves*>& multiple_optimal_deferred;
    mutable tbb::concurrent_vector<Profitable_Moves*>::iterator iter;
    mutable bool valid;
    multiple_moves_resolver(tbb::concurrent_vector<Profitable_Moves*>& multiple_optimal_deferred):multiple_optimal_deferred(multiple_optimal_deferred),iter(multiple_optimal_deferred.begin()),valid(true){}
    multiple_moves_resolver(const multiple_moves_resolver& other):multiple_optimal_deferred(other.multiple_optimal_deferred),iter(other.iter),valid(true){
        other.valid=false;
    }
    Profitable_Moves* operator() (tbb::flow_control) const;

};