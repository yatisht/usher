#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include "Twice_Bloom_Filter.hpp"
void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output){
    Twice_Bloom_Filter filter;
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            filter.insert(m.position);
        }
    }
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            if(filter.query(m.position)){
                output.push_back(n);
                break;
            }
        }
    }
}
