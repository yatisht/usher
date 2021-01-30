#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <queue>
#include <tuple>
#include <utility>
static MAT::Mutation reverse_mutation(MAT::Mutation to_reverse) {
    auto old_par_nuc = to_reverse.par_nuc;
    to_reverse.par_nuc = to_reverse.mut_nuc;
    to_reverse.mut_nuc = old_par_nuc;
    return to_reverse;
}


static MAT::Mutations_Collection* BFS(MAT::Node* src,MAT::Node* excluded, int radius,std::vector<Dst_Mut>& out){
    struct queue_content{
        MAT::Node* node;
        MAT::Node* reached_from; //take advantage of the tree, no loops
        int dist;
        MAT::Mutations_Collection& mutations;
    };
    MAT::Mutations_Collection* to_search=new MAT::Mutations_Collection;
    std::queue<queue_content> bfs_queue;
    MAT::Mutations_Collection empty;
    bfs_queue.push({src,excluded,radius,empty});
    
#define bfs_add_node(node,flag) \
MAT::Mutations_Collection mutation_path;\
mutations.merge_out(node->mutations, mutation_path,flag);\
out.push_back({node,mutation_path});\
bfs_queue.push({node,src,dist,mutation_path});\
to_search->merge(node->mutations, 1);

    while (!bfs_queue.empty()) {
        MAT::Node* src=bfs_queue.front().node;
        MAT::Node* excluded=bfs_queue.front().reached_from;
        int dist=bfs_queue.front().dist-1;
        MAT::Mutations_Collection& mutations=bfs_queue.front().mutations;
        bfs_queue.pop();
        if (dist<0) {
            break;
        }
        if (src->parent!=excluded&&src->parent&&src->parent->parent) {
            bfs_add_node(src->parent,MAT::Mutations_Collection::INVERT_MERGE)
        }
        for(MAT::Node* c:src->children){
            if(c==excluded) continue;
            bfs_add_node(c, MAT::Mutations_Collection::MERGE);
        }
    }
    
}
Possible_Moves* Neighbors_Finder::operator()(MAT::Node* src)const{
    Possible_Moves* result=new Possible_Moves;
    result->src=src;
    result->to_search=BFS(src->parent,src,radius,result->dst);
    
    tbb::concurrent_hash_map<MAT::Node*, tbb::concurrent_vector<int>>::const_accessor temp;
    bool search_pos=repeatedly_mutating_loci.find(temp,src);
    if(search_pos){
        for(int pos:temp->second){
            if(result->to_search->find(pos)!=result->to_search->end()){
                postponed.push_back(src);
                return nullptr;
            }
        }
    }
    
    return result;
}