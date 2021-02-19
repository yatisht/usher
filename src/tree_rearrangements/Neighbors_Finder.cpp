#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <queue>
#include <tuple>
#include <utility>

static void BFS(MAT::Node* src,MAT::Node* ori_src, int radius,std::vector<MAT::Node*>& out){
    struct queue_content{
        MAT::Node* node;
        MAT::Node* reached_from; //take advantage of the tree, no loops
        int dist;
    };
    std::queue<queue_content> bfs_queue;
    MAT::Node* excluded=ori_src;
    bfs_queue.push({src,excluded,radius});

#define bfs_add_node(node) \
out.emplace_back(node);\
bfs_queue.push({node,src,dist});\

    while (!bfs_queue.empty()) {
        MAT::Node* src=bfs_queue.front().node;
        MAT::Node* excluded=bfs_queue.front().reached_from;
        int dist=bfs_queue.front().dist-1;
        bfs_queue.pop();
        if (dist<0) {
            break;
        }
        if (src->parent&&src->parent->parent&&src->parent!=excluded) {
            bfs_add_node(src->parent)
        }
        for(MAT::Node* c:src->children){
            if(c==excluded) continue;
            bfs_add_node(c);
        }
    }

}
Possible_Moves* Neighbors_Finder::operator()(MAT::Node* src)const{
    Possible_Moves* result=new Possible_Moves;
    result->src=src;
    BFS(src->parent,src,radius,result->dsts);
    return result;
}