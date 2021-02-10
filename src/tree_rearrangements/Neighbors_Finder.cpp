#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <cstddef>
#include <queue>
#include <tuple>
#include <utility>

static void BFS(MAT::Node* src,MAT::Node* ori_src, int radius,tbb::flow::interface11::internal::multifunction_output<struct Possible_Move *>& out){
    struct queue_content{
        MAT::Node* node;
        MAT::Node* reached_from; //take advantage of the tree, no loops
        int dist;
    };
    std::queue<queue_content> bfs_queue;
    MAT::Node* excluded=ori_src;
    bfs_queue.push({src,excluded,radius});
    
#define bfs_add_node(node) \
out.try_put(new Possible_Move{src,node});\
bfs_queue.push({node,src,dist});\

    while (!bfs_queue.empty()) {
        MAT::Node* src=bfs_queue.front().node;
        MAT::Node* excluded=bfs_queue.front().reached_from;
        int dist=bfs_queue.front().dist-1;
        bfs_queue.pop();
        if (dist<0) {
            break;
        }
        if (src->parent!=excluded&&src->parent&&src->parent->parent) {
            bfs_add_node(src->parent)
        }
        for(MAT::Node* c:src->children){
            if(c==excluded) continue;
            bfs_add_node(c);
        }
    }
    
}
void Neighbors_Finder::operator()(MAT::Node* src, Neighbors_Finder_t::output_ports_type& out)const{
    BFS(src->parent,src,radius,std::get<0>(out));
}