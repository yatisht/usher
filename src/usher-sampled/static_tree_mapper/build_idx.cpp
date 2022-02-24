#include "index.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <vector>
#include <algorithm>
#include <array>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <stack>
#include <sys/types.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>
#include "src/usher-sampled/usher.hpp"
namespace MAT = Mutation_Annotated_Tree;
struct Temp_Idx_Tree_Node{
    MAT::Node* covering_node;
    int dfs_start_idx;
    int dfs_end_idx;
    bool is_self;
    int out_idx;
    std::vector<Temp_Idx_Tree_Node*> children;
};
typedef std::vector<Temp_Idx_Tree_Node*> Temp_Tree_Node_Coll_t;
struct temp_tree_build_comp {
    bool operator()(const Temp_Idx_Tree_Node* a,const Temp_Idx_Tree_Node* b)const {
        auto a_idx=a->covering_node->dfs_index;
        auto b_idx=b->covering_node->dfs_index;
        if (a_idx<b_idx) {
            return true;
        } else if (a_idx==b_idx&&b->is_self) {
            return true;
        }
        return false;
    }
};
static void check_node(const Temp_Idx_Tree_Node* to_check){
    //check sorting
    auto idx=to_check->covering_node->dfs_index;
    if (idx>to_check->dfs_start_idx) {
        fprintf(stderr, "Not covered\n");
        raise(SIGTRAP);
    }
    idx=to_check->dfs_start_idx;
    for (const auto child : to_check->children) {
        if (idx>child->dfs_start_idx) {
            fprintf(stderr, "Children out of order\n");
            raise(SIGTRAP);
        }
        idx=child->dfs_end_idx;
    }
    if (idx>to_check->dfs_end_idx) {
        fprintf(stderr, "children Not covered\n");
        raise(SIGTRAP);
    }
    idx=to_check->dfs_end_idx;
    if (idx>to_check->covering_node->dfs_end_index) {
        fprintf(stderr, "Not covered end\n");
        raise(SIGTRAP);
    }

}
static void set_idx(Temp_Idx_Tree_Node* to_set){
    to_set->dfs_start_idx=to_set->children.front()->dfs_start_idx;
    to_set->dfs_end_idx=to_set->children.back()->dfs_end_idx;
}
static void bulid_idx_tree(Temp_Tree_Node_Coll_t& in,std::vector<index_ele>& output){
    if (in.empty()) {
        output.push_back(index_ele{INT_MAX,INT_MAX,-1});
        return;
    }else if (in.size()==1) {
        output.push_back(index_ele{in[0]->dfs_start_idx,in[0]->dfs_end_idx,-1});
        output.push_back(index_ele{INT_MAX,INT_MAX,-1});
        return;
    }
    std::make_heap(in.begin(),in.end(),temp_tree_build_comp());
    int node_count=in.size();
    while (in.size()>1) {
        std::pop_heap(in.begin(),in.end(),temp_tree_build_comp());
        auto curr=in.back();
        in.pop_back();
        std::vector<Temp_Idx_Tree_Node*> children;
        while ((!in.empty())&&in.front()->covering_node==curr->covering_node) {
            auto front=in.front();
            children.push_back(front);        
            if (front->is_self) {
                if (front->dfs_start_idx==front->covering_node->dfs_index) {
                    fprintf(stderr, "Mult self node\n");
                    raise(SIGTRAP);
                }
            }
            std::pop_heap(in.begin(),in.end(),temp_tree_build_comp());
            in.pop_back();
        }
        auto par_node=curr->covering_node->parent;
        if (!par_node) {
            par_node=curr->covering_node;
        }
        bool is_curr_self=curr->is_self &&
            curr->dfs_start_idx == curr->covering_node->dfs_index;
        if (is_curr_self) {
            if (!curr->children.empty()) {
                Temp_Idx_Tree_Node *temp_node = new Temp_Idx_Tree_Node;
                temp_node->is_self = false;
                temp_node->children = std::move(curr->children);
                set_idx(temp_node);
                children.push_back(temp_node);
            }
        } else {
            children.push_back(curr);
        }
        std::vector<Temp_Idx_Tree_Node *> children_expanded;
        children_expanded.reserve(children.size());
        std::sort(children.begin(), children.end(),
                  [](Temp_Idx_Tree_Node *first, Temp_Idx_Tree_Node *second) {
                      return first->children.size() < second->children.size();
                  });
        int count = children.size();
        for (const auto& child : children) {
            if (child->is_self) {
                children_expanded.push_back(child);
            }else if (child->children.size()==1) {
                children_expanded.push_back(child->children[0]);
                delete child;
            }else if ((child->children.size()-1)+count<8) {
                children_expanded.insert(children_expanded.end(),child->children.begin(),child->children.end());
                count+=(child->children.size()-1);
                delete child;
            }else {
                children_expanded.push_back(child);
            }
        }
        std::sort(children_expanded.begin(),children_expanded.end(),[](Temp_Idx_Tree_Node* first,Temp_Idx_Tree_Node* second){
            return first->dfs_start_idx<second->dfs_start_idx;
        });
        if (is_curr_self) {
            curr->children=std::move(children_expanded);
            auto old_curr=curr;
            check_node(old_curr);
            curr=new Temp_Idx_Tree_Node{curr->covering_node->parent,old_curr->dfs_start_idx,old_curr->dfs_end_idx,false};
            node_count++;
            curr->children.push_back(old_curr);
        }else {
            curr=new Temp_Idx_Tree_Node;
            curr->children=std::move(children_expanded);
            curr->covering_node=par_node;
            curr->is_self=false;
            set_idx(curr);
        }
        check_node(curr);
        in.push_back(curr);
        std::push_heap(in.begin(),in.end(),temp_tree_build_comp());
    }
    if (in.empty()) {
        raise(SIGTRAP);
    }
    in.front()->out_idx=-1;
    std::vector<Temp_Idx_Tree_Node*> to_push{in.front()};
    std::vector<Temp_Idx_Tree_Node*> next_to_push;
    output.reserve(node_count);
    while (!to_push.empty()) {
        for (auto ele : to_push) {
            if (ele->out_idx != -1) {
                output[ele->out_idx].children_idx = output.size();
            }
            for (auto child : ele->children) {
                if ((!output.empty())&&output.back().dfs_idx_end!=INT_MAX&&output.back().dfs_idx_end>=child->dfs_start_idx) {
                    fprintf(stderr, "Output out of order\n");
                    raise(SIGTRAP);
                }
                output.push_back(index_ele{child->dfs_start_idx, child->dfs_end_idx,
                                    -1});
                if (child->children.empty()) {
                    delete child;
                } else {
                    child->out_idx = output.size() - 1;
                    next_to_push.push_back(child);
                }
            }
            delete ele;
        }
        output.emplace_back(index_ele{INT_MAX,INT_MAX,-1});
        to_push.clear();
        to_push.swap(next_to_push);
    }
    output.emplace_back(index_ele{INT_MAX,INT_MAX,-1});

}
Traversal_Info build_idx(MAT::Tree& tree){
    Traversal_Info out;
    auto dfs=tree.depth_first_expansion();
    out.traversal_track.reserve(dfs.size());
    std::vector<std::array<Temp_Tree_Node_Coll_t,4>> idx_node_vec(MAT::Mutation::refs.size());
    out.indexes.resize(MAT::Mutation::refs.size());
    for (auto node : dfs) {
        out.traversal_track.emplace_back(traversal_track_elem{(int)node->mutations.size(),(int)node->dfs_end_index,node->mutations.data()});
        for (const auto& mut : node->mutations) {
            auto new_temp_node=new Temp_Idx_Tree_Node{node,static_cast<int>(node->dfs_index),static_cast<int>(node->dfs_end_index),true};
            idx_node_vec[mut.get_position()][one_hot_to_two_bit(mut.get_mut_one_hot())].push_back(new_temp_node);
        }
    }
    tbb::parallel_for(tbb::blocked_range<size_t>(0,MAT::Mutation::refs.size()),[&idx_node_vec,&out](tbb::blocked_range<size_t>r){
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            for (int nuc_idx=0; nuc_idx<4; nuc_idx++) {
                /*if (idx==24047&&nuc_idx==2) {
                    raise(SIGTRAP);
                }*/
                bulid_idx_tree(idx_node_vec[idx][nuc_idx], out.indexes[idx][nuc_idx]);
            }
        }
    });
    out.tree_height=tree.get_max_level();
    return out;
}