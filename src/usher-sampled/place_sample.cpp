#include "mapper.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for.h>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#define DELETED_NODE_THRESHOLD 1000
static void update_possible_descendant_alleles(
    const MAT::Mutations_Collection &mutations_to_set,
    MAT::Node *node) {
    std::unordered_map<int, uint8_t> alleles;
    alleles.reserve(mutations_to_set.size());
    for (auto &mut : mutations_to_set) {
        alleles.emplace(mut.get_position(), mut.get_mut_one_hot());
    }
    while (!alleles.empty() && node) {
        for (auto &mut : node->mutations) {
            auto iter = alleles.find(mut.get_position());
            if (iter != alleles.end()) {
                if ((mut.get_descendant_mut() & iter->second) ==
                    iter->second) {
                    alleles.erase(iter);
                } else {
                    mut.set_descendant_mut(mut.get_descendant_mut()|iter->second);
                }
            }
        }
        node->bfs_index++;
        node = node->parent;
    }
    while (node) {
        node->bfs_index++;
        node = node->parent;
    }
}
static void gather_par_mutation_step(std::unordered_map<int, int>& to_find,MAT::Mutations_Collection& upstream, MAT::Mutations_Collection& output){
    for (const auto& mut : upstream) {
        auto iter=to_find.find(mut.get_position());
        if (iter!=to_find.end()) {
            output[iter->second].set_par_one_hot(mut.get_mut_one_hot());
            to_find.erase(iter);
        }
    }
}
static void gather_par_mut(std::unordered_map<int, int>& to_find,MAT::Node*& node, MAT::Mutations_Collection& output){
    while (node&&(!to_find.empty())) {
        gather_par_mutation_step(to_find,node->mutations,output);
        node=node->parent;
    }
    for(auto& temp: to_find){
        output[temp.second].set_par_one_hot(output[temp.second].get_ref_one_hot());
    }
}

static void discretize_mutations(std::vector<To_Place_Sample_Mutation> &in,
                                 MAT::Mutations_Collection &shared_mutations,
                                 MAT::Node *parent_node,
                                 MAT::Mutations_Collection &out) {
    out.reserve(in.size());
    std::unordered_map<int, int> par_nuc_idx;
    assert(in.back().position==INT_MAX);
    for (size_t idx=0;idx<(in.size()-1) ; idx++) {
        const auto & mut=in[idx];
        if (mut.mut_nuc == 0xf) {
            for (int pos = mut.position; pos <= mut.get_end_range(); pos++) {
                par_nuc_idx.emplace(pos, out.size());
                out.push_back(MAT::Mutation(mut.chrom_idx, pos, 0, 0xf));
                out.back().set_descendant_mut(0xf);
            }
        } else {
            out.push_back(MAT::Mutation(mut.chrom_idx, mut.position,
                                        mut.par_nuc, mut.mut_nuc));
            out.back().set_descendant_mut(mut.mut_nuc);
        }
    }
    gather_par_mutation_step(par_nuc_idx, shared_mutations,out);
    gather_par_mut(par_nuc_idx, parent_node, out);
}
static MAT::Node* add_children(MAT::Node* target_node,MAT::Node* sample_node){
    MAT::Node* deleted_node=nullptr;
    if ((target_node->children.size()+1)>=target_node->children.capacity()) {
        MAT::Node* new_target_node=new MAT::Node(*target_node);
        for(auto child:new_target_node->children){
            child->parent=new_target_node;
        }
        new_target_node->children.reserve(4*new_target_node->children.size());
        auto& parent_children=target_node->parent->children;
        auto iter=std::find(parent_children.begin(),parent_children.end(),target_node);
        *iter=new_target_node;
        deleted_node=target_node;
        target_node=new_target_node;
    }
    target_node->children.push_back(sample_node);
    sample_node->parent=target_node;
    return deleted_node;
}
static  MAT::Node *update_main_tree(Main_Tree_Target &target,
                                         std::string &&sample_string) {
    // Split branch?
    MAT::Node* deleted_node=nullptr;
    MAT::Node *sample_node = new MAT::Node;
    sample_node->level=target.target_node->level;
    sample_node->identifier = std::move(sample_string);
    if (sample_node->identifier=="s484400s") {
        //raise(SIGTRAP);
    }
    discretize_mutations(target.sample_mutations, target.shared_mutations, target.parent_node, sample_node->mutations);
    int sample_node_mut_count = 0;
    for (const auto &mut : sample_node->mutations) {
        if (!(mut.get_par_one_hot() & mut.get_mut_one_hot())) {
            sample_node_mut_count++;
        }
        assert(mut.get_position());
    }
    sample_node->branch_length = sample_node_mut_count;
    if (target.splited_mutations.empty() && (!target.target_node->is_leaf())) {
        deleted_node=add_children(target.target_node, sample_node);
    } else if (target.shared_mutations.empty() &&
               (!target.target_node->is_leaf())) {
        deleted_node=add_children(target.parent_node, sample_node);
    } else {
        MAT::Node* new_target_node=new MAT::Node;
        new_target_node->identifier=target.target_node->identifier;
        new_target_node->level=target.target_node->level;
        new_target_node->children.reserve(4*target.target_node->children.size());
        new_target_node->children=target.target_node->children;
        new_target_node->mutations = std::move(target.splited_mutations);
        for (auto child : new_target_node->children) {
            child->parent=new_target_node;
        }
        int target_node_mut_count = 0;
        for (const auto &mut : target.target_node->mutations) {
            if (!(mut.get_mut_one_hot() & mut.get_par_one_hot())) {
                target_node_mut_count++;
            }
        }
        new_target_node->branch_length = target_node_mut_count;
        MAT::Node *split_node = new MAT::Node;
        new_target_node->parent = split_node;
        sample_node->parent = split_node;
        split_node->identifier = "";
        split_node->level=target.target_node->level;
        split_node->parent = target.parent_node;
        split_node->mutations = std::move(target.shared_mutations);
        split_node->children.reserve(4);
        split_node->children.push_back(new_target_node);
        split_node->children.push_back(sample_node);
        split_node->branch_length = split_node->mutations.size();
        auto iter =
            std::find(target.parent_node->children.begin(),
                      target.parent_node->children.end(), target.target_node);
        if (iter==target.parent_node->children.end()||*iter!=target.target_node) {
            std::raise(SIGTRAP);
        }
        *iter = split_node;
        deleted_node=target.target_node;
        deleted_node->parent=nullptr;
        target.target_node=new_target_node;
    }
    update_possible_descendant_alleles(sample_node->mutations, sample_node->parent);
    #ifndef NDEBUG
    check_descendant_nuc(sample_node);
    check_descendant_nuc(target.target_node);
    check_descendant_nuc(sample_node->parent);
    #endif
    return deleted_node;
}

struct Finder{
    MAT::Tree& tree;
    std::vector<Sample_Muts>& to_place;
    std::tuple<std::vector<Main_Tree_Target>, int,size_t>* operator()(size_t idx){
        auto output=new std::tuple<std::vector<Main_Tree_Target>, int,size_t>;
        std::get<2>(*output)= idx;
        const std::vector<MAT::Mutation>& sample_mutations =to_place[idx].muts;
        std::vector<To_Place_Sample_Mutation> condensed_muts;
        convert_mut_type(sample_mutations,condensed_muts);
        //auto start_time=std::chrono::steady_clock::now();
        auto main_tree_out=place_main_tree(condensed_muts, tree);
        //auto duration=std::chrono::steady_clock::now()-start_time;
        //fprintf(stderr, "Search took %ld \n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
        
        std::get<0>(*output)=std::move(std::get<0>(main_tree_out));
        std::get<1>(*output)=std::get<1>(main_tree_out);
        return output;
    }    
};

struct Placer{
    std::vector<Sample_Muts>& to_place;
    std::vector<MAT::Node*>& deleted_nodes;
    size_t& count;
    std::chrono::steady_clock::time_point start_time;
#ifndef NDEBUG
    Original_State_t &ori_state;
    MAT::Tree& tree;
#endif
    size_t operator()(std::tuple<std::vector<Main_Tree_Target>, int,size_t>* in){
        //auto start_time=std::chrono::steady_clock::now();
        auto & search_result=std::get<0>(*in);
        auto idx=std::get<2>(*in);
        for (const auto& placement : search_result) {
            if (placement.parent_node!=placement.target_node->parent) {
                fprintf(stderr, "Redoing %s \n",to_place[idx].sample_name.c_str());
                return idx;
            }
        }
        std::string sample_string = to_place[idx].sample_name;
#ifndef NDEBUG
    Mutation_Set new_set;
    std::vector<To_Place_Sample_Mutation> condensed_muts_copy;
    const std::vector<MAT::Mutation>& sample_mutations =to_place[idx].muts;
    convert_mut_type(sample_mutations,condensed_muts_copy);
    new_set.reserve(sample_mutations.size());
    for (const auto &mut : sample_mutations) {
        new_set.insert(mut);
    }
    ori_state.emplace(sample_string, new_set);
#endif
auto& selected_target=search_result[0];
    int min_level=0;
    for (auto& target : search_result) {
        if (target.target_node->level<min_level) {
            min_level=target.target_node->level;
            selected_target=target;
        }
    }
    fprintf(stderr, "Sample: %s\t%d\t%zu\n",sample_string.c_str(),std::get<1>(*in),search_result.size());
    auto deleted_node=update_main_tree(selected_target, std::move(sample_string));
    if (deleted_node) {
        deleted_nodes.push_back(deleted_node);
    }
    //auto duration=std::chrono::steady_clock::now()-start_time;
    //fprintf(stderr, "Placement took %ld \n",std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());
    count++;
    if (count%100==0) {
        auto duration=std::chrono::steady_clock::now()-start_time;
        auto per_sample=duration/count;
        fprintf(stderr, "placed %zu samples, took %ld msec per sample \n",count,std::chrono::duration_cast<std::chrono::milliseconds>(per_sample).count());

    }
    #ifndef NDEBUG
    if (idx%100==0) {
        check_samples(tree.root, ori_state, &tree);
    }
    #endif
    return SIZE_MAX;
    }
};
typedef tbb::flow::multifunction_node<size_t,tbb::flow::tuple<size_t>> Pusher_Node_T;
struct Pusher{
    size_t& idx;
    size_t max_idx;
    std::vector<MAT::Node*> deleted_nodes;
    void operator()(size_t idx_in,Pusher_Node_T::output_ports_type& out ){
        if (idx_in==SIZE_MAX&&(deleted_nodes.size()<DELETED_NODE_THRESHOLD)) {
            if (idx<max_idx) {
                std::get<0>(out).try_put(idx);
                idx++;
            }
            return;
        }
        std::get<0>(out).try_put(idx_in);
    }
};
void place_sample(std::vector<Sample_Muts> &sample_to_place, MAT::Tree &main_tree,int batch_size
#ifndef NDEBUG
                  ,
                  Original_State_t &ori_state
#endif
) {
    size_t idx=0;
    size_t idx_max=sample_to_place.size();
    std::vector<MAT::Node*> deleted_nodes;
    deleted_nodes.reserve(DELETED_NODE_THRESHOLD);
    auto start_time=std::chrono::steady_clock::now();
    size_t count=0;
    while (idx < idx_max) {
        tbb::flow::graph g;
        Pusher_Node_T init(g, 1, Pusher{idx, idx_max, deleted_nodes});
        tbb::flow::function_node<
            size_t, std::tuple<std::vector<Main_Tree_Target>, int, size_t> *>
            searcher(g, tbb::flow::unlimited,
                     Finder{main_tree,sample_to_place});
        tbb::flow::make_edge(std::get<0>(init.output_ports()),searcher);
        tbb::flow::function_node<std::tuple<std::vector<Main_Tree_Target>, int, size_t> *,size_t>
            placer(g, 1,
                     Placer{sample_to_place,deleted_nodes,count,start_time
                     #ifndef NDEBUG
                     ,ori_state,main_tree
                     #endif
                     });
        tbb::flow::make_edge(searcher,placer);
        tbb::flow::make_edge(placer,init);
        for(int temp=0;temp<batch_size;temp++){
            init.try_put(SIZE_MAX);
        }
        g.wait_for_all();
        for(auto node:deleted_nodes){
            delete node;
        }
        deleted_nodes.clear();
    }
#ifndef NDEBUG
    /*std::vector<Sampled_Tree_Node *> output;
    sample_tree_dfs(sampled_tree_root, output);
    check_sampled_tree(main_tree, output, sampling_radius);
    fprintf(stderr, "%zu samples \n", ori_state.size());*/
    //
#endif
}