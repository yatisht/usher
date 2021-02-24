#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <utility>
#include <vector>
struct FS_result_container{
    size_t start_node_idx;
    const Mutation_Annotated_Tree::Mutation* mutation;
    std::vector<std::shared_ptr<Fitch_Sankoff_Result>*> targets;
};
/*
struct Node_Collate_Sorter{
    bool operator()(const FS_result_container& first, const FS_result_container& second){
        return first.start_node_idx<second.start_node_idx;
    }
};
*/
struct Mutation_Collate_Sorter{
    bool operator()(const FS_result_container* first, const FS_result_container* second){
        return *(first->mutation)<*(second->mutation);
    }
};
struct FS_targe_t{
    size_t LCA_id;
    std::shared_ptr<Fitch_Sankoff_Result>* move;
    bool operator<(const FS_targe_t& other)const{
        return LCA_id<other.LCA_id;
    }
};
typedef std::unordered_map<MAT::Mutation, std::vector<FS_targe_t>,Mutation_Pos_Only_Hash, Mutation_Pos_Only_Comparator> Pending_FS_t;

/*
static void distribute(MAT::Mutations_Collection& src, std::vector<std::shared_ptr<Fitch_Sankoff::States_Type>>& dst,MAT::Node* this_node){
    assert(src.size()==dst.size());
    for(size_t idx=0;idx<src.size();idx++){
        #ifndef NDEBUG
        dst[idx]->emplace_back(src[idx].mut_nuc,this_node);
        #else
        dst[idx]->push_back(src[idx].mut_nuc);
        #endif
        assert(src[idx].mut_nuc==get_genotype(this_node,src[idx]));
    }
}
static std::vector<char>  get_original_states(std::vector<MAT::Node*> dfs_ordered_nodes, const std::pair<size_t, size_t> &range,MAT::Mutations_Collection target,std::vector<std::shared_ptr<Fitch_Sankoff::States_Type>>& result){
    std::vector<char> start_node_parent_states(target.size());
    MAT::Node* start_node=dfs_ordered_nodes[range.first];
    result.reserve(target.size());
    for(size_t i=0;i<target.size();i++){
        Mutation_Annotated_Tree::Mutation &m = target[i];
        char parent_state=(start_node->parent)?get_genotype(start_node->parent, m):m.ref_nuc;
        m.mut_nuc=parent_state;
        start_node_parent_states[i]=parent_state;
        auto mut_iter=start_node->mutations.find(m.position);
        if(mut_iter!=start_node->mutations.end()){
            m.mut_nuc=mut_iter->mut_nuc;
        }
        result.emplace_back(new Fitch_Sankoff::States_Type);
        result.back()->reserve(range.second-range.first);    
    }
    distribute(target, result,start_node);

    for (size_t idx=range.first+1; idx<range.second; idx++) {
        size_t loci_idx=0;
        auto par_idx=dfs_ordered_nodes[idx]->parent->index-range.first;
        for(Mutation_Annotated_Tree::Mutation& m : target){
            Fitch_Sankoff::States_Type& temp=*(result[loci_idx]);
            m.mut_nuc=temp[par_idx];
            loci_idx++;
        }
        dfs_ordered_nodes[idx]->mutations.batch_find(target);
        distribute(target, result,dfs_ordered_nodes[range.first+result[0]->size()]);
    }
    return start_node_parent_states;
}
*/
MAT::Node* get_mutation_path(MAT::Mutations_Collection& mutations,MAT::Node* src, MAT::Node* dst,std::vector<MAT::Node*>& path){
    std::unordered_set<MAT::Node*> src_to_root;
    std::vector<MAT::Node*> src_to_root_path;
    assert(src!=dst);
    MAT::Node* LCA=src;
    while (LCA) {
        src_to_root.insert(LCA);
        src_to_root_path.push_back(LCA);
        LCA=LCA->parent;
    }

    LCA=dst;
    std::vector<MAT::Node*> dst_to_root_path;
    while (!src_to_root.count(LCA)) {
        dst_to_root_path.push_back(LCA);
        LCA=LCA->parent;
    }
    assert(dst_to_root_path.empty()||dst_to_root_path.back()->parent==LCA);

    mutations.clear();
    for(MAT::Node* n:src_to_root_path){
        if (n==LCA) {
            break;
        }
        path.push_back(n);
        mutations.merge(n->mutations, MAT::Mutations_Collection::KEEP_SELF);
    }
    path.push_back(LCA);
    for (auto iter=dst_to_root_path.rbegin(); iter<dst_to_root_path.rend(); iter++) {
        path.push_back(*iter);
        mutations.merge((*iter)->mutations, MAT::Mutations_Collection::KEEP_SELF);
    }

    return LCA;
}

static void plan_FS(Pending_FS_t& mut_LCA,std::vector<FS_result_container>& out,const std::vector<MAT::Node*>& dfs_ordered_nodes){
    for(auto& mut_targets:mut_LCA){
        auto& targets=mut_targets.second;
        auto& mut=mut_targets.first;
        std::sort(targets.begin(),targets.end());
        out.push_back(FS_result_container{targets[0].LCA_id,&mut,std::vector<std::shared_ptr<Fitch_Sankoff_Result>*>({targets[0].move})});
        for (auto iter=targets.begin()+1; iter<targets.end(); iter++) {
            size_t this_node_idx=iter->LCA_id;
            size_t last_node_idx=out.back().start_node_idx;
            if(this_node_idx==last_node_idx&&check_grand_parent(dfs_ordered_nodes[this_node_idx], dfs_ordered_nodes[last_node_idx])){
                //Will have result ready from FS of last_node
                out.back().targets.push_back(iter->move);
            }else {
                out.push_back(FS_result_container{this_node_idx,&mut,std::vector<std::shared_ptr<Fitch_Sankoff_Result>*>({iter->move})});
            }
        }
    }
}
/*
static void find_parental_states(std::vector<FS_result_container>& fs_results,std::vector<MAT::Node*> dfs_ordered_nodes,std::vector<Fitch_Sankoff_Result>& repo){
    repo=std::vector<Fitch_Sankoff_Result>(fs_results.size());
    std::sort(fs_results.begin(),fs_results.end(),Node_Collate_Sorter());
    std::vector<std::vector<FS_result_container*>> collated;
    size_t last_node_idx=-1;
    for(FS_result_container& this_node:fs_results){
        if(last_node_idx!=this_node.start_node_idx){
            collated.push_back(std::vector<FS_result_container*>({&this_node}));
        }else{
            collated.back().push_back(&this_node);
        }
    }
    tbb::parallel_for(tbb::blocked_range<size_t>(0,collated.size()),[&collated,&dfs_ordered_nodes,&repo](const tbb::blocked_range<size_t>& r){
        for (size_t batch_idx=r.begin(); batch_idx<r.end(); batch_idx++) {
            auto& this_batch=collated[batch_idx];
            std::sort(this_batch.begin(), this_batch.end(),Mutation_Collate_Sorter());
            MAT::Mutations_Collection mut_collection;
            mut_collection.reserve(this_batch.size());
            for(FS_result_container* this_mut:this_batch){
                mut_collection.push_back(*(this_mut->mutation));
            }
            std::pair<size_t, size_t> range=Fitch_Sankoff::dfs_range(dfs_ordered_nodes[this_batch[0]->start_node_idx], dfs_ordered_nodes);
            std::vector<std::shared_ptr<Fitch_Sankoff::States_Type>> original_states;
            std::vector<char> LCA_parent_states=get_original_states(dfs_ordered_nodes, range, mut_collection,original_states);
            for (size_t i=0; i<this_batch.size(); i++) {
                Fitch_Sankoff_Result* this_result=&repo[this_batch[i]->idx];
                this_result->LCA_parent_state=LCA_parent_states[i];
                this_result->original_state.swap(original_states[i]);
                this_result->range=range;
                for (Fitch_Sankoff_Result** out_cont:this_batch[i]->targets ) {
                    *(out_cont)=this_result;
                }
            }
        }
    });
}
*/
static void execute_FS(std::vector<FS_result_container>& fs_results,std::vector<MAT::Node*> dfs_ordered_nodes,const Original_State_t& ori){
    tbb::parallel_for(tbb::blocked_range<size_t>(0,fs_results.size()),[&](const tbb::blocked_range<size_t>& range){
        for (size_t container_idx=range.begin(); container_idx<range.end(); container_idx++) {
            Fitch_Sankoff_Result* this_result=new Fitch_Sankoff_Result;
            std::shared_ptr<Fitch_Sankoff_Result> this_ptr(this_result);
            MAT::Node* this_node=dfs_ordered_nodes[fs_results[container_idx].start_node_idx];
            this_result->mutation=*(fs_results[container_idx].mutation);
            this_result->LCA_parent_state=get_genotype(this_node->parent, this_result->mutation);
            this_result->range=Fitch_Sankoff::dfs_range(this_node, dfs_ordered_nodes);
            if(fs_results[container_idx].targets.size()>1){
            Fitch_Sankoff::sankoff_backward_pass(this_result->range, dfs_ordered_nodes,this_result->scores,ori,this_result->mutation,this_result->LCA_parent_state);
            }else{
                assert(this_result->scores.empty());
            }
            for(std::shared_ptr<Fitch_Sankoff_Result>* store:fs_results[container_idx].targets){
                *store=this_ptr;
            }
        }
    });
}

Candidate_Moves* Parsimony_Score_Calculator::operator()(Possible_Moves* in)const{
    Pending_FS_t mut_LCA;
    Candidate_Moves* result=new Candidate_Moves;
    result->src=in->src;
    std::vector<Move_info>& to_push=result->moves;
    to_push.reserve(in->dsts.size());
    for(MAT::Node* dst:in->dsts){
        MAT::Mutations_Collection mutations;
        std::vector<MAT::Node* >path;
        MAT::Node* LCA=get_mutation_path(mutations, in->src, dst,path);
        assert(LCA!=in->src);
        assert(LCA);
        if(mutations.empty()){
            continue;
        }
        to_push.push_back(Move_info{dst,LCA,std::vector<MAT::Node*>(),std::vector<std::shared_ptr<Fitch_Sankoff_Result>>(mutations.size())});
        Move_info& move=to_push.back();
        move.path.swap(path);

        for(size_t i=0;i<mutations.size();i++){
            const Mutation_Annotated_Tree::Mutation& mutation=mutations[i];
            size_t LCA_idx=move.LCA->index;
            auto ins_result=mut_LCA.insert(std::make_pair(mutation, std::vector<FS_targe_t>()));
            /*
            auto iter=mut_LCA.find(mutation);
            if(iter==mut_LCA.end()){
                auto ins_result=mut_LCA.insert(std::make_pair(mutation, std::vector<FS_targe_t>()));
                iter=ins_result.first;
            }*/
            ins_result.first->second.push_back(FS_targe_t{LCA_idx,&move.FS_results[i]});
        }
    }
    std::vector<FS_result_container> fs_results;
    plan_FS(mut_LCA, fs_results, dfs_ordered_nodes);
    /*
    for(const FS_result_container& t:fs_results){
        if(t.targets.size()>1){
            printf("%zu \n",t.targets.size());
        }
    }
    */
    execute_FS(fs_results, dfs_ordered_nodes,original_states);
    delete in;
    return result;
}