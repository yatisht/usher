#include "ripples.hpp"
void prep_output(const Pruned_Sample &sample,
                 Ripples_Mapper_Output_Interface &out, size_t node_size) {
    const auto &init_mut = sample.sample_mutations;
    out.mut_count_out.resize(node_size * (init_mut.size() + 1));
    out.mut_out.resize(node_size);
    out.is_sibling.resize(node_size);
    std::vector<Ripples_Mapper_Mut> &init = out.mut_out[0];
    init.reserve(init_mut.size() + 1);
    for (size_t mut_idx = 0; mut_idx < init_mut.size(); mut_idx++) {
        init.emplace_back(init_mut[mut_idx], mut_idx);
        out.mut_count_out[mut_idx * node_size].set(mut_idx, true);
    }
    out.mut_count_out[(init_mut.size()) * node_size].set(init_mut.size(), true);
    init.emplace_back();
}
struct Mapper_Op_Common{
    Ripples_Mapper_Output_Interface &out;
    const size_t stride;
    const size_t mut_count;
    Mut_Count_t& get_mut_count_ref(size_t node_idx,unsigned short mut_idx){
        return out.mut_count_out[mut_idx*stride+node_idx];
    }
    Mut_Count_t& get_final_mut_count_ref(size_t node_idx){
        return out.mut_count_out[mut_count*stride+node_idx];
    }
    std::vector<Ripples_Mapper_Mut>& get_mut(size_t node_idx){
        return out.mut_out[node_idx];
    }
    void set_sibling(size_t node_idx,bool is_sibling){
        out.is_sibling[node_idx]=is_sibling;
    }
};
std::pair<const MAT::Node*,MAT::Mutation> get_parent_mut(const MAT::Node* node,int pos){
    node=node->parent;
    while (node!=nullptr)
    {
        for(auto mut:node->mutations){
            if(mut.position==pos){
                return std::make_pair(node,mut);
            }
        }
        node=node->parent;
    }
    return std::make_pair(nullptr,MAT::Mutation());
}
struct Mapper_Op : public tbb::task {
    const std::vector<Ripples_Mapper_Mut> &parent_muts;
    Mapper_Op_Common cfg;
    const MAT::Node *node;
    Mapper_Op(const std::vector<Ripples_Mapper_Mut> &parent_muts,
              Mapper_Op_Common cfg,
              const MAT::Node *node)
        : parent_muts(parent_muts), cfg(cfg), node(node) {}
    void only_from_parent(const Ripples_Mapper_Mut& mut,size_t this_idx,unsigned short& mut_accumulated,std::vector<Ripples_Mapper_Mut>& this_mut_out){
        auto is_valid = mut.valid();
                if (mut.mut_idx != mut.NULL_MUT_IDX) {
                    cfg.get_mut_count_ref(this_idx,mut.mut_idx).set(
                        mut_accumulated, is_valid);
                }
                if (is_valid) {
                    mut_accumulated++;
                }
                this_mut_out.push_back(mut);
    }
    tbb::task *execute() override {
        auto iter = parent_muts.begin();
        auto end = parent_muts.end();
        unsigned short mut_accumulated = 0;
        auto this_idx = node->dfs_idx;
        auto &this_mut_out = cfg.get_mut(this_idx);
        bool is_sibling=false;
        bool have_not_shared=false;
        if (node->identifier=="node_100005")
        {
            //raise(SIGTRAP);
        }
        bool has_common=false;
        if(this_idx!=0){
        for (const auto &this_mut : node->mutations) {
            while (iter->position < this_mut.position) {
                only_from_parent(*iter,this_idx,mut_accumulated,this_mut_out);
                iter++;
            }
            if (iter->position == this_mut.position) {
                if (iter->mut_idx != iter->NULL_MUT_IDX) {
                    // Position being counted
                    // Need descendant to know no matter whether it is valid
                    this_mut_out.emplace_back(*iter, this_mut.mut_nuc);
                    auto &mut_for_descendant = this_mut_out.back();
                    // May be splited
                    bool have_to_be_valid =
                        mut_for_descendant.valid() && iter->valid();
                        cfg.get_mut_count_ref(this_idx,iter->mut_idx).set(
                        mut_accumulated, have_to_be_valid);
                    if(!mut_for_descendant.valid()){
                        has_common=true;
                    }
                    if((!iter->valid())&&mut_for_descendant.valid()){
                        is_sibling=true;
                    }
                    if (have_to_be_valid) {
                        mut_accumulated++;
                    }
                } else {
                    assert(iter->valid());
                    if (this_mut.mut_nuc != iter->dest_mut) {
                        this_mut_out.emplace_back(*iter, this_mut.mut_nuc);
                        mut_accumulated++;
                    }
                }
                iter++;
            } else {
                //Only at this node, back mutations for descendant, but maybe split
                if (this_mut.mut_nuc != this_mut.ref_nuc) {
                    this_mut_out.emplace_back(this_mut);
                    is_sibling=true;
                }
            }
        }
        while (iter < end) {
                only_from_parent(*iter,this_idx,mut_accumulated,this_mut_out);
                iter++;
        }
        //Artifically increase mutation by one for spliting?
        if(node->children.empty()||!has_common){
            is_sibling=true;
        }
        cfg.get_final_mut_count_ref(this_idx).set(mut_accumulated,false);
        }
        cfg.set_sibling(this_idx,is_sibling);
        auto cont=new (allocate_continuation()) tbb::empty_task;
        cont->set_ref_count(node->children.size());
        for(const auto child:node->children){
            cont->spawn(* new(cont->allocate_child()) Mapper_Op(this_mut_out,cfg,child));
        }
        return node->children.empty()?cont:nullptr;
    }
};
void ripples_mapper(const Pruned_Sample &sample,
                    Ripples_Mapper_Output_Interface &out,
                    const std::vector<MAT::Node *> &dfs_ordered_nodes,
                    const MAT::Node *root) {
    auto temp=get_parent_mut(root,0);
    auto node_size = dfs_ordered_nodes.size();
    auto mut_size = sample.sample_mutations.size();
    prep_output(sample, out, node_size);
    Mapper_Op_Common cfg{out,node_size,mut_size};
    tbb::task::spawn_root_and_wait(*new( tbb::task::allocate_root())Mapper_Op(out.mut_out[0],cfg,root) );
}
