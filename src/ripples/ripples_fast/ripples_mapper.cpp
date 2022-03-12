#include "ripples.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <cassert>
#include <csignal>
#include <stack>
#include <vector>
struct Mapper_Cont:public tbb::task {
    std::vector<Ripples_Mapper_Mut> muts;
    tbb::task* execute() override {
        return nullptr;
    }
};
void prep_output(const Pruned_Sample &sample,std::vector<Ripples_Mapper_Mut> &init,
                 Ripples_Mapper_Output_Interface &out, size_t node_size) {
    const auto &init_mut = sample.sample_mutations;
    auto out_size=node_size * (init_mut.size() + 1);
    out.mut_count_out.reserve(out_size+8);
    out.mut_count_out.resize(out_size);
    out.is_sibling.resize(node_size);
    init.reserve(init_mut.size() + 1);
    for (size_t mut_idx = 0; mut_idx < init_mut.size(); mut_idx++) {
        init.emplace_back(init_mut[mut_idx], mut_idx);
        out.mut_count_out[mut_idx * node_size].set(mut_idx, true);
    }
    out.mut_count_out[(init_mut.size()) * node_size].set(init_mut.size(), true);
    init.emplace_back();
}
struct Mapper_Op_Common {
    Ripples_Mapper_Output_Interface &out;
    const size_t stride;
    const size_t mut_count;
    const std::vector<int>& idx_map;
    const std::vector<bool>& do_parallel;
    const std::vector<Mapper_Info> &traversal_track;
    const int skip_idx;
    const unsigned short tree_height;
    Mut_Count_t& get_mut_count_ref(size_t node_idx,unsigned short mut_idx) {
        return out.mut_count_out[mut_idx*stride+node_idx];
    }
    Mut_Count_t& get_final_mut_count_ref(size_t node_idx) {
        return out.mut_count_out[mut_count*stride+node_idx];
    }
    void set_sibling(size_t node_idx,bool is_sibling) {
        out.is_sibling[node_idx]=is_sibling;
    }
};
std::pair<const MAT::Node*,MAT::Mutation> get_parent_mut(const MAT::Node* node,int pos) {
    node=node->parent;
    while (node!=nullptr) {
        for(auto mut:node->mutations) {
            if(mut.position==pos) {
                return std::make_pair(node,mut);
            }
        }
        node=node->parent;
    }
    return std::make_pair(nullptr,MAT::Mutation());
}
static void only_from_parent(const Ripples_Mapper_Mut& mut,size_t this_idx,unsigned short& mut_accumulated,std::vector<Ripples_Mapper_Mut>& this_mut_out,Mapper_Op_Common& cfg) {
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
void set_mut_count(int this_idx,
                   std::vector<Ripples_Mapper_Mut> &this_mut_out,
                   const MAT::Mutation *this_mut_start,
                   const MAT::Mutation *this_mut_end, bool is_leaf,
                   const std::vector<Ripples_Mapper_Mut> &parent_muts,
                   Mapper_Op_Common &cfg) {
    auto iter = parent_muts.begin();
    auto end = parent_muts.end();
    unsigned short mut_accumulated = 0;
    bool has_unique=false;
    bool has_common=false;
    bool is_sibling=false;
    auto reserve_size=(this_mut_end-this_mut_start)+parent_muts.size();
    if (reserve_size<0) {
        raise(SIGTRAP);
    }
    this_mut_out.reserve(reserve_size);
    if(this_idx!=0) {
        for (; this_mut_start<this_mut_end; this_mut_start++) {
            const auto& this_mut=*this_mut_start;
            while (iter->position < this_mut.position) {
                only_from_parent(*iter,this_idx,mut_accumulated,this_mut_out,cfg);
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
                    if(!mut_for_descendant.valid()) {
                        has_common=true;
                    }
                    if(mut_for_descendant.valid()) {
                        has_unique=true;
                    }
                    if (have_to_be_valid) {
                        mut_accumulated++;
                    }
                } else {
                    assert(iter->valid());
                    if (this_mut.mut_nuc != iter->dest_mut) {
                        this_mut_out.emplace_back(*iter, this_mut.mut_nuc);
                        mut_accumulated++;
                    } else {
                        has_common=true;
                    }
                }
                iter++;
            } else {
                //Only at this node, back mutations for descendant, but maybe split
                if (this_mut.mut_nuc != this_mut.ref_nuc) {
                    this_mut_out.emplace_back(this_mut);
                    has_unique=true;
                } else {
                    has_common=true;
                }
            }
        }
        while (iter < end) {
            only_from_parent(*iter,this_idx,mut_accumulated,this_mut_out,cfg);
            iter++;
        }
        //Artifically increase mutation by one for spliting?
        is_sibling=((has_unique && !is_leaf && (has_common)) || \
                    (is_leaf && (has_common)) || (!has_unique && !is_leaf ));
        cfg.get_final_mut_count_ref(this_idx).set(mut_accumulated,false);
    }
    cfg.set_sibling(this_idx,is_sibling);
}
void merge_mutation_only(std::vector<Ripples_Mapper_Mut>& this_mut_out,const MAT::Node* node,const std::vector<Ripples_Mapper_Mut>& parent_muts) {
    auto iter = parent_muts.begin();
    auto end = parent_muts.end();
    this_mut_out.reserve(node->mutations.size()+parent_muts.size());
    for (const auto &this_mut : node->mutations) {
        while (iter->position < this_mut.position) {
            this_mut_out.push_back((*iter));
            iter++;
        }
        if (iter->position == this_mut.position) {
            if (iter->mut_idx != iter->NULL_MUT_IDX) {
                this_mut_out.emplace_back(*iter, this_mut.mut_nuc);
            } else {
                assert(iter->valid());
                if (this_mut.mut_nuc != iter->dest_mut) {
                    this_mut_out.emplace_back(*iter, this_mut.mut_nuc);
                }
            }
            iter++;
        } else {
            //Only at this node, back mutations for descendant, but maybe split
            if (this_mut.mut_nuc != this_mut.ref_nuc) {
                this_mut_out.emplace_back(this_mut);
            }
        }
    }
    while (iter < end) {
        this_mut_out.push_back((*iter));
        iter++;
    }
}
static void
serial_mapper(const std::vector<Ripples_Mapper_Mut> &parent_muts,
              Mapper_Op_Common &cfg,
              int start_idx) {
    const std::vector<Mapper_Info> &traversal_track=cfg.traversal_track;
    std::vector<std::vector<Ripples_Mapper_Mut>> stack;
    int end_idx=traversal_track[start_idx].sibling_start_idx;
    auto base_level = traversal_track[start_idx].level;
    stack.reserve(cfg.tree_height - base_level);
    if (cfg.tree_height<base_level) {
        raise(SIGTRAP);
    }
    stack.emplace_back();
    set_mut_count(start_idx, stack.back(), traversal_track[start_idx].begin,
                  traversal_track[start_idx].end,
                  traversal_track[start_idx].is_leaf,
                  stack.size()==1 ? parent_muts : stack[stack.size()-2], cfg);

    for (auto cur_idx=start_idx+1; cur_idx<end_idx; cur_idx++) {
        if (cur_idx==cfg.skip_idx) {
            cur_idx=traversal_track[cur_idx].sibling_start_idx;
        }
        if (cur_idx>=end_idx) {
            break;
        }
        auto next_level=traversal_track[cur_idx].level;
#ifndef NDEBUG
        auto cur_level=base_level+stack.size();
#endif
        assert(next_level<=cur_level);
        if (next_level<=base_level) {
            raise(SIGTRAP);
        }
        stack.resize(next_level-base_level);
        stack.emplace_back();
        set_mut_count(cur_idx, stack.back(), traversal_track[cur_idx].begin,
                      traversal_track[cur_idx].end,
                      traversal_track[cur_idx].is_leaf,
                      stack.size()==1 ? parent_muts : stack[stack.size()-2], cfg);
    }
}
struct Mapper_Op : public tbb::task {
    const std::vector<Ripples_Mapper_Mut> &parent_muts;
    Mapper_Op_Common& cfg;
    const MAT::Node *node;
    Mapper_Op(const std::vector<Ripples_Mapper_Mut> &parent_muts,
              Mapper_Op_Common& cfg,
              const MAT::Node *node)
        : parent_muts(parent_muts), cfg(cfg), node(node) {}
    tbb::task *execute() override {
        auto dfs_idx=node->dfs_idx;
        auto this_idx = cfg.idx_map[dfs_idx];
        if(this_idx<0||this_idx==cfg.skip_idx) {
            return nullptr;
        }
        if (!cfg.do_parallel[dfs_idx]) {
            serial_mapper(parent_muts,cfg,this_idx);
            return nullptr;
        }

        auto cont=new (allocate_continuation()) Mapper_Cont;
        auto &this_mut_out = cont->muts;

        if (node->identifier=="node_7194") {
            //raise(SIGTRAP);
        }
        //if (this_idx>=0)
        //{
        set_mut_count(this_idx, this_mut_out, node->mutations.data(),
                      node->mutations.data() + node->mutations.size(),
                      node->children.empty(), parent_muts, cfg);
        //}else{
        // merge_mutation_only(this_mut_out,node,parent_muts);
        //}

        cont->set_ref_count(node->children.size());
        if (this_idx == 0) {
            for (const auto child : node->children) {
                cont->spawn(*new (cont->allocate_child())
                            Mapper_Op(parent_muts, cfg, child));
            }
        } else {

            for (const auto child : node->children) {
                cont->spawn(*new (cont->allocate_child())
                            Mapper_Op(this_mut_out, cfg, child));
            }
        }
        return node->children.empty()?cont:nullptr;
    }
};
void ripples_mapper(const Pruned_Sample &sample,
                    Ripples_Mapper_Output_Interface &out,
                    size_t node_size,
                    const std::vector<int>& idx_map,
                    const std::vector<bool>& do_parallel,
                    const std::vector<Mapper_Info> &traversal_track,
                    const unsigned short tree_height,
                    const MAT::Node *root,
                    const MAT::Node *skip_node) {
    auto temp=get_parent_mut(root,0);
    auto mut_size = sample.sample_mutations.size();
    std::vector<Ripples_Mapper_Mut> root_muts;
    prep_output(sample, root_muts,out, node_size);
    Mapper_Op_Common cfg{out,
                         node_size,
                         mut_size,
                         idx_map,
                         do_parallel,
                         traversal_track,
                         idx_map[skip_node->dfs_idx],
                         tree_height};
    tbb::task::spawn_root_and_wait(*new (tbb::task::allocate_root())
                                   Mapper_Op(root_muts, cfg, root));
}
