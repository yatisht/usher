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
#include <tuple>
#include <unordered_map>
#include <vector>
#include "src/usher-sampled/usher.hpp"
#include "src/usher-sampled/mapper.hpp"
#include "src/usher-sampled/place_sample.hpp"
#define INDEX_END_POSITION INT_MAX
namespace MAT = Mutation_Annotated_Tree;
union fixed_tree_search_mutation{
    struct{
        int position;
        uint8_t chrom;
        uint8_t mut_nuc;
        union{
        uint16_t range;
        uint8_t par_nuc;
        };
        int index_tree_idx;
    };
    std::array<int, 3> index_tree_idx_array;
    int get_end_range() const{
        if (mut_nuc==0xf) {
            return position+range;
        }else {
            return position;
        }
    }
};
int update_idx(int& idx, int dfs_idx,int dfs_end_idx, const std::vector<index_ele>& dfs_elem){
    if (idx==INDEX_END_POSITION) {
        return INDEX_END_POSITION;
    }
    int seek_distance=0;
    for (; dfs_elem[idx].dfs_idx_end<dfs_idx;idx++) {
        seek_distance++;
    }
    if (seek_distance>8) {
        fprintf(stderr, "Seek distance %d\n",seek_distance);
    }
    if (dfs_elem[idx].dfs_idx>dfs_end_idx) {
        return INDEX_END_POSITION;
    }
    if (dfs_elem[idx].dfs_idx<dfs_idx) {
        if (dfs_elem[idx].dfs_idx_end<dfs_end_idx) {
            fprintf(stderr, "Not covering\n");
            raise(SIGTRAP);
        }
        if (dfs_elem[idx].children_idx==-1) {
            return INDEX_END_POSITION;
        }else {
            int ignored=dfs_elem[idx].children_idx;
            return update_idx(ignored, dfs_idx,dfs_end_idx,dfs_elem);
        }
    }
    return idx;
}
static bool check_more_bit_set(int in){
    return in&(in-1);
}
struct Output_No_Lock{
    std::vector<Main_Tree_Target>& targets;
    int best_par_score;
};
void register_target(Main_Tree_Target &target, int this_score,Output_No_Lock &output
    ,int dfs_idx,const std::vector<MAT::Node*>& dfs_ordered_nodes) {
    if ((dfs_idx==0)&&target.shared_mutations.empty()) {
        return;
    }
    if (output.best_par_score > this_score) {
        output.best_par_score = this_score;
        output.targets.clear();
    }
    if (output.best_par_score == this_score) {
        target.target_node=dfs_ordered_nodes[dfs_idx];
        target.parent_node=dfs_ordered_nodes[dfs_idx]->parent;
        output.targets.push_back(std::move(target));
    }
}
static int
add_existing_mut(std::vector<fixed_tree_search_mutation>::iterator &iter,
                 std::vector<fixed_tree_search_mutation> &output,
                 size_t dfs_idx, size_t dfs_end_idx, const Traversal_Info &in,
                 uint8_t par_nuc,Main_Tree_Target& sibling_out ) {
    
    if (iter->mut_nuc==0xf) {
        output.push_back(*iter);
        output.back().par_nuc=par_nuc;
        To_Place_Sample_Mutation mut(iter->position,iter->chrom,iter->mut_nuc);
        mut.range=iter->range;
        sibling_out.sample_mutations.push_back(mut);    
        iter++;
        return 0;
    }
    To_Place_Sample_Mutation mut(iter->position,iter->chrom,iter->mut_nuc,par_nuc);
    sibling_out.sample_mutations.push_back(mut);    
    if (check_more_bit_set(iter->mut_nuc)) {
        bool not_found=!(iter->mut_nuc&par_nuc);
        output.push_back(*iter);
        auto position=output.back().position;
        output.back().par_nuc=par_nuc;
        if (iter->index_tree_idx!=INDEX_END_POSITION) {
                output.back().index_tree_idx =
                    update_idx(iter->index_tree_idx, dfs_idx, dfs_end_idx,
                       in.indexes[position][0]);
                not_found &= (output.back().index_tree_idx == INDEX_END_POSITION);
        }else {
            output.back().index_tree_idx =INDEX_END_POSITION;
        }
        iter++;
        output.emplace_back();
        for (int nuc_idx=1; nuc_idx<4; nuc_idx++) {
            if (iter->index_tree_idx_array[nuc_idx-1]!=INDEX_END_POSITION) {
                output.back().index_tree_idx_array[nuc_idx-1] =
                    update_idx(iter->index_tree_idx_array[nuc_idx-1], dfs_idx, dfs_end_idx,
                       in.indexes[position][nuc_idx]);
                not_found &= (output.back().index_tree_idx_array[nuc_idx-1] == INDEX_END_POSITION);
            }else {
                output.back().index_tree_idx_array[nuc_idx-1] =INDEX_END_POSITION;
            }
        }
        iter++;
        return not_found?1:0;
    }else {
        output.push_back(*iter);
        output.back().par_nuc=par_nuc;
        output.back().index_tree_idx=update_idx(iter->index_tree_idx, dfs_idx, dfs_end_idx, in.indexes[iter->position][one_hot_to_two_bit(iter->mut_nuc)]);
        iter++;
        return ((!(output.back().mut_nuc&output.back().par_nuc))
                &&(output.back().index_tree_idx==INDEX_END_POSITION))?1:0;
    }
}
static void ins_mut(const MAT::Mutations_Collection& to_insert,std::unordered_map<int, uint8_t>& avail){
    for (const auto& mut : to_insert) {
        avail.emplace(mut.get_position(),mut.get_mut_one_hot());
    }
}
static void compare(std::unordered_map<int, uint8_t> pos,const std::vector<To_Place_Sample_Mutation>& ref){
        bool trap=false;
    for (const auto& ref_nuc : ref) {
        if (ref_nuc.mut_nuc==0xf) {
            for (int idx=ref_nuc.position; idx<=ref_nuc.get_end_range(); idx++) {
                auto iter=pos.find(idx);
                if (iter==pos.end()||iter->second!=0xf) {
                    fprintf(stderr, "%d to N Not found\n",idx);
                    trap=true;
                }else {
                    pos.erase(iter);
                }
            }
        }else {
            auto iter=pos.find(ref_nuc.position);
                if (iter==pos.end()||iter->second!=ref_nuc.mut_nuc) {
                    fprintf(stderr, "%d to %d Not found\n",ref_nuc.position,ref_nuc.mut_nuc);
                    trap=true;
                }else {
                    pos.erase(iter);
                }
        }
    }
    for (const auto& remaining : pos) {
        if (MAT::Mutation::refs[remaining.first]!=remaining.second) {
            fprintf(stderr, "Extra mutation at %d to %d\n",remaining.first,remaining.second);
            trap=true;
        }
    }
    if (trap) {
        raise(SIGTRAP);
    }
}
static void check_output(const Main_Tree_Target& to_check,const std::vector<To_Place_Sample_Mutation>& ref){
    std::unordered_map<int, uint8_t> pos;
    pos.reserve(ref.size());
    for (const auto& mut : to_check.sample_mutations) {
        if (mut.position==INT_MAX) {
            continue;
        }
        if (mut.mut_nuc==0xf) {
            auto end_range=mut.get_end_range();
            if (end_range<mut.position) {
                fprintf(stderr, "inversion\n");
                raise(SIGTRAP);
            }
            for (int idx=mut.position; idx<=end_range; idx++) {
                pos.emplace(idx,0xf);
            }
        }else {
            pos.emplace(mut.position,mut.mut_nuc);
        }
    }
    ins_mut(to_check.shared_mutations, pos);
    auto node=to_check.parent_node;
    while (node) {
        ins_mut(node->mutations, pos);
        node=node->parent;
    }
    compare(pos, ref);
}
static void check_descendant(const std::vector<fixed_tree_search_mutation>& to_check,const std::vector<To_Place_Sample_Mutation>& ref,const MAT::Node* node){
        std::unordered_map<int, uint8_t> pos;
    pos.reserve(ref.size());
    for (int idx=0;idx<to_check.size();idx++) {
        const auto& mut=to_check[idx];
        if (mut.mut_nuc==0xf) {
            for (int idx=mut.position; idx<=mut.get_end_range(); idx++) {
                pos.emplace(idx,0xf);
            }
        }else {
            pos.emplace(mut.position,mut.mut_nuc);
            if (check_more_bit_set(mut.mut_nuc)) {
                idx++;
            }
        }
    }
    while (node) {
        ins_mut(node->mutations, pos);
        node=node->parent;
    }
    compare(pos, ref);
}
static int merge_mutations(std::vector<fixed_tree_search_mutation>& parent,const MAT::Mutation* node_mut, size_t node_mut_count,
    std::vector<fixed_tree_search_mutation>& descendant_output,
    int dfs_idx, int dfs_end_idx,const Traversal_Info& in,Output_No_Lock& to_out
    ,const std::vector<MAT::Node*>& dfs_ordered_nodes,int last_lower_bound
){
    Main_Tree_Target sibling_out;
    int lower_bound=0;
    int parsimony_score=0; 
    auto iter=parent.begin();
    auto end=parent.end();
    for (int idx=0; idx<node_mut_count; idx++) {
        const auto& this_mut=node_mut[idx];
        while (iter->get_end_range()<this_mut.get_position()) {
            if (!(iter->mut_nuc==0xf||(iter->mut_nuc&iter->par_nuc))) {
                parsimony_score++;
            }
            lower_bound+=add_existing_mut(iter, descendant_output, dfs_idx, dfs_end_idx, in, iter->par_nuc,sibling_out);
        }
        if (iter->position<=this_mut.get_position()&&iter->mut_nuc==0xf) {
            continue;
        }
        if (iter->position==this_mut.get_position()) {
            if (iter->mut_nuc==0xf) {
                raise(SIGTRAP);
            }
            if (iter->mut_nuc==this_mut.get_mut_one_hot()) {
                sibling_out.shared_mutations.push_back(this_mut);
                if (check_more_bit_set(iter->mut_nuc)) {
                    fprintf(stderr, "Mult bit set\n");
                    raise(SIGTRAP);
                }
                if (iter->index_tree_idx==INDEX_END_POSITION) {
                    fprintf(stderr, "pos %d bound,idx %zu \n",iter->position,node_mut+idx);
                    raise(SIGTRAP);
                }
                iter++;
                continue;
            }else {
                if (iter->mut_nuc&this_mut.get_mut_one_hot()) {
                    sibling_out.shared_mutations.push_back(this_mut);
                    auto coinciding_two_bit=one_hot_to_two_bit(iter->mut_nuc&this_mut.get_mut_one_hot());
                    auto trap=false;
                    if (coinciding_two_bit==0) {
                        trap=iter->index_tree_idx==INDEX_END_POSITION;
                    }else {
                        trap=(iter+1)->index_tree_idx_array[coinciding_two_bit-1]==INDEX_END_POSITION;
                    }
                    if (trap) {
                        fprintf(stderr, "pos %d bound, idx %zu, co2bit %d \n",iter->position,node_mut+idx,coinciding_two_bit);
                        raise(SIGTRAP);
                    }
                }else {
                    parsimony_score++;
                    sibling_out.splited_mutations.push_back(this_mut);
                }
                lower_bound+=add_existing_mut(iter, descendant_output, dfs_idx, dfs_end_idx, in, this_mut.get_mut_one_hot(),sibling_out);
            }
        }else {
            sibling_out.splited_mutations.push_back(this_mut);
            fixed_tree_search_mutation mut;
            mut.position=this_mut.get_position();
            mut.chrom=this_mut.get_chromIdx();
            mut.par_nuc=this_mut.get_mut_one_hot();
            mut.mut_nuc=this_mut.get_par_one_hot();
            if (this_mut.get_descendant_mut()&this_mut.get_par_one_hot()) {
                int temp=0;
                mut.index_tree_idx=update_idx(temp, dfs_idx, dfs_end_idx, 
                    in.indexes[this_mut.get_position()][one_hot_to_two_bit(this_mut.get_par_one_hot())]);
            }else {
                mut.index_tree_idx=INDEX_END_POSITION;
                lower_bound++;
            }
            descendant_output.push_back(mut);
        }
    }
    while (iter<end) {
        if (!(iter->mut_nuc==0xf||(iter->mut_nuc&iter->par_nuc))) {
            parsimony_score++;
        }
        lower_bound+=add_existing_mut(iter, descendant_output, dfs_idx, dfs_end_idx, in, iter->par_nuc,sibling_out);
    }
    if (parsimony_score<last_lower_bound) {
        fprintf(stderr, "Bound par %d bound %d\n",parsimony_score,last_lower_bound);
        raise(SIGTRAP);
    }
    //check_descendant(descendant_output, ref, dfs_ordered_nodes[dfs_idx]);
    register_target(sibling_out, parsimony_score, to_out,dfs_idx,dfs_ordered_nodes);
    return lower_bound;
}
move_type* place_sample_fixed_idx(const Traversal_Info &in,
                  Sample_Muts* to_search,
                  const std::vector<MAT::Node*>& dfs_ordered_nodes) {
    struct stack_content{
        std::vector<fixed_tree_search_mutation> last_mut;
        int dfs_end_idx;
        int lower_bound;
    };
    const auto &mutations=to_search->muts;
    auto res=new move_type;
    std::get<1>(*res)=to_search;
    Output_No_Lock output{std::get<0>(*res),INT_MAX};
    std::vector<stack_content> stack;
    stack.reserve(in.tree_height);
    stack.emplace_back();
    auto& base_stack=stack.back().last_mut;
    auto dfs_end_idx=in.traversal_track[0].dfs_end_idx;
    stack.back().dfs_end_idx=dfs_end_idx;
    //debug
    stack.back().lower_bound=0;
    
    for (const auto& mut : mutations) {
        fixed_tree_search_mutation mut_out;
        mut_out.chrom=mut.chrom_idx;
        mut_out.position=mut.position;
        mut_out.mut_nuc=mut.mut_nuc;    
        if (mut.mut_nuc==0xf) {
            mut_out.range=mut.range;
            base_stack.push_back(mut_out);
        }else if (check_more_bit_set(mut.mut_nuc)) {
            mut_out.par_nuc=mut.par_nuc;
            auto position=mut.position;
            bool not_found=!(mut.par_nuc&mut.mut_nuc);
            if (1 & mut.mut_nuc) {
                int idx = 0;
                mut_out.index_tree_idx =
                update_idx(idx, 0, dfs_end_idx, in.indexes[position][0]);
                not_found&=(mut_out.index_tree_idx==INDEX_END_POSITION);
            }else {
                mut_out.index_tree_idx=INDEX_END_POSITION;
            }
            base_stack.push_back(mut_out);
            fixed_tree_search_mutation next;
            for (int nuc_idx = 1; nuc_idx < 4; nuc_idx++) {
                if ((1 << nuc_idx) & mut.mut_nuc) {
                    int idx = 0;
                    next.index_tree_idx_array[nuc_idx - 1] = update_idx(
                        idx, 0, dfs_end_idx, in.indexes[position][nuc_idx]);
                    not_found&=(next.index_tree_idx_array[nuc_idx - 1]==INDEX_END_POSITION);
                } else {
                    next.index_tree_idx_array[nuc_idx - 1] = INDEX_END_POSITION;
                }
            }
            if (not_found) {
                stack.back().lower_bound++;
            }
            base_stack.push_back(next);
        }else {
            mut_out.par_nuc=mut.par_nuc;
            auto nuc_two_bit=one_hot_to_two_bit(mut.mut_nuc);
             int idx=0;
            mut_out.index_tree_idx=update_idx(idx, 0, dfs_end_idx,in.indexes[mut.position][nuc_two_bit]);
            if (mut_out.index_tree_idx==INDEX_END_POSITION) {
                stack.back().lower_bound++;
            }
            base_stack.push_back(mut_out);
        }
    }
    fixed_tree_search_mutation mut_out;
    mut_out.position=INT_MAX;
    mut_out.range=0;
    mut_out.mut_nuc=0xf;
    base_stack.push_back(mut_out);

    for (int cur_idx=1; cur_idx<dfs_end_idx; cur_idx++) {
        const auto & curr_node=in.traversal_track[cur_idx];
        while (stack.back().dfs_end_idx<curr_node.dfs_end_idx) {
            stack.pop_back();
        }
        /*if (stack.back().node!=dfs_ordered_nodes[cur_idx]->parent) {
            fprintf(stderr, "parent mismatch\n");
            raise(SIGTRAP);
        }*/
        stack.emplace_back();
        stack.back().dfs_end_idx=curr_node.dfs_end_idx;


        stack.back().lower_bound=merge_mutations(stack[stack.size()-2].last_mut, curr_node.mutation_start,
            curr_node.mutation_size, stack.back().last_mut, cur_idx, 
            curr_node.dfs_end_idx, in, output,dfs_ordered_nodes,stack[stack.size()-2].lower_bound);
    }
    for(const auto& result:std::get<0>(*res)){
        check_output(result, to_search->muts);
    }
    return res;
}
