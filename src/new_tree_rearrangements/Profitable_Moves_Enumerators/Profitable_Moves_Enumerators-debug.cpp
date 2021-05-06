#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
#include "Profitable_Moves_Enumerators.hpp"
#include "../check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "../Fitch_Sankoff.hpp"
#include <cstddef>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

void update_par_cur_nuc(MAT::Mutations_Collection::const_iterator parent_mutation_iter,dbg_iter& debug_iter,dbg_iter& debug_end){
    rewind_mutations(parent_mutation_iter->get_position(), debug_iter,debug_end);
    if (debug_iter!=debug_end&&debug_iter->position==parent_mutation_iter->get_position()) {
                debug_iter->par_nuc.push_back(parent_mutation_iter->get_par_one_hot());
                debug_iter->major_allele.push_back(parent_mutation_iter->get_all_major_allele());
                debug_iter->mutation_score_change.push_back(0);
                debug_iter->count_change.emplace_back(*parent_mutation_iter);
                debug_iter->count_change.back().set_change(0, 0,0,true);
                debug_iter++;
    }
}
void update_dbg_vector_score_only(
    int position, dbg_iter &debug_iter,dbg_iter &debug_end,
    int par_score_change) {
    rewind_mutations(position, debug_iter,debug_end);
    assert(debug_iter->position == position);
    nuc_one_hot parent_nuc = debug_iter->par_nuc.back();
    debug_iter->par_nuc.push_back(parent_nuc);
    debug_iter->major_allele.push_back(parent_nuc);
    debug_iter->mutation_score_change.push_back(par_score_change);
    debug_iter->count_change.push_back(debug_iter->count_change.back());
    debug_iter->count_change.back().set_change(0, 0,0,true);
    debug_iter++;
}
struct Mutation_Count{
    int position;
    int count;
    std::unordered_map<std::string, std::vector<MAT::Node*>> mut_count;
    Mutation_Count(){}
    Mutation_Count(int position):position(position),count(0){}
    //Mutation_Count(const MAT::Mutation& mut,int count):position(mut.get_position()),count(count){}
};

static void count_mutations(MAT::Node* start, std::vector<Mutation_Count>& mutations_count){
    size_t mut_idx=0;
    /*if (start->identifier=="Wuhan_HBCDC-HB-04_2020") {
    fputc('ab', stderr);
    }*/
    for(const MAT::Mutation& m:start->mutations){
        while (mut_idx!=mutations_count.size()&&mutations_count[mut_idx].position<m.get_position()) {
            mut_idx++;
        }if (mut_idx!=mutations_count.size()&&mutations_count[mut_idx].position==m.get_position()&&m.is_valid()) {
            mutations_count[mut_idx].count++;
            auto res=mutations_count[mut_idx].mut_count.emplace(start->parent->identifier,0);
            res.first->second.push_back(start);
        }
    }
    for(auto child:start->children){
        count_mutations(child, mutations_count);
    }
}


#ifdef DEBUG_CHECK_MAJOR_ALLELE_EACH_STEP
bool get_children_alleles_count(const MAT::Node *this_node,
                                       const MAT::Node *override_node1,
                                       const MAT::Node *override_node2,
                                       std::array<int, 4> &count,
                                       int position,nuc_one_hot this_allele,int& mutation_count) {
    int par_allele_count = 0;
    mutation_count=0;
    int overiden_par_allele_count=0;
    bool matched=false;
    for (MAT::Node *child : this_node->children) {
        auto mut_iter = child->mutations.find(position);
        auto end_iter=child->mutations.end();
                    if(mut_iter!=end_iter){
                if(mut_iter->is_valid()){
                    mutation_count++;
                }
        }
        if ((child == override_node1) || (child == override_node2)) {
             if (mut_iter == end_iter||(!mut_iter->is_valid())){
                 overiden_par_allele_count++;
             }
             matched=true;
            continue;
        }
        if (mut_iter == end_iter) {
            par_allele_count++;
        } else {
            nuc_one_hot majority_nuc =
                mut_iter->get_all_major_allele();
            for (int i = 0; i < 4; i++) {
                if (majority_nuc & (1 << i)) {
                    count[i]++;
                }
            }
        }
    }


    count[one_hot_to_two_bit(this_allele)] += par_allele_count;
    assert(mutation_count==this_node->children.size()-overiden_par_allele_count-count[one_hot_to_two_bit(this_allele)]);
    return matched;
}

static std::pair<nuc_one_hot,int> get_majority_allele(const MAT::Node *this_node,
                                const MAT::Node *override_node1,
                                const MAT::Node *override_node2,
                                nuc_one_hot override_nuc, int position,
                                nuc_one_hot this_allele,nuc_one_hot override_nuc2=0) {
    std::array<int, 4> count{0, 0, 0, 0};
    int ori_mutation_count;
    bool override_match=
        get_children_alleles_count(this_node, override_node1, override_node2,
                                 count, position,this_allele,ori_mutation_count);
    for (int i = 0; i < 4; i++) {
        if (override_nuc & (1 << i)) {
            count[i]++;
        }
    }
    if(override_nuc2){
    for (int i = 0; i < 4; i++) {
        if (override_nuc2 & (1 << i)) {
            count[i]++;
        }
    }
    }
    int max_count = *std::max_element(count.begin(), count.end());
    char major_allele = 0;
    for (int i = 0; i < 4; i++) {
        if (count[i] == max_count) {
            major_allele |= (1 << i);
        }
    }
    size_t new_child_size=this_node->children.size();
    if (override_match&&(!override_nuc)) {
        new_child_size--;
    }else if((!override_match)&&override_nuc){
        new_child_size++;
    }
    int new_mutation_count=new_child_size-max_count;
    return std::make_pair(major_allele,new_mutation_count-ori_mutation_count);
}
static void check_mutation(const MAT::Mutation& mutation,
               uint8_t major_alleles_to_check,
               const Mutation_Count_Change_Collection &inserted_count,
               uint8_t major_alleles_ref) {
    assert(major_alleles_ref == major_alleles_to_check);
    nuc_one_hot ori_majority_allele =mutation.get_all_major_allele();
    if (ori_majority_allele == major_alleles_ref) {
        assert(inserted_count.empty()||inserted_count.back().get_position() < mutation.get_position());
    } else {
        assert(!inserted_count.empty());
        const Mutation_Count_Change& last_inserted_count=inserted_count.back();
        assert(last_inserted_count.get_position() == mutation.get_position());
        uint8_t expected_incremented=major_alleles_ref & (~ori_majority_allele);
        assert(last_inserted_count.get_incremented() ==expected_incremented);
        uint8_t expected_decremented=(~major_alleles_ref) & ori_majority_allele;
        assert(last_inserted_count.get_decremented()==expected_decremented);
    }
}
#endif

void test_allele_out_init(
    const MAT::Node *this_node,const MAT::Node *altered_node,
    const MAT::Mutation &mutation,int score_change,
    nuc_one_hot major_alleles_to_check, nuc_one_hot altered_allele,
    const Mutation_Count_Change_Collection &last_inserted_count,
    std::vector<state_change_hist_dbg>& debug
){
    bool is_add_leaf=this_node->is_leaf()&&(altered_node==nullptr);
    int position=mutation.get_position();
    nuc_one_hot par_nuc=mutation.get_par_one_hot();
    #ifdef DEBUG_CHECK_MAJOR_ALLELE_EACH_STEP
    nuc_one_hot second_allele=is_add_leaf?mutation.get_mut_one_hot():(nuc_one_hot)0;
    auto temp =
        get_majority_allele(this_node,altered_node, 0, altered_allele,
                            mutation.get_position(), mutation.get_mut_one_hot(),second_allele);
    nuc_one_hot major_alleles_ref=temp.first;
    check_mutation(mutation, major_alleles_to_check, last_inserted_count, major_alleles_ref);
    if (is_add_leaf) {
       // assert(!mutation.get_tie_one_hot());
        temp.second++;
        /*if (this_node->mutations.find(mutation.get_position())!=this_node->mutations.end()&&mutation.is_valid()) {
            score_change++;
        }*/
    }
    assert(score_change==temp.second);
    #endif
    assert(debug.empty()||debug.back().position<position);
    debug.emplace_back(position,par_nuc,major_alleles_to_check,score_change,(last_inserted_count.empty()||last_inserted_count.back().get_position()!=mutation.get_position())?Mutation_Count_Change(mutation):last_inserted_count.back());

}

int rewind_mutations(int target_position,dbg_iter& debug,dbg_iter& end){
            int rewinded=0;
            while (debug!=end&&debug->position<target_position) {
                nuc_one_hot parent_nuc = debug->par_nuc.back();
                debug->par_nuc.push_back(parent_nuc);
                debug->major_allele.push_back(parent_nuc);
                debug->mutation_score_change.push_back(0);
                debug->count_change.push_back(debug->count_change.back());
                debug->count_change.back().set_change(0, 0,0,true);
                rewinded++;
                debug++;
            }
            return rewinded;
}

void test_allele_count_out(
    const MAT::Node *this_node,
    const MAT::Mutation &mutation,
    nuc_one_hot major_alleles_to_check, nuc_one_hot altered_allele,
    const Mutation_Count_Change_Collection &last_inserted_count, int score_change,
    dbg_iter& debug,dbg_iter& debug_end,
    const std::vector<MAT::Node *>& node_stack) {

    rewind_mutations(mutation.get_position(), debug,debug_end);
    altered_allele = debug->major_allele.back();
    nuc_one_hot par_nuc=mutation.get_par_one_hot();
    #ifdef DEBUG_CHECK_MAJOR_ALLELE_EACH_STEP
    auto temp =
        get_majority_allele(this_node, node_stack.back(), 0, altered_allele,
                            mutation.get_position(),mutation.get_mut_one_hot());
    nuc_one_hot major_alleles_ref=temp.first;
    check_mutation(mutation, major_alleles_to_check, last_inserted_count,
              major_alleles_ref);
    assert(score_change==temp.second);
    #endif
    debug->par_nuc.push_back(par_nuc);
    debug->major_allele.push_back(
        major_alleles_to_check);
    debug->mutation_score_change.push_back(score_change);
    debug->count_change.push_back((last_inserted_count.empty()||last_inserted_count.back().get_position()!=mutation.get_position())?Mutation_Count_Change(mutation):last_inserted_count.back());
    debug++;
}

void test_allele_count_out_LCA(
    const MAT::Node *LCA_node,
    const MAT::Node *src_branch,bool is_src_terminal,
    const MAT::Node *dst_branch,uint8_t dst_terminal_state,
    const MAT::Mutation &mutation,
    nuc_one_hot major_alleles_to_check,int score_change,
    const Mutation_Count_Change_Collection &last_inserted_count,
    std::vector<state_change_hist_dbg>& debug,
    std::vector<LCA_merged_states>::const_iterator& in,const std::vector<LCA_merged_states>::const_iterator& end) {

    int position=mutation.get_position();
    while (in!=end&&in->position<position) {
        assert(debug.empty()||debug.back().position<in->position);
        debug.emplace_back(in->position,in->par_allele,in->par_allele,0,Mutation_Count_Change());
        in++;
    }

    nuc_one_hot altered_allele_src=0;
    nuc_one_hot altered_allele_dst=0;
    int parsimony_score_adjustment=0;
    if (in==end||in->position!=position) {
        assert(in==end||in->position>position);
        if(!(is_src_terminal|dst_terminal_state)){
            return;
        }
        if (!is_src_terminal) {
            src_branch=0;
        }else {
            altered_allele_src=0;
            dst_branch=0;
        }
    }else{
        altered_allele_src=in->major_allele_src;
        altered_allele_dst=in->major_allele_dst;
        if (!altered_allele_dst) {
            dst_branch=nullptr;
        }
        if (!is_src_terminal&&!altered_allele_src) {
            src_branch=nullptr;
            parsimony_score_adjustment=1;
        }
    }
    if (dst_terminal_state) {
        altered_allele_dst=dst_terminal_state;
    }else if(dst_branch&&dst_branch->parent!=LCA_node) {
        altered_allele_dst=0;
    }
    
    nuc_one_hot par_nuc=mutation.get_par_one_hot();
    #ifdef DEBUG_CHECK_MAJOR_ALLELE_EACH_STEP
    auto temp =
        get_majority_allele(LCA_node, src_branch, dst_branch, altered_allele_src,
                            mutation.get_position(), mutation.get_mut_one_hot(),altered_allele_dst);
    if (dst_terminal_state) {
        temp.second++;
    }
    temp.second+=parsimony_score_adjustment;
    nuc_one_hot major_alleles_ref=temp.first;
    check_mutation(mutation, major_alleles_to_check, last_inserted_count,
              major_alleles_ref);
    assert(score_change==temp.second);
    #endif
    assert(debug.empty()||debug.back().position<position);
    debug.emplace_back(position,par_nuc,major_alleles_to_check,score_change,(last_inserted_count.empty()||last_inserted_count.back().get_position()!=mutation.get_position())?Mutation_Count_Change(mutation):last_inserted_count.back());
    if (in!=end&&in->position==position) {
        in++;
    }
}

void extract_last(const std::vector<std::pair<int, std::vector<nuc_one_hot>>>& in,std::vector<std::pair<int, nuc_one_hot>>& out ){
    for( auto e:in){
        out.emplace_back(e.first,e.second.back());
    }
}

static void remove_child(MAT::Node* child_to_remove){
        auto parent=child_to_remove->parent;
        auto iter=std::find(parent->children.begin(),parent->children.end(),child_to_remove);
        assert(iter!=parent->children.end());
        parent->children.erase(iter);
}
static void check_change_with_full_fitch_sankoff(const state_change_hist_dbg& to_check,const std::vector<MAT::Node*>& node_stack, const MAT::Tree& new_tree, const std::vector<uint8_t>& boundary1_major_allele, const std::unordered_map<std::string, std::vector<MAT::Node*>>& mutation_counts_old,const std::vector<std::vector<MAT::Node*>>& mutation_count_ref,MAT::Node* added,std::unordered_set<std::string>& checked_nodes,int& total_difference,std::vector<std::pair<MAT::Node* ,int>>& reported_diffenence){
    const std::vector<nuc_one_hot>& major_alleles_to_check=to_check.major_allele;
    const std::vector<int>& mutation_count_difference_to_check=to_check.mutation_score_change;
    assert(!major_alleles_to_check.empty());
    assert(!mutation_count_difference_to_check.empty());
    auto allele_iter=major_alleles_to_check.begin();
    auto count_diff_iter=mutation_count_difference_to_check.begin();
    auto node_iter=node_stack.begin();
    auto end=node_stack.end();
    while (node_iter!=end) {
        MAT::Node* node=*node_iter;
        const std::string& node_identifier=node->identifier;
        MAT::Node* node_in_new_tree;
        int mutation_count_old=0;
        if (added&&node_iter==node_stack.begin()) {
            node_in_new_tree=added;
        }else{
            node_in_new_tree=new_tree.get_node(node_identifier);
            auto iter=mutation_counts_old.find(node_identifier);
            mutation_count_old=iter==mutation_counts_old.end()?0:iter->second.size();
        }
        int node_idx_in_new_tree=node_in_new_tree->bfs_index;
        nuc_one_hot major_allele_ref=boundary1_major_allele[node_idx_in_new_tree]&0xf;
        int new_mutation_count=mutation_count_ref[node_idx_in_new_tree].size();
        int mutation_count_difference_ref=new_mutation_count-mutation_count_old;
        int mutation_count_difference_to_check=*count_diff_iter;
        if (mutation_count_difference_to_check!=mutation_count_difference_ref||*allele_iter!=major_allele_ref) {
            nuc_one_hot allele_to_check=*allele_iter;
            auto iter=node->mutations.find(to_check.position);
            nuc_one_hot state;
            int idx=node_iter-node_stack.begin();
            if (iter==node->mutations.end()) {
                fputs("\nis par\n", stderr);
                state=get_parent_state(node, to_check.position);
            }else {
                const auto& corresponding_input_change=to_check.count_change[idx];
                if (corresponding_input_change.get_incremented()==0&&corresponding_input_change.get_decremented()==0) {
                    allele_to_check=iter->get_all_major_allele();
                }
            }
            std::vector<std::pair<MAT::Node*, MAT::Mutation*>> ori_mutation_map;
            MAT::Node* replaced_node=0;
            std::string id_replaced=idx?node_stack[idx-1]->identifier:"";
            for(int i=0;i<node->children.size();i++){
                MAT::Node* child=node->children[i];
                if (child->identifier==id_replaced) {
                    replaced_node=child;
                }
                auto iter=child->mutations.find(to_check.position);
                if (iter==child->mutations.end()) {
                }else {
                    ori_mutation_map.emplace_back(child,&(*iter));
                }
            }
            fprintf(stderr,"major_allele is %d, ref %d at position %d\n", (uint8_t)allele_to_check,(uint8_t)major_allele_ref,to_check.position);
        assert(allele_to_check==major_allele_ref);
        assert(mutation_count_difference_to_check==mutation_count_difference_ref);
        }
        total_difference+=mutation_count_difference_to_check;
        reported_diffenence.emplace_back(node_in_new_tree,mutation_count_difference_to_check);
        checked_nodes.emplace(node_in_new_tree->identifier);
        node_iter++;
        allele_iter++;
        count_diff_iter++;
    }
}
static void confirm_no_change(
    const std::vector<uint8_t> &boundary1_major_allele,
    const std::unordered_map<std::string, std::vector<MAT::Node *>>
        &mutation_counts_old,
    const std::vector<std::vector<MAT::Node *>> &mutation_count_ref,
    int position, MAT::Node *node, size_t &node_idx_in_new_tree,
    const std::string &node_identifier,bool is_added=false) {
    nuc_one_hot major_allele_ref =
        boundary1_major_allele[node_idx_in_new_tree] & 0xf;
    int new_mutation_count = mutation_count_ref[node_idx_in_new_tree].size();
    int mutation_count_old = 0;
    if(!is_added){
    auto iter = mutation_counts_old.find(node_identifier);
    if (iter != mutation_counts_old.end()) {
        mutation_count_old = iter->second.size();
    }}

    auto old_mut_iter = node->mutations.find(position);
    nuc_one_hot old_state;
    if (old_mut_iter == node->mutations.end()) {
        old_state = get_parent_state(node, position);
    } else {
        old_state =
            old_mut_iter->get_all_major_allele();
    }

    assert(old_state = major_allele_ref);
    assert(mutation_count_old == new_mutation_count);
}
static void check_no_change_whole_tree(const std::vector<MAT::Node*>& all_nodes, const std::unordered_set<std::string>& checked_nodes,const std::vector<uint8_t> &boundary1_major_allele,
    const std::unordered_map<std::string, std::vector<MAT::Node *>>
        &mutation_counts_old,
    const std::vector<std::vector<MAT::Node *>> &mutation_count_ref,
    int position,MAT::Tree* old_tree){
    for (size_t idx=0; idx<all_nodes.size(); idx++) {
        MAT::Node* node=all_nodes[idx];
        std::string& id=node->identifier;
        if (checked_nodes.count(id)) {
            continue;
        }
        MAT::Node* old_node=old_tree->get_node(id);
        confirm_no_change(boundary1_major_allele, mutation_counts_old, mutation_count_ref, position, old_node, idx, id);
    }
}
static void check_no_change_with_full_fitch_sankoff(
    const std::vector<MAT::Node *> &node_stack, const MAT::Tree &new_tree,
    const std::vector<uint8_t> &boundary1_major_allele,
    const std::unordered_map<std::string, std::vector<MAT::Node *>>
        &mutation_counts_old,
    const std::vector<std::vector<MAT::Node *>> &mutation_count_ref,
    int position, MAT::Node *added) {
    for (auto node : node_stack) {
        size_t node_idx_in_new_tree;
        const std::string &node_identifier = node->identifier;
        if (added && node == node_stack[0]) {
            node_idx_in_new_tree = added->bfs_index;
        } else {
            node_idx_in_new_tree = new_tree.get_node(node_identifier)->bfs_index;
        }
        confirm_no_change(boundary1_major_allele, mutation_counts_old,
                  mutation_count_ref, position, node, node_idx_in_new_tree,
                  node_identifier,added&&node==node_stack[0]);
    }
}
int get_parsimmony_score_dumb(MAT::Node* ancestor,MAT::Node* LCA,MAT::Node* src, MAT::Node* dst,const std::vector<state_change_hist_dbg>& debug_from_src,const std::vector<MAT::Node*> node_stack_from_src,const std::vector<state_change_hist_dbg>& debug_from_dst,const std::vector<MAT::Node*> node_stack_from_dst,const std::vector<state_change_hist_dbg>& debug_above_LCA,const std::vector<MAT::Node*> node_stack_above_LCA){
            MAT::Tree new_tree;
            MAT::Node *new_ancestor = new Mutation_Annotated_Tree::Node(*ancestor,nullptr,&new_tree,false);
            new_tree.curr_internal_node=ancestor->tree->curr_internal_node;
            new_tree.root=new_ancestor;
            MAT::Node *new_src = new_tree.get_node(src->identifier);
            MAT::Node *new_LCA = new_tree.get_node(LCA->identifier);
            MAT::Node *original_src_parent_in_new_tree = new_src->parent;
            MAT::Node *new_dst = new_tree.get_node(dst->identifier);
            remove_child(new_src);
            MAT::Node *added;
            if (dst==LCA) {
                added=new_tree.create_node(std::to_string(++new_tree.curr_internal_node),new_dst);
                MAT::Node* new_src_branch_node=new_tree.get_node(node_stack_from_src.back()->identifier);
                auto iter=std::find(new_LCA->children.begin(),new_LCA->children.end(),new_src_branch_node);
                assert(iter!=new_LCA->children.end());
                new_LCA->children.erase(iter);
                added->children.push_back(new_src_branch_node);
                new_src_branch_node->parent=added;
                added->children.push_back(new_src);
                new_src->parent=added;
            }else {
                MAT::Node* dst_parent=new_dst->parent;
                added=new_tree.create_node(std::to_string(++new_tree.curr_internal_node),dst_parent);
                auto iter=std::find(dst_parent->children.begin(),dst_parent->children.end(),new_dst);
                dst_parent->children.erase(iter);
                added->children.push_back(new_dst);
                new_dst->parent=added;
                added->children.push_back(new_src);
                new_src->parent=added;
            }
            std::vector<MAT::Node*> new_bfs_ordered_nodes=new_tree.breadth_first_expansion();
            std::vector<Mutation_Count> original_mutation_count;
            original_mutation_count.reserve(debug_above_LCA.size());
            for(const auto& mut:debug_above_LCA){
                original_mutation_count.emplace_back(mut.position);
            }
            count_mutations(ancestor, original_mutation_count);
            int parsimony_score_change=0;
            auto from_src_iter=debug_from_src.begin();
            auto from_dst_iter=debug_from_dst.begin();
            auto from_src_end=debug_from_src.end();
            auto from_dst_end=debug_from_dst.end();
            
            for(size_t idx=0;idx<original_mutation_count.size();idx++){
                int position=original_mutation_count[idx].position;
                while (from_src_iter!=from_src_end&&from_src_iter->position<position) {
                    from_src_iter++;
                }
                while (from_dst_iter!=from_dst_end&&from_dst_iter->position<position) {
                    from_dst_iter++;
                }
                std::unordered_set<std::string>  checked_nodes;
                nuc_one_hot LCA_parent_state=get_parent_state(ancestor, position);
                std::vector<uint8_t> boundary1_major_allele(new_bfs_ordered_nodes.size()+8);
                std::vector<uint8_t> boundary2_allele(new_bfs_ordered_nodes.size()+8);
                MAT::Mutation mut(position);
                const auto non_ref_muts=mutated_positions[mut];
                FS_backward_pass(new_bfs_ordered_nodes, boundary1_major_allele,  *non_ref_muts, MAT::Mutation::refs[position]);
                int total_difference=0;
                std::vector<uint8_t> states_out(new_bfs_ordered_nodes.size());
                std::vector<std::pair<MAT::Node* ,int>> reported_diffenence;
                std::vector<std::vector<MAT::Node*>> mutation_count(new_bfs_ordered_nodes.size());
                int new_parsimony_score=FS_forward_assign_states_only(new_bfs_ordered_nodes, boundary1_major_allele, LCA_parent_state, states_out,mutation_count);
                if (from_src_iter!=from_src_end&&from_src_iter->position==position) {
                    check_change_with_full_fitch_sankoff(*from_src_iter, node_stack_from_src, new_tree, boundary1_major_allele,original_mutation_count[idx].mut_count ,mutation_count,0,checked_nodes,total_difference,reported_diffenence);
                }else {
                    check_no_change_with_full_fitch_sankoff(node_stack_from_src, new_tree, boundary1_major_allele, original_mutation_count[idx].mut_count ,mutation_count, position,0);
                }
                if (from_dst_iter!=from_dst_end&&from_dst_iter->position==position) {
                    check_change_with_full_fitch_sankoff(*from_dst_iter, node_stack_from_dst, new_tree, boundary1_major_allele,original_mutation_count[idx].mut_count ,mutation_count,added,checked_nodes,total_difference,reported_diffenence);
                }else {
                    check_no_change_with_full_fitch_sankoff(node_stack_from_dst, new_tree, boundary1_major_allele, original_mutation_count[idx].mut_count ,mutation_count, position,added);
                }
                std::vector<MAT::Node*> have_change;
                
                /*for (int i=0; i<new_bfs_ordered_nodes.size(); i++) {
                    if (mutation_count[i].size()) {
                        have_change.push_back(new_bfs_ordered_nodes[i]);
                    }
                }*/
                check_change_with_full_fitch_sankoff(debug_above_LCA[idx], node_stack_above_LCA, new_tree, boundary1_major_allele,original_mutation_count[idx].mut_count ,mutation_count,dst==LCA?added:nullptr,checked_nodes,total_difference,reported_diffenence);

                int ref_difference=new_parsimony_score-original_mutation_count[idx].count;
                if(ref_difference!=total_difference){
                    check_no_change_whole_tree(new_bfs_ordered_nodes, checked_nodes, boundary1_major_allele, original_mutation_count[idx].mut_count, mutation_count, position, dst->tree);
                    for (int i=0; i<new_bfs_ordered_nodes.size(); i++) {
                            for(auto n:mutation_count[i]){
                                fprintf(stderr,"Mutation at idx %zu, contibuted by %zu\n",new_bfs_ordered_nodes[i]->bfs_index,n->bfs_index);
                            }
                    }
                    assert(false);
                }
                parsimony_score_change+=(new_parsimony_score-original_mutation_count[idx].count);
            }
            new_tree.delete_nodes();
            return parsimony_score_change;
}

void prep_LCA_checker(
    const std::vector<state_change_hist_dbg> &from_src,
    const std::vector<state_change_hist_dbg> &from_dst,
    std::vector<LCA_merged_states> &out) {
    
    #define src_ele from_src[from_src_idx]
    #define dst_ele from_dst[from_dst_idx]
    #define src_pos (src_ele.position)
    #define dst_pos (dst_ele.position)
    #define dst_bound_check (from_dst_idx < from_dst.size())
    #define push_out(ele,src_or_dst) \
    if((!out.empty())&&ele.position==out.back().position){\
        assert(ele.par_nuc.back()==out.back().par_allele);\
    }else {\
        assert(out.empty()||ele.position>out.back().position);\
        out.emplace_back(ele.position,ele.par_nuc.back());\
    }\
    out.back().src_or_dst=ele.major_allele.back();
    #define push_src push_out(src_ele,major_allele_src)
    #define push_dst push_out(dst_ele,major_allele_dst)
    
    size_t from_dst_idx = 0;
    for (size_t from_src_idx = 0; from_src_idx<from_src.size();from_src_idx++ ) {
        while ( dst_bound_check &&dst_pos <=src_pos) {
            push_dst;
            from_dst_idx++;
        }
        push_src;

    }
    while (dst_bound_check) {
            push_dst;
            from_dst_idx++;
    }
    #undef src_ele
    #undef dst_ele
    #undef src_pos
    #undef dst_pos
    #undef dst_bound_check
    #undef push_out
}
#endif