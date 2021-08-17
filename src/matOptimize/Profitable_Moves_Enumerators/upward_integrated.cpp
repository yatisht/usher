#include <algorithm>
#include <climits>
#include <cstdlib>
#include "process_individual_mutation.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
void check_parsimony_score_change_above_LCA(MAT::Node *LCA, int &parsimony_score_change,
        Mutation_Count_Change_Collection &parent_added,
        const std::vector<MAT::Node *> &node_stack_from_src,
        std::vector<MAT::Node *> &node_stack_above_LCA,
        Mutation_Count_Change_Collection &parent_of_parent_added,
        MAT::Node *ancestor) ;
void output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,
                   const std::vector<MAT::Node *> &node_stack_from_src,
                   std::vector<MAT::Node *> &node_stack_from_dst,
                   std::vector<MAT::Node *> &node_stack_above_LCA,int radius_left);
struct Split_LCA_outputter{
    MAT::Node* src;
    const std::vector<MAT::Node *> &node_stack_from_src;
    output_t &output;
    int radius;
    void operator()(Mutation_Count_Change_Collection& allele_count_change_from_splitting_LCA, int par_score_change,MAT::Node* LCA){
        Mutation_Count_Change_Collection parent_of_parent_added;
        parent_of_parent_added.reserve(allele_count_change_from_splitting_LCA.size());
        std::vector<MAT::Node *> node_stack_above_LCA;
        check_parsimony_score_change_above_LCA(LCA, par_score_change, allele_count_change_from_splitting_LCA, node_stack_from_src,
                                           node_stack_above_LCA, parent_of_parent_added, LCA);
        std::vector<MAT::Node*> ignored;
        output_result(src, LCA, LCA, par_score_change, output,
                  node_stack_from_src, ignored, node_stack_above_LCA,radius);
    }
};
/*
1. accumulate mutation vector of src to what is needed if placed as children of parent
2. whether it is profitable to just place src as children of parent node
3. change to major allele set at this node
            parent node
            /           \   
            ? (2.split)   1. mutation vector if placed here
          /
        this node (3. change to major allele set, also parsimony score change)
        /
    ...
    /
    src
*/
template<typename T>
class src_op{};
template<>
class src_op<MAT::Mutations_Collection>{
    MAT::Mutations_Collection::const_iterator iter;
    MAT::Mutations_Collection::const_iterator end;
    public:
    src_op(const MAT::Mutations_Collection& in):iter(in.begin()),end(in.end()){}
    /**
     * Process allele count change up to positon
     * @param[in] position For position and makor allele set
     * @param[out] mut_count_change_out majar allele set change output
     * @param[inout] parsimony score change for mutations preceding position
     * @return 
     */
    nuc_one_hot next_pos(const MAT::Mutation& parent_variable,Mutation_Count_Change_Collection& mut_count_change_out,int& par_score_change){
        while (iter!=end&&iter->get_position()<parent_variable.get_position()) {
            if (iter->is_valid()) {
                par_score_change--;
            }
            iter++;
        }
        if (iter!=end&&iter->get_position()==parent_variable.get_position()) {
            Mutation_Count_Change temp(*iter,iter->get_all_major_allele(), 0);
            //same as above
            par_score_change--;
            iter++;
            return decrement_mutation_count(
                mut_count_change_out, parent_variable, temp,
                par_score_change);
        }
        Mutation_Count_Change temp(parent_variable,parent_variable.get_mut_one_hot(),0);
        //decrement_mutation_count manipulate parsimony score change only on how such change change the number of children able to have major allele state (as if src branch node is not removed), so need to decrement artificailly to account for src->parent have one less children
        par_score_change--;
        return decrement_mutation_count(mut_count_change_out, parent_variable, temp,par_score_change);
    }
    void finish(Mutation_Count_Change_Collection& mut_count_change_out,int& par_score_change){
        while (iter!=end) {
            if(iter->is_valid()){
                par_score_change--;
            }
            iter++;
        }
        //add sentinel
        mut_count_change_out.emplace_back();
    }
};
template<>
class src_op<Mutation_Count_Change_Collection>{
    Mutation_Count_Change_Collection::const_iterator iter;
    Mutation_Count_Change_Collection::const_iterator end;
    public:
    src_op(const Mutation_Count_Change_Collection& in):iter(in.begin()),end(in.end()){}
    nuc_one_hot next_pos(const MAT::Mutation& parent_variable,Mutation_Count_Change_Collection& mut_count_change_out,int& par_score_change){
        while (iter->get_position()<parent_variable.get_position()) {
            iter++;
        }
        if (iter->get_position()==parent_variable.get_position()) {
            return decrement_increment_mutation_count(parent_variable, *iter, mut_count_change_out, par_score_change);
        }
        return parent_variable.get_all_major_allele();
    }
    void finish(Mutation_Count_Change_Collection& mut_count_change_out,int& par_score_change){
        //add sentinel
        mut_count_change_out.emplace_back();
    }
};
template< typename T> 
class Mut_Merger{};
template<>
class Mut_Merger<MAT::Mutations_Collection>{
    MAT::Mutations_Collection::const_iterator iter;
    MAT::Mutations_Collection::const_iterator end;
    Mutation_Count_Change_Collection& mutations_out;
    Mut_Merger(const MAT::Mutations_Collection& in,MAT::Node* this_node,int children_par_score,Mutation_Count_Change_Collection& mutations_out):iter(in.begin()),end(in.end()),mutations_out(mutations_out){
        mutations_out.reserve(in.size()+this_node->mutations.size());
    }
    void operator()(const MAT::Mutation& mut){
        if (!mut.is_valid()) {
            return;
        }
        while (iter!=end&&iter->get_position()<mut.get_position()) {
            mutations_out.emplace_back(*iter,0,iter->get_all_major_allele());
            iter++;
        }
        if (iter!=end&&iter->get_position()==mut.get_position()) {
            auto new_par_nuc=mut.get_par_one_hot();
            if (new_par_nuc!=iter->get_all_major_allele()) {
                mutations_out.emplace_back(*iter,0,iter->get_all_major_allele());
                mutations_out.back().set_par_nuc(new_par_nuc);
            }
            iter++;
        }else {
            mutations_out.emplace_back(mut,0,mut.get_mut_one_hot());
        }
    }
    void finish(Split_LCA_outputter& output){
        while (iter!=end) {
            mutations_out.emplace_back(*iter,0,iter->get_all_major_allele());
            iter++;
        }
        //add sentinel
        mutations_out.emplace_back();
    }
};
template<>
class Mut_Merger<Mutation_Count_Change_Collection>{
    Mutation_Count_Change_Collection::const_iterator iter;
    Mutation_Count_Change_Collection::const_iterator end;
    bool have_not_shared;
    int par_score_change_split_LCA;
    MAT::Node* LCA;
    Mutation_Count_Change_Collection allele_count_change_out;
    Mutation_Count_Change_Collection& mutations_out;
    Mut_Merger(const Mutation_Count_Change_Collection& in,MAT::Node* this_node,int children_par_score,Mutation_Count_Change_Collection& mutations_out):iter(in.begin()),end(in.end()),have_not_shared(false),par_score_change_split_LCA(children_par_score),LCA(this_node),mutations_out(mutations_out){
        allele_count_change_out.reserve(in.size()+this_node->mutations.size());
        mutations_out.reserve(in.size()+this_node->mutations.size());
    }
    void operator()(const MAT::Mutation& mut,nuc_one_hot major_allele){
        while (iter->get_position()<mut.get_position()) {
            mutations_out.emplace_back(*iter);
            have_not_shared|=add_node_split(*iter,allele_count_change_out, par_score_change_split_LCA);
            iter++;
        }
        if (iter->get_position()==mut.get_position()) {
            //accumulate mutation for placement
            if (mut.is_valid()) {
                auto new_par_nuc=mut.get_par_one_hot();
                if (new_par_nuc!=iter->get_incremented()) {
                    mutations_out.emplace_back(*iter);
                    mutations_out.back().set_par_nuc(new_par_nuc);
                }
            }
            //split LCA
            have_not_shared|=add_node_split(mut,iter->get_incremented(),major_allele, allele_count_change_out, par_score_change_split_LCA);
            iter++;
        }else {
            if (mut.is_valid()) {
                mutations_out.emplace_back(mut,0,mut.get_mut_one_hot());
            }
            have_not_shared|=add_node_split(mut,allele_count_change_out, par_score_change_split_LCA);
        }
    }
    void finish(Split_LCA_outputter& output){
        while (iter!=end) {
            allele_count_change_out.emplace_back(*iter);
            have_not_shared|=add_node_split(*iter, allele_count_change_out, par_score_change_split_LCA);
            iter++;
        }
        //sentinel included
        if (have_not_shared) {
            output(allele_count_change_out,par_score_change_split_LCA,LCA);
        }
    }
};

template<typename T>
int upward_integrated(const MAT::Node* this_node,const T& src_mut_under_this,const T& src_branch_change,Mutation_Count_Change_Collection& this_allele_change_out,Mutation_Count_Change_Collection& src_under_parent,Split_LCA_outputter& output,int child_par_score_change){
    int par_score_change=0;
    src_op<T> prev_count_change(src_branch_change);
    this_allele_change_out.reserve(this_node->mutations.size()+src_branch_change.size());
    Mut_Merger<T> merger(src_mut_under_this,this_node,child_par_score_change,src_under_parent);
    for( const auto& mut:this_node->mutations){
        auto major_allele=prev_count_change.next_pos(mut,this_allele_change_out,par_score_change);
        merger(mut,major_allele);
    }
    prev_count_change.finish(this_allele_change_out,par_score_change);
    merger.finalize(output);
}