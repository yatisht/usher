#ifndef PROFITABLE_MOVES_ENUMERATOR
#define PROFITABLE_MOVES_ENUMERATOR
#include "../tree_rearrangement_internal.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;

class Mutation_Count_Change {
    int position;
    uint8_t chromIdx;
    nuc_one_hot decremented_allele;
    nuc_one_hot incremented_allele;
    bool was_valid;
    nuc_one_hot ori_state;
    nuc_one_hot par_state;
    nuc_one_hot new_state;

  public:
  static const char VALID_MASK=1;
  static const char END_MASK=2;
  Mutation_Count_Change(){
              decremented_allele=0;
        incremented_allele=0;
  }
    Mutation_Count_Change(const MAT::Mutation &pos) {
        position = pos.get_position();
        chromIdx = pos.get_chromIdx();
        was_valid=pos.is_valid();
        ori_state=pos.get_all_major_allele();
        par_state=pos.get_par_one_hot();
        decremented_allele=0;
        incremented_allele=0;
    }
    operator bool() const{
        return was_valid;
    }
    void set_valid(bool valid){
        was_valid=valid;
    }
    int get_position() const { return position; }
    nuc_one_hot get_decremented() const { assert(decremented_allele!=0xff); return decremented_allele; }
    nuc_one_hot get_incremented() const { assert(incremented_allele!=0xff); return incremented_allele; }
    nuc_one_hot get_ori_state() const{return ori_state;}
    void set_change(nuc_one_hot decremented, nuc_one_hot incremented,nuc_one_hot new_state) {
        decremented_allele = decremented;
        incremented_allele = incremented;
        this->new_state=new_state;
        assert(new_state=(ori_state|incremented)&(~decremented));
    }
    nuc_one_hot get_new_state() const{
        return new_state;
    }
    nuc_one_hot get_par_state() const {
        return par_state;
    }
    int get_default_change_internal()const {
        if(par_state&incremented_allele){
            assert(was_valid);
            return -1;
        }
        if(par_state&decremented_allele){
            assert(!was_valid);
            return 1;
        }
        return 0;
    }
    int get_default_change_terminal()const {
        if(incremented_allele&&(!(par_state&incremented_allele))){
            return 1;
        }
        if((!(decremented_allele&par_state))&&decremented_allele){
            return -1;
        }
        return 0;
    }
    int get_default_change(const Mutation_Count_Change& other)const {
        nuc_one_hot incremented=incremented_allele&other.incremented_allele;
        nuc_one_hot decremented=decremented_allele&other.decremented_allele;
        assert(position==other.position);
        assert(par_state==other.par_state);
        if(par_state&incremented){
            assert(was_valid);
            return -1;
        }
        if(par_state&decremented){
            assert(!was_valid);
            return 1;
        }
        return 0;
    }
    bool operator<(const Mutation_Count_Change &rhs) const {
        return position < rhs.position;
    }
    void set_par_nuc(nuc_one_hot par_nuc){
        par_state=par_nuc;
    }
};

static bool operator<(const MAT::Mutation &lhs, const Mutation_Count_Change &rhs) {
    return lhs.get_position() < rhs.get_position();
}

typedef std::vector<Mutation_Count_Change> Mutation_Count_Change_Collection;

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
struct state_change_hist_dbg{
    int position;
    std::vector<nuc_one_hot> par_nuc;
    std::vector<nuc_one_hot> major_allele;
    std::vector<int> mutation_score_change;
    std::vector<Mutation_Count_Change> count_change;
    state_change_hist_dbg(){}
    state_change_hist_dbg(int position,nuc_one_hot par_nuc,nuc_one_hot major_allele,int mutation_score_change,const Mutation_Count_Change& count_change):position(position),par_nuc({par_nuc}),major_allele({major_allele}),mutation_score_change({mutation_score_change}),count_change({count_change}){}
};
typedef std::vector<state_change_hist_dbg>::iterator dbg_iter;
struct LCA_merged_states{
    int position;
    nuc_one_hot major_allele_src;
    nuc_one_hot major_allele_dst;
    nuc_one_hot par_allele;
    LCA_merged_states(){}
    LCA_merged_states(int position,nuc_one_hot par_allele):position(position),major_allele_src(0),major_allele_dst(0),par_allele(par_allele){}
};
void test_allele_out_init(
    const MAT::Node *this_node,const MAT::Node *altered_node,
    const MAT::Mutation &mutation,int score_change,
    nuc_one_hot major_alleles_to_check, nuc_one_hot altered_allele,
    const Mutation_Count_Change_Collection &last_inserted_count,
    std::vector<state_change_hist_dbg>& debug
);

void test_allele_count_out(
    const MAT::Node *this_node,
    const MAT::Mutation &mutation,
    nuc_one_hot major_alleles_to_check, nuc_one_hot altered_allele,
    const Mutation_Count_Change_Collection &last_inserted_count, int score_change,
    dbg_iter& debug,dbg_iter& debug_end,
    const std::vector<MAT::Node *>& node_stack);

void test_allele_count_out_LCA(
    const MAT::Node *LCA_node,
    const MAT::Node *src_branch,bool is_src_terminal,
    const MAT::Node *dst_branch,uint8_t dst_terminal_state,
    const MAT::Mutation &mutation,
    nuc_one_hot major_alleles_to_check,int score_change,
    const Mutation_Count_Change_Collection &last_inserted_count,
    std::vector<state_change_hist_dbg>& debug,
    std::vector<LCA_merged_states>::const_iterator& in,const std::vector<LCA_merged_states>::const_iterator& end);

int get_parsimmony_score_dumb(MAT::Node* LCA,MAT::Node* src, MAT::Node* dst,const std::vector<state_change_hist_dbg>& debug_from_src,const std::vector<MAT::Node*> node_stack_from_src,const std::vector<state_change_hist_dbg>& debug_from_dst,const std::vector<MAT::Node*> node_stack_from_dst,const std::vector<state_change_hist_dbg>& debug_above_LCA,const std::vector<MAT::Node*> node_stack_above_LCA);

void prep_LCA_checker(
    const std::vector<state_change_hist_dbg> &from_src,
    const std::vector<state_change_hist_dbg> &from_dst,
    std::vector<LCA_merged_states> &out);

int rewind_mutations(int target_position,dbg_iter& debug,dbg_iter& end);
void update_par_cur_nuc(MAT::Mutations_Collection::const_iterator parent_mutation_iter,dbg_iter& debug_iter,dbg_iter& debug_end);
void update_dbg_vector_score_only(
    int position, dbg_iter &debug_iter,dbg_iter &debug_end,
    int par_score_change);

void set_LCA_par_score_change(
    const Mutation_Count_Change &change,
    std::vector<state_change_hist_dbg> &debug,
    std::vector<LCA_merged_states>::const_iterator &in,
    const std::vector<LCA_merged_states>::const_iterator &begin,
    const std::vector<LCA_merged_states>::const_iterator &end,
    int score_change);
#endif

/*
struct New_Tie_t {
    int position;
    uint8_t tied_alleles;
    uint8_t ori_par_allele;
    // short score_change;
};
typedef std::vector<New_Tie_t> New_Tie_Collection_t;


void
register_new_tied_mutations(New_Tie_Collection_t &new_unresolved_ties,
                            nuc_one_hot major_alleles, const MAT::Mutation &mut,
                            int parsimony_score_change) ;
                            
                             void resolve_remaining_ties(
    New_Tie_Collection_t::const_iterator &parent_unresolved_iter,
    const New_Tie_Collection_t::const_iterator &parent_unresolved_end,
    int &parsimony_score_change) ;
    
     void
resolve_ties(New_Tie_Collection_t::const_iterator &parent_unresolved_iter,
             const New_Tie_Collection_t::const_iterator &parent_unresolved_end,
             int this_mutation_position, nuc_one_hot major_allele,
             New_Tie_Collection_t &remain_unresolved,
             int &parsimony_score_change);
*/
#endif