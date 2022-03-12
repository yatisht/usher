#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <src/matOptimize/mutation_annotated_tree.hpp>
#include <array>
#include <unordered_map>
#include <vector>
class Min_Back_Mut_FS_Score_PerNode_T{
    std::array<unsigned int, 16> last_mutation;
    public:
    std::array<unsigned int, 4> parsimony_score;
    std::array<unsigned int, 4> back_mutation_count;
    Min_Back_Mut_FS_Score_PerNode_T(){
        parsimony_score.fill(0);
        back_mutation_count.fill(0);
        last_mutation.fill(0);
    }
    unsigned int& get_last_mut(uint8_t surface_nuc,uint8_t nuc_beneath){
        return last_mutation[4*surface_nuc+nuc_beneath];
    }
};
class Min_Back_Mut_FS_Score_PerNode_Choice{
    uint8_t choices;
    public:
    Min_Back_Mut_FS_Score_PerNode_Choice():choices(0){}
    uint8_t get_choice(uint8_t par_nuc_2bit) const{
        return (choices>>(2*par_nuc_2bit))&3;
    }
    void set_choice(uint8_t par_nuc_2bit, uint8_t child_nuc_2bit){
        /*if (par_nuc_2bit>3||child_nuc_2bit>3) {
            fprintf(stderr, "Not 2but\n");
            raise(SIGTRAP);
        }*/
        choices|=(child_nuc_2bit<<(2*par_nuc_2bit));
        /*if (get_choice(par_nuc_2bit)!=child_nuc_2bit) {
            raise(SIGTRAP);
        }*/
    }
    uint8_t set_leaf_choice(uint8_t leaf_nuc_one_hot){
        uint8_t default_nuc=__builtin_ctz(leaf_nuc_one_hot);
        for (int nuc_idx=0; nuc_idx<4; nuc_idx++) {
            if (leaf_nuc_one_hot&(1<<nuc_idx)) {
                set_choice(nuc_idx, nuc_idx);
            }else {
                set_choice(nuc_idx, default_nuc);
            }
        }
        return default_nuc;
    }
};
namespace MAT = Mutation_Annotated_Tree;
typedef std::vector<std::pair<long, nuc_one_hot>> mutated_t;
void find_best_child_nuc(Min_Back_Mut_FS_Score_PerNode_T &child_score_output,
               int this_nuc_idx, int &best_child_nuc,
               int &best_parsimony_score, int &best_back_mutation_count) {
    for (int child_nuc = 0; child_nuc < 4; child_nuc++) {
        int parsimony;
        int back_mutation_count;
        if (child_nuc == this_nuc_idx) {
            // child nucleotide is the same as the current node
            // No parsimony score, back mutation increase
            parsimony = child_score_output.parsimony_score[child_nuc];
            back_mutation_count =
                child_score_output.back_mutation_count[child_nuc];
        } else {
            // different
            // Parsimony score increase
            parsimony = child_score_output.parsimony_score[child_nuc] + 1;
            // back mutation count increase by this_nuc->child_nuc->....
            // ->child_nuc->this_nuc
            back_mutation_count =
                child_score_output.back_mutation_count[child_nuc] +
                child_score_output.get_last_mut(child_nuc, this_nuc_idx);
        }
        if (parsimony < best_parsimony_score ||
            ((parsimony == best_parsimony_score) &&
             (back_mutation_count < best_back_mutation_count))) {
            best_parsimony_score = parsimony;
            best_child_nuc = child_nuc;
            best_back_mutation_count = back_mutation_count;
        }
    }
}
static Min_Back_Mut_FS_Score_PerNode_T
backward_pass(MAT::Node *node, mutated_t::const_iterator &iter,
              std::vector<Min_Back_Mut_FS_Score_PerNode_Choice> &state_output,
              uint8_t ref_nuc) {
    Min_Back_Mut_FS_Score_PerNode_T output;
    for (const auto child : node->children) {
        auto child_idx=child->dfs_index;
        auto& child_state_output=state_output[child->dfs_index];
        if (child->is_leaf()) {
            auto child_nuc=ref_nuc;
            while (iter->first<child_idx) {
                iter++;
            }
            if (iter->first==child_idx) {
                child_nuc=iter->second;
            }
            auto default_nuc=child_state_output.set_leaf_choice(child_nuc);
            for (int nuc_idx=0; nuc_idx<4; nuc_idx++) {
                if (!(child_nuc&(1<<nuc_idx))) {
                    output.parsimony_score[nuc_idx]++;
                    output.get_last_mut(nuc_idx, default_nuc)++;
                }
            }
        }else {
            auto child_score_output=backward_pass(child, iter, state_output, ref_nuc);
            for (int this_nuc_idx=0; this_nuc_idx<4; this_nuc_idx++) {
                int best_child_nuc=0;
                int best_parsimony_score=INT_MAX;
                int best_back_mutation_count=INT_MAX;
                find_best_child_nuc(child_score_output, this_nuc_idx, best_child_nuc,
                          best_parsimony_score, best_back_mutation_count);
                //Commit changes
                output.parsimony_score[this_nuc_idx]+=best_parsimony_score;
                output.back_mutation_count[this_nuc_idx]+=best_back_mutation_count;
                child_state_output.set_choice(this_nuc_idx, best_child_nuc);
                if (best_child_nuc==this_nuc_idx) {
                    //No mutation at this node so last mutation visible is transparent
                    for (int beneath_nuc_idx=0; beneath_nuc_idx<4; beneath_nuc_idx++) {
                        output.get_last_mut(this_nuc_idx, beneath_nuc_idx)+=child_score_output.get_last_mut(this_nuc_idx, beneath_nuc_idx);
                    }
                }else {
                    output.get_last_mut(this_nuc_idx, best_child_nuc)++;
                }
            }
        }
    }
    return output;
}
static void forward_pass(MAT::Node* node, const std::vector<Min_Back_Mut_FS_Score_PerNode_Choice> &state_output,uint8_t par_nuc,
    const MAT::Mutation& mut_template,std::vector<std::vector<MAT::Mutation>>& mutation_output,int& mutations_added){
    auto this_idx=node->dfs_index;
    auto this_nuc=state_output[this_idx].get_choice(par_nuc);
    if (this_nuc!=par_nuc) {
        MAT::Mutation mut_copy(mut_template);
        mut_copy.set_par_mut(1<<par_nuc, 1<<this_nuc);
        mut_copy.set_auxillary(1<<this_nuc, 0);
        if (!(mut_copy.get_mut_one_hot()&&mut_copy.get_par_one_hot())) {
            fprintf(stderr, "Not setting par mut, par_mut :%d, mut_nuc %d \n",par_nuc,this_nuc);
            raise(SIGTRAP);
        }
        mutation_output[this_idx].push_back(mut_copy);
        mutations_added++;
    }
    for (auto child : node->children) {
        forward_pass(child,state_output,this_nuc,mut_template,mutation_output,mutations_added);
    }
}

void Min_Back_Fitch_Sankoff(MAT::Node* root_node,const MAT::Mutation& mut_template,
    std::vector<std::vector<MAT::Mutation>>& mutation_output,mutated_t& positions,size_t dfs_size){
    mutated_t::const_iterator mut_iter=positions.begin();
    std::vector<Min_Back_Mut_FS_Score_PerNode_Choice> best_nucleotide(dfs_size);
    auto root_scores=backward_pass(root_node, mut_iter,best_nucleotide, mut_template.get_ref_one_hot());
    int root_best_nuc=0;
    int root_best_par=INT_MAX;
    int root_best_back=INT_MAX;
    auto ref_nuc_two_bit=__builtin_ctz(mut_template.get_ref_one_hot());
    find_best_child_nuc(root_scores,ref_nuc_two_bit , root_best_nuc, root_best_par, root_best_back);
    //fprintf(stderr, "parsimony at %d:%d\n",mut_template.get_position(),root_best_par);
    best_nucleotide[0].set_choice(ref_nuc_two_bit, root_best_nuc);
    int mutations_added=0;
    forward_pass(root_node,best_nucleotide,ref_nuc_two_bit,mut_template,mutation_output,mutations_added);
    if (mutations_added!=root_best_par) {
        fprintf(stderr, "Mutation mismatch %d backward vs %d after forward\n",root_best_par,mutations_added);
        raise(SIGTRAP);
    }
}

static int count_back_mutations_helper(const MAT::Node* root, std::unordered_map<int, MAT::Mutation> parent_mutations){
    int back_mutation_count=0;
    for (const auto& mut : root->mutations) {
        auto emplace_result=parent_mutations.emplace(mut.get_position(),mut);
        if (!emplace_result.second) {
            if (mut.get_mut_one_hot()==emplace_result.first->second.get_par_one_hot()) {
                back_mutation_count++;
            }
            emplace_result.first->second=mut;
        }
    }
    for (const auto child : root->children) {
        back_mutation_count+=count_back_mutations_helper(child, parent_mutations);
    }
    return back_mutation_count;
}
int count_back_mutation(const MAT::Tree& tree){
    return count_back_mutations_helper(tree.root, std::unordered_map<int, MAT::Mutation>());
}
