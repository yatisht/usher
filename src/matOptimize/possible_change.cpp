#include "Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#include <utility>
#include <vector>

#include "Profitable_Moves_Enumerators/process_individual_mutation.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"

std::vector<std::array<std::vector<node_info>,4>> addable_idxes;
namespace MAT = Mutation_Annotated_Tree;
short default_decrement_effect[4];
short default_increment_effect[4];
struct initer {
    initer(){
    for (int nuc_idx = 0; nuc_idx < 4; nuc_idx++) {
        default_increment_effect[nuc_idx] = 0;
        default_decrement_effect[nuc_idx] = 0;
        for (int ptr_idx = 0; ptr_idx < 16; ptr_idx++) {
            if (ptr_idx & (1 << nuc_idx)) {
                default_increment_effect[nuc_idx] |= (1 << ptr_idx);
                default_decrement_effect[nuc_idx] |= (1 << ptr_idx);
            }
        }
    }
    }
};
initer init;

struct Sensitive_Alleles {
    int postion;
    uint16_t decrement_effect;
    uint16_t increment_effect;
};

bool decrement_sensitive(uint8_t decrement, uint8_t major_allele,
                         uint8_t boundary_allele, int parent_increment_effect,
                         int parent_decrement_effect) {
    auto not_decremented = major_allele & (~decrement);
    if (not_decremented) {
        // No Par score change at this position, but did remove major alleles
        // that will increase parsimony score of ancestors?
        uint8_t decremented_major_allele = major_allele & decrement;
        return parent_decrement_effect & (1 << decremented_major_allele);
    }
    // decremented all major allele. see whether some new major alllel from
    // boundary alleles can offset
    return !(parent_increment_effect & (1 << (boundary_allele & (~decrement))));
}

bool increment_sensitive(uint8_t increment, uint8_t major_allele,
                         uint8_t boundary_allele, int parent_increment_effect,
                         int parent_decrement_effect) {
    auto major_allele_incremented = increment & major_allele;
    if (major_allele_incremented) {
        // Reduced parsimony score at this node, make sure not decrease it at
        // parent node from removal of major allele not decremented
        return !(parent_decrement_effect &
                 (1 << (major_allele & (~increment))));
    }
    // Incrementing boundary allele, no par score change at this node, but may
    // parent
    return parent_increment_effect & (1 << (boundary_allele & increment));
}
uint8_t move_bits(uint32_t in) {
    return (3 & (in >> 1)) | (4 & (in >> 2)) | (8 & (in >> 5));
}
#define ROT(x, t) ((x << t) | (x >> (8 - t)))
static std::pair<uint16_t, uint16_t>
update_sensitve_allele(int parent_decrement_effect, int parent_increment_effect,
                       MAT::Mutation &out) {
    // Effect of decrementing alleles
    // only decrementing major allele is interesting
    auto major_allele = out.get_all_major_allele();
    auto boundary_allele = out.get_boundary1_one_hot();
    uint32_t increment_effect = (1 << 15);
    uint32_t decrement_effect = (1 << 15);
    for (uint32_t ptr = 1; ptr < 15; ptr++) {
        increment_effect |=
            (increment_sensitive(ptr, major_allele, boundary_allele,
                                 parent_increment_effect,
                                 parent_decrement_effect)
             << ptr);
    }
    for (uint32_t ptr = 1; ptr < 15; ptr++) {
        decrement_effect |=
            (decrement_sensitive(ptr, major_allele, boundary_allele,
                                 parent_increment_effect,
                                 parent_decrement_effect)
             << ptr);
    }
    out.set_sensitive_change(move_bits(decrement_effect),
                             move_bits(increment_effect));
    return std::make_pair(decrement_effect, increment_effect);
}
static std::pair<uint16_t, uint16_t>
update_sensitve_allele(const Sensitive_Alleles &par_mut, MAT::Mutation &out) {
    return update_sensitve_allele(par_mut.decrement_effect,
                                  par_mut.increment_effect, out);
}
static std::pair<uint16_t, uint16_t>
update_sensitve_allele(MAT::Mutation &out) {
    auto par_idx = one_hot_to_two_bit(out.get_par_one_hot());
    return update_sensitve_allele(default_decrement_effect[par_idx],
                                  default_increment_effect[par_idx], out);
}
static void filter_output(const MAT::Mutation &mut,
                          std::pair<uint16_t, uint16_t> effect,
                          std::vector<Sensitive_Alleles> &out) {
    auto mut_idx = one_hot_to_two_bit(mut.get_mut_one_hot());
    if (effect.first != default_decrement_effect[mut_idx] ||
        effect.second != default_increment_effect[mut_idx]) {
        out.push_back(
            Sensitive_Alleles{mut.get_position(), effect.first, effect.second});
    }
}

typedef std::vector<std::array<tbb::concurrent_vector<node_info>, 4>> pos_tree_t;
struct Walker : public tbb::task {
    std::vector<Sensitive_Alleles> sensitive_locus;
    MAT::Node *root;
    pos_tree_t& pos_tree;
    // Sensitive locus is for this node
    Walker(MAT::Node *root,pos_tree_t& pos_tree) : root(root),pos_tree(pos_tree) {}
    void register_change(const MAT::Node* node, const MAT::Mutation& mut){
        auto new_alleles=mut.get_sensitive_increment()&(~mut.get_par_one_hot());
        if (new_alleles) {
            for (int idx=0; idx<4; idx++) {
                if (new_alleles&(1<<idx)) {
                    pos_tree[mut.get_position()][idx].emplace_back(node_info{node->dfs_index,node->level});
                }
            }
        }
    }
    void register_change(const Sensitive_Alleles& last, uint16_t new_increment,int position,const MAT::Node* node){
        auto new_alleles=new_increment&(~last.increment_effect);
        if (new_alleles) {
            for (int idx=0; idx<4; idx++) {
                if (new_alleles&(1<<(two_bit_to_one_hot(idx)))) {
                    pos_tree[position][idx].emplace_back(node_info{node->dfs_index,node->level});
                }
            }
        }
    }
    void merge(MAT::Node *to_set, std::vector<Sensitive_Alleles> *output) {
        auto iter = sensitive_locus.begin();
        std::pair<uint16_t, uint16_t> effect;
        for (auto &mut : to_set->mutations) {
            while (iter->postion < mut.get_position()) {
                iter++;
            }
            if (mut.get_position() == iter->postion) {
                effect = update_sensitve_allele(*iter, mut);
                register_change(*iter,effect.second,mut.get_position(),to_set);
            } else {
                effect = update_sensitve_allele(mut);
                register_change(to_set,mut);
            }
            if (output) {
                filter_output(mut, effect, *output);
            }
        }
        if (output) {
            output->push_back(Sensitive_Alleles{INT_MAX});
        }
    }
    tbb::task *execute() override {
        auto continuation = new (allocate_continuation()) tbb::empty_task;
        std::vector<tbb::task *> tasks;
        tasks.reserve(root->children.size());
        for (auto child : root->children) {
            if (child->is_leaf()) {
                merge(child, nullptr);
            } else {

                auto child_task =
                    new (continuation->allocate_child()) Walker(child,pos_tree);
                merge(child, &child_task->sensitive_locus);
                tasks.push_back(child_task);
            }
        }
        continuation->set_ref_count(tasks.size());
        for (auto task : tasks) {
            continuation->spawn(*task);
        }
        return tasks.empty() ? continuation : nullptr;
    }
};
void output_addable_idxes(pos_tree_t& in){
    addable_idxes=std::vector<std::array<std::vector<node_info>,4>>(MAT::Mutation::refs.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,MAT::Mutation::refs.size()),[&in](tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
            for (int nu_idx=0; nu_idx<4; nu_idx++) {
                addable_idxes[idx][nu_idx]=std::vector<node_info>(in[idx][nu_idx].begin(),in[idx][nu_idx].end());
                std::sort(addable_idxes[idx][nu_idx].begin(),addable_idxes[idx][nu_idx].end());
            }
        }
    });
}
void adjust_all(MAT::Tree &tree) {
    fprintf(stderr, "start\n");
    auto start = std::chrono::steady_clock::now();
    pos_tree_t pos_tree(MAT::Mutation::refs.size());
    auto task_root = new (tbb::task::allocate_root()) Walker(tree.root,pos_tree);
    for (auto &mut : tree.root->mutations) {
        auto effect = update_sensitve_allele(mut);
        filter_output(mut, effect, task_root->sensitive_locus);
    }
    task_root->sensitive_locus.push_back(Sensitive_Alleles{INT_MAX});
    tbb::task::spawn_root_and_wait(*task_root);
    output_addable_idxes(pos_tree);
    size_t max_change=0;
    size_t total=0;
    for (const auto& pos_nuc : addable_idxes) {
        for (int nuc_idx=0; nuc_idx<4; nuc_idx++) {
            auto this_size=pos_nuc[nuc_idx].size();
            max_change=std::max(max_change,this_size);
            total+=this_size;
        }
    }
    fprintf(stderr, "Total %zu,max %zu \n",total,max_change);
    fprintf(stderr, "Updating sensitive alleles take %ld sec",
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start)
                .count());
}
#ifdef TEST
struct change_log {
    MAT::Mutation *relevant_mut;
    int score_change;
    MAT::Node *node;
};
static void check(MAT::Node *start_node, nuc_one_hot incremented,
                  const MAT::Mutation &mut) {
    Mutation_Count_Change mut_count_change(mut, 0, incremented);
    std::vector<change_log> log;
    int change_so_far = 0;
    int expected_change = mut.get_sensitive_increment() & incremented ? -1 : 0;
    while (start_node) {
        auto iter = start_node->mutations.find(mut.get_position());
        if (iter == start_node->mutations.end()) {
            auto change = mut_count_change.get_default_change_internal();
            log.emplace_back(change_log{nullptr, change, start_node});
            change_so_far += change;
            break;
        } else {
            Mutation_Count_Change_Collection ignored;
            int change = 0;
            decrement_increment_mutation_count(*iter, mut_count_change, ignored,
                                               change);
            log.emplace_back(change_log{&(*iter), change, start_node});
            change_so_far += change;
            if (ignored.empty()) {
                break;
            } else {
                mut_count_change = ignored[0];
            }
        }
        start_node = start_node->parent;
    }

    assert(expected_change == change_so_far);
}
static void check_each_node(MAT::Node *node) {
    for (const auto &mut : node->mutations) {
        uint8_t nuc_to_check =
            mut.get_all_major_allele() | mut.get_boundary1_one_hot();
        for (int idx = 0; idx < 4; idx++) {
            uint8_t temp = 1 << idx;
            if (temp & nuc_to_check) {
                check(node, temp, mut);
            }
        }
    }
}
struct checker {
    std::vector<MAT::Node *> to_check;
    void operator()(const tbb::blocked_range<size_t> &range) const {
        for (size_t idx = range.begin(); idx < range.end(); idx++) {
            check_each_node(to_check[idx]);
        }
    }
};

int main(int argc, char **argv) {
    MAT::Tree tree;
    tree.load_detatiled_mutations(argv[1]);
    adjust_all(tree);
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
                      checker{bfs_ordered_nodes});
}
#endif