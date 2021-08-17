#include "Profitable_Moves_Enumerators.hpp"
#include "process_each_node.hpp"
#include "split_node_helpers.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
static void
add_remaining_dst_to_LCA_nodes(MAT::Node *cur, const MAT::Node *LCA,
                               std::vector<MAT::Node *> &dst_stack) {
    while (cur != LCA) {
        // assert(dst_stack.empty()||dst_stack.back()!=cur);
        dst_stack.push_back(cur);
        cur = cur->parent;
    }
}

void output_result(MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
                   int parsimony_score_change, output_t &output,
                   const std::vector<MAT::Node *> &node_stack_from_src,
                   std::vector<MAT::Node *> &node_stack_from_dst,
                   std::vector<MAT::Node *> &node_stack_above_LCA,
                   int radius_left);
struct Not_LCA_Outputer {
    int par_score_from_src_removal;
    MAT::Node *LCA;
    MAT::Node *src;
    const Mutation_Count_Change_Collection &allele_count_change_from_src;
    const std::vector<MAT::Node *> &node_stack_from_src;
    output_t &out;
    int src_side_max_improvement;
    int radius;
    void operator()(Mutation_Count_Change_Collection &parent_added,
                    MAT::Node *dst_node, int par_score_change) {
        auto parsimony_score_change =
            par_score_change + par_score_from_src_removal;
        std::vector<MAT::Node *> node_stack_from_dst;
        Mutation_Count_Change_Collection parent_of_parent_added;
        parent_of_parent_added.reserve(parent_added.size());
        node_stack_from_dst.push_back(dst_node);
        auto this_node = dst_node;
        while (this_node != LCA) {
            parent_of_parent_added.clear();
            get_intermediate_nodes_mutations(
                this_node, node_stack_from_dst.back(), parent_added,
                parent_of_parent_added, parsimony_score_change);
            node_stack_from_dst.push_back(this_node);
            parent_added.swap(parent_of_parent_added);
            if (parent_added.empty()) {
                add_remaining_dst_to_LCA_nodes(this_node->parent, LCA,
                                               node_stack_from_dst);
                break;
            }
            this_node = this_node->parent;
            auto parsimony_score_lower_bound =
                parsimony_score_change - src_side_max_improvement;
            if (parsimony_score_lower_bound > 0) {
                for (const auto &mut_count_change : parent_added) {
                    if (mut_count_change.get_incremented()) {
                        parsimony_score_lower_bound--;
                    }
                }
                if (parsimony_score_lower_bound > 0) {
                    return;
                }
            }
        }
        std::vector<MAT::Node *> node_stack_above_LCA;
        parent_of_parent_added.reserve(parent_added.size() +
                                       allele_count_change_from_src.size());
        // Adjust LCA node and above

        bool is_src_terminal = src->parent == LCA;
        if ((!(allele_count_change_from_src.empty() && parent_added.empty())) ||
            is_src_terminal) {
            get_LCA_mutation(
                LCA, is_src_terminal ? src : node_stack_from_src.back(),
                is_src_terminal, allele_count_change_from_src, parent_added,
                parent_of_parent_added, parsimony_score_change);
        }
        node_stack_above_LCA.push_back(LCA);
        parent_added.swap(parent_of_parent_added);
        check_parsimony_score_change_above_LCA(
            LCA, parsimony_score_change, parent_added, node_stack_from_src,
            node_stack_above_LCA, parent_of_parent_added, LCA->parent);
        output_result(src, dst_node, LCA, parsimony_score_change, out,
                      node_stack_from_src, node_stack_from_dst,
                      node_stack_above_LCA, radius);
    }
};
template <typename T> class src_op {};
template <> class src_op<MAT::Mutations_Collection> {
    MAT::Mutations_Collection::const_iterator iter;
    MAT::Mutations_Collection::const_iterator end;
    Mutation_Count_Change_Collection &mut_change_out;

  public:
    src_op(const MAT::Mutations_Collection &in,
           Mutation_Count_Change_Collection &out,MAT::Node* this_node)
        : iter(in.begin()), end(in.end()), mut_change_out(out) {
              mut_change_out.reserve(this_node->mutations.size()+in.size());
          }
    void operator()(const MAT::Mutation &mut) {
        if (!mut.is_valid()) {
            return;
        }
        while (iter != end && iter->get_position() < mut.get_position()) {
            mut_change_out.emplace_back(*iter, 0, iter->get_all_major_allele());
            iter++;
        }
        if (iter != end && iter->get_position() == mut.get_position()) {
            mut_change_out.emplace_back(*iter, 0, iter->get_all_major_allele());
            mut_change_out.back().set_par_nuc(mut.get_mut_one_hot());
            iter++;
        } else {
            mut_change_out.emplace_back(mut, 0, mut.get_par_one_hot());
            mut_change_out.back().set_par_nuc(mut.get_mut_one_hot());
        }
    }
    void finalize(Not_LCA_Outputer& out) {
        while (iter != end) {
            mut_change_out.emplace_back(*iter, 0, iter->get_all_major_allele());
            iter++;
        }
        // add sentinel
        mut_change_out.emplace_back();
    }
};
template <> class src_op<Mutation_Count_Change_Collection> {
    Mutation_Count_Change_Collection::const_iterator iter;
    Mutation_Count_Change_Collection::const_iterator end;
    Mutation_Count_Change_Collection &mut_change_out;
    Mutation_Count_Change_Collection allele_count_change_if_add;
    int par_score_change_from_add;
    bool have_not_shared;
    MAT::Node* this_node;
  public:
    src_op(const Mutation_Count_Change_Collection &in,
           Mutation_Count_Change_Collection &out,MAT::Node* this_node)
        : iter(in.begin()), end(in.end()), mut_change_out(out),
          par_score_change_from_add(0), have_not_shared(false),this_node(this_node) {
              mut_change_out.reserve(this_node->mutations.size()+in.size());
              allele_count_change_if_add.reserve(this_node->mutations.size()+in.size());
          }
    void operator()(const MAT::Mutation &mut) {
        while (iter->get_position() < mut.get_position()) {
            have_not_shared |= add_node_split(*iter, allele_count_change_if_add,
                                              par_score_change_from_add);
            mut_change_out.emplace_back(*iter);
            iter++;
        }
        if (iter->get_position() == mut.get_position()) {
            mut_change_out.emplace_back(*iter);
            mut_change_out.back().set_par_nuc(mut.get_mut_one_hot());
            have_not_shared |= add_node_split(
                mut, mut.get_all_major_allele(), iter->get_incremented(),
                allele_count_change_if_add, par_score_change_from_add);
            iter++;
        } else {
            if (mut.is_valid()) {
                mut_change_out.emplace_back(mut, 0, mut.get_par_one_hot());
                mut_change_out.back().set_par_nuc(mut.get_mut_one_hot());
            }
            have_not_shared |= add_node_split(mut, allele_count_change_if_add,
                                              par_score_change_from_add);
        }
    }
    void finalize(Not_LCA_Outputer& out) {
        while (iter != end) {
            mut_change_out.emplace_back(*iter);
            have_not_shared |= add_node_split(*iter, allele_count_change_if_add,
                                              par_score_change_from_add);
            iter++;
        }
        out(allele_count_change_if_add,this_node,par_score_change_from_add);
    }
};
template<typename T>
void downward_integrated(
    MAT::Node *this_node, const T &from_parent,
    Mutation_Count_Change_Collection &mut_count_change_out,Not_LCA_Outputer& out) {
        src_op<T> op(from_parent,mut_count_change_out,this_node);
        for (const auto& mut : this_node->mutations) {
            op(mut);
        }
        op.finalize(out);
}