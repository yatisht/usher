#include "ripples_util.hpp"

inline MAT::Node *get_node_cstr(MAT::Tree &tree, char *name) {
    return tree.get_node(std::string(name));
}

Ripple_Result_Pack* Ripple_Pipeline::operator()(MAT::Node* node_to_consider) const {
    fprintf(stderr, "At node id: %s\n",
            node_to_consider->identifier.c_str());

    int orig_parsimony = (int)node_to_consider->mutations.size();

    Pruned_Sample pruned_sample(node_to_consider);
    // Find mutations on the node to prune
    auto node_to_root = T.rsearch(node_to_consider->identifier, true);
    for (auto curr : node_to_root) {
        for (auto m : curr->mutations) {
            pruned_sample.add_mutation(m);
        }
    }
    auto parsimony_threshold=orig_parsimony - parsimony_improvement;
    if (parsimony_threshold<0) {
        return (new Ripple_Result_Pack{std::vector<Recomb_Interval>(),node_to_consider,orig_parsimony});
    }
    //==== new mapper
    Ripples_Mapper_Output_Interface mapper_out;
    ripples_mapper(pruned_sample, mapper_out, nodes_to_seach.size(),index_map,do_parallel, traversal_track,tree_height,T.root,node_to_consider);
    //==== END new mapper
    tbb::concurrent_vector<Recomb_Interval> valid_pairs_con;
    ripplrs_merger(pruned_sample, index_map,nodes_to_seach, nodes_to_seach.size(),
                   parsimony_threshold, T,
                   valid_pairs_con, mapper_out, num_threads, branch_len,
                   min_range, max_range);
    std::vector<Recomb_Interval> temp(std::vector<Recomb_Interval>(valid_pairs_con.begin(),valid_pairs_con.end()));
    std::sort(temp.begin(),temp.end(),interval_sorter());
    return (new Ripple_Result_Pack{combine_intervals(temp),node_to_consider,orig_parsimony});
}

void Ripple_Finalizer::operator()(Ripple_Result_Pack* result) const {
    // print combined pairs
    auto & valid_pairs=result->intervals;
    auto node_to_consider=result->node_to_consider;
    auto orig_parsimony=result->orig_parsimony;
    for (auto p : valid_pairs) {
        std::string end_range_high_str =
            (p.end_range_high == 1e9) ? "GENOME_SIZE"
            : std::to_string(p.end_range_high);
        auto donor_adj_parsimony=p.d.node_parsimony+!p.d.is_sibling;
        auto acceptor_adj_parsimony=p.a.node_parsimony+!p.a.is_sibling;
        fprintf(
            recomb_file,
            "%s\t(%i,%i)\t(%i,%s)\t%s\t%c\t%i\t%s\t%c\t%i\t%i\t%i\t%i\n",
            node_to_consider->identifier.c_str(), p.start_range_low,
            p.start_range_high, p.end_range_low, end_range_high_str.c_str(),
            p.d.node->identifier.c_str(), p.d.is_sibling?'y':'n',
            donor_adj_parsimony, p.a.node->identifier.c_str(),
            p.a.is_sibling?'y':'n', acceptor_adj_parsimony, orig_parsimony,
            std::min(
        {orig_parsimony, donor_adj_parsimony, acceptor_adj_parsimony}),
        p.d.parsimony + p.a.parsimony);
        fflush(recomb_file);
    }

    if (!valid_pairs.empty()) {
        fprintf(desc_file, "%s\t", node_to_consider->identifier.c_str());
        for (auto l : T.get_leaves(node_to_consider->identifier)) {
            fprintf(desc_file, "%s,", l->identifier.c_str());
        }
        fprintf(desc_file, "\n");
        fflush(desc_file);
        fprintf(stderr, "Done %zu/%zu branches [RECOMBINATION FOUND!]\n\n",
                ++num_done, total_size);
    } else {
        fprintf(stderr, "Done %zu/%zu branches\n\n", ++num_done,
                total_size);
    }
    delete result;
}
