#include "mutation_annotated_tree.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <stack>
#include <queue>
#include <tbb/parallel_sort.h>
#include <tbb/task.h>
#include <utility>
#include <vector>
#include <fstream>
// Uses one-hot encoding if base is unambiguous
void check_leaves(const Mutation_Annotated_Tree::Tree& T);
// A:1,C:2,G:4,T:8
using Mutation_Annotated_Tree::Node;
tbb::concurrent_unordered_map<std::string, uint8_t>  Mutation_Annotated_Tree::Mutation::chromosome_map;
std::vector<std::string>  Mutation_Annotated_Tree::Mutation::chromosomes;
std::mutex Mutation_Annotated_Tree::Mutation::ref_lock;
std::vector<nuc_one_hot> Mutation_Annotated_Tree::Mutation::refs;



/* === Tree === */
std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;
    Node* node;
    if (nid == "") {
        if (root == NULL) {
            return traversal;
        }
        node=root;
    }else {
        node=get_node(nid);
    }
    size_t idx=0;
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        curr_node->bfs_index=idx++;
        traversal.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }

    return traversal;
}

static void depth_first_expansion_helper(Mutation_Annotated_Tree::Node* node, std::vector<Mutation_Annotated_Tree::Node*>& vec, size_t& index,size_t level) {
#ifdef DETAIL_DEBUG_NO_LOOP
    assert(std::find(vec.begin(),vec.end(),node)==vec.end());
#endif
    vec.push_back(node);
    node->level=level;
    //assert(vec.size()-1==index);
    node->dfs_index=index;
    index++;
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec,index,level+1);
    }
    node->dfs_end_index=index-1;
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::depth_first_expansion(Mutation_Annotated_Tree::Node* node) const {
    TIMEIT();
    std::vector<Node*> traversal;
    if (node == NULL) {
        node = root;
    }
    size_t index=0;
    if (node == NULL) {
        return traversal;
    }
    depth_first_expansion_helper(node, traversal,index,0);
    return traversal;
}

size_t Mutation_Annotated_Tree::Tree::get_parsimony_score() {
    size_t score = 0;
    auto dfs = depth_first_expansion();
    for (auto n: dfs) {
        for (const auto& mut:n->mutations) {
            score+=mut.is_valid();
        }
    }
    return score;
}
static size_t level_helper(const Node* node) {
    size_t level = 0;
    for (auto child : node->children) {
        level=std::max(level,level_helper(child));
    }
    return level+1;
}
size_t Mutation_Annotated_Tree::Tree::get_max_level() {
    max_level=level_helper(root);
    return max_level;
}
void Mutation_Annotated_Tree::Tree::uncondense_leaves() {
    for (auto cn = condensed_nodes.begin(); cn != condensed_nodes.end(); cn++) {

        auto n = all_nodes[cn->first];
        auto par = (n->parent != NULL) ? n->parent : n;

        size_t num_samples = cn->second.size();

        if (num_samples > 0) {
            rename_node(n->node_id, cn->second[0]);
        }

        for (size_t s = 1; s < num_samples; s++) {
            auto new_node=create_node(cn->second[s]);
            par->add_child(new_node);
        }
    }
    condensed_nodes.clear();
}

void Node::delete_this() {
    for(Node* n:children) {
        n->delete_this();
    }
    delete this;
}
void Mutation_Annotated_Tree::Tree::delete_nodes() {
    if (!root) {
        return;
    }
    root->delete_this();
    root=nullptr;
}
static void get_leaves_helper(const Node* root, std::vector<Node*>& out) {
    for(auto child:root->children) {
        if (child->is_leaf()) {
            out.push_back(child);
        } else {
            get_leaves_helper(child, out);
        }
    }
}
std::vector<Node*> Mutation_Annotated_Tree::Tree::get_leaves(const Node* root) const {
    std::vector<Node*> out;
    if (root==nullptr) {
        root=this->root;
    }
    get_leaves_helper(root, out);
    return out;
}
Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::Tree::copy_tree(){
    Mutation_Annotated_Tree::Tree out;
    out.node_names=node_names;
    out.node_name_to_idx_map=node_name_to_idx_map;
    out.node_idx=node_idx;
    out.max_level=max_level;
    out.condensed_nodes=condensed_nodes;
    out.curr_internal_node=curr_internal_node;
    out.root=new Mutation_Annotated_Tree::Node(*root,nullptr,&out);
    return out;
}
std::string Mutation_Annotated_Tree::Tree::get_clade_assignment (const Node* n, int clade_id, bool include_self) const {
    const Node* anc=n;
    while (anc) {
        if (include_self) {
            if ((int)anc->clade_annotations.size() > clade_id && anc->clade_annotations[clade_id] != "") {
                return anc->clade_annotations[clade_id];
            }
        }
        include_self=true;
        auto curr_node=anc;
        auto par_node=anc->parent;
        while (curr_node!=get_node(curr_node->node_id)) {
            curr_node=get_node(curr_node->node_id);
            par_node=curr_node->parent;
        }
        anc=par_node;
    }
    return "UNDEFINED";
}
void Mutation_Annotated_Tree::Tree::populate_ignored_range() {
    auto leaves=breadth_first_expansion();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,leaves.size()),[&leaves](const tbb::blocked_range<size_t>& range) {
        for (auto idx=range.begin(); idx<range.end(); idx++) {
            auto leaf=leaves[idx];
            leaf->populate_ignored_range();
        }
    });
    fprintf(stderr, "populated ignored range\n");
}
void Mutation_Annotated_Tree::get_sample_mutation_paths(
    Mutation_Annotated_Tree::Tree *T, std::vector<Node *> samples,
    std::string mutation_paths_filename) {
    FILE *mutation_paths_file = fopen(mutation_paths_filename.c_str(), "w");

    for (size_t s = 0; s < samples.size(); s++) {
        auto sample_node = samples[s];

        // Stack for last-in first-out ordering
        std::stack<std::string> mutation_stack;
        std::string curr_node_mutation_string;

        // Mutations on the added sample
        auto curr_node_mutations = sample_node->mutations;
        if (curr_node_mutations.size() > 0) {
            curr_node_mutation_string =
                T->get_node_name(sample_node->node_id) + ":";
            size_t num_mutations = curr_node_mutations.size();
            for (size_t k = 0; k < num_mutations; k++) {
                curr_node_mutation_string +=
                    curr_node_mutations[k].get_string();
                if (k < num_mutations - 1) {
                    curr_node_mutation_string += ',';
                } else {
                    curr_node_mutation_string += ' ';
                }
            }
            mutation_stack.push(curr_node_mutation_string);
        }
        auto anc_node = sample_node;
        // Mutations on the ancestors of added sample
        while (anc_node) {
            curr_node_mutations = anc_node->mutations;
            if (curr_node_mutations.size() > 0) {
                curr_node_mutation_string =
                    T->get_node_name(anc_node->node_id) + ":";
                size_t num_mutations = curr_node_mutations.size();
                for (size_t k = 0; k < num_mutations; k++) {
                    curr_node_mutation_string +=
                        curr_node_mutations[k].get_string();
                    if (k < num_mutations - 1) {
                        curr_node_mutation_string += ',';
                    } else {
                        curr_node_mutation_string += ' ';
                    }
                }
                mutation_stack.push(curr_node_mutation_string);
            }
            anc_node = anc_node->parent;
        }

        fprintf(mutation_paths_file, "%s\t",
                T->get_node_name(sample_node->node_id).c_str());
        while (mutation_stack.size()) {
            fprintf(mutation_paths_file, "%s", mutation_stack.top().c_str());
            mutation_stack.pop();
        }
        fprintf(mutation_paths_file, "\n");
    }
    fclose(mutation_paths_file);
}
static bool is_anncestor_after_dfs(const Node *ancestral, const Node *to_test) {
    if (to_test->dfs_index > ancestral->dfs_index &&
        to_test->dfs_index <= ancestral->dfs_end_index) {
        return true;
    }
    return false;
}
static std::vector<const Node *> get_ancestor(const Node *node, bool include_self) {
    std::vector<const Node *> output;
    if (include_self) {
        output.push_back(node);
    }
    node = node->parent;
    while (node) {
        output.push_back(node);
        node = node->parent;
    }
    return output;
}
Mutation_Annotated_Tree::Tree
Mutation_Annotated_Tree::get_subtree(const Mutation_Annotated_Tree::Tree &tree,
                                     const std::vector<Node *> &samples,
                                     bool keep_clade_annotations) {
    TIMEIT();
    Tree subtree;

    // Set of leaf and internal nodes corresponding to the subtree
    tbb::concurrent_unordered_set<size_t> subtree_nodes;
    // Maintain a set of all ancestors of a sample for each sample
    std::vector<tbb::concurrent_unordered_set<const Node *>> all_ancestors(
        samples.size());

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t k = r.begin(); k < r.end(); ++k) {
                subtree_nodes.insert(samples[k]->node_id);
                auto anc = samples[k];
                while (anc) {
                    all_ancestors[k].insert(anc);
                    anc = anc->parent;
                }
            }
        },
        ap);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, samples.size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                for (size_t j = i + 1; j < samples.size(); ++j) {
                    auto anc = samples[i];
                    while (anc) {
                        if (all_ancestors[j].find(anc) !=
                            all_ancestors[j].end()) {
                            subtree_nodes.insert(anc->node_id);
                            break;
                        }
                        anc = anc->parent;
                    }
                }
            }
        },
        ap);

    auto dfs = tree.depth_first_expansion();
    size_t num_annotations = 0;
    if (keep_clade_annotations) {
        num_annotations = tree.get_num_annotations();
    }
    check_leaves(tree);

    std::stack<Node *> last_subtree_node;
    for (const Mutation_Annotated_Tree::Node* n : dfs) {
        // If the node is in subtree_nodes, it should be added to the subtree
        if (subtree_nodes.find(n->node_id) != subtree_nodes.end()) {
            Node *subtree_parent = NULL;
            if (last_subtree_node.size() > 0) {
                while (!is_anncestor_after_dfs(last_subtree_node.top(), n)) {
                    last_subtree_node.pop();
                }
                subtree_parent = last_subtree_node.top();
            }
            // Add as root of the subtree
            Node *new_node;
            if (subtree_parent == NULL) {
                // for root node, need to size the annotations vector
                new_node = subtree.create_node(
                    tree.get_node_name(n->node_id), num_annotations);
                subtree.root=new_node;
                // need to assign any clade annotations which would belong to
                // that root as well
                for (size_t k = 0; k < num_annotations; k++) {
                    if (n->clade_annotations[k] != "") {
                        new_node->clade_annotations[k] =
                            n->clade_annotations[k];
                    }
                }
                std::vector<const Node *> root_to_node = get_ancestor(n, true);
                std::reverse(root_to_node.begin(), root_to_node.end());
                // root_to_node.emplace_back(n);

                for (auto curr : root_to_node) {
                    for (auto m : curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            // Add to the parent identified
            else {
                new_node = subtree.create_node(
                    tree.get_node_name(n->node_id), num_annotations);
                subtree_parent->children.push_back(new_node);
                new_node->parent = subtree_parent;
                auto par_to_node = get_ancestor(n, true);
                std::reverse(par_to_node.begin(), par_to_node.end());
                par_to_node.erase(par_to_node.begin(),
                                  std::find(par_to_node.begin(),
                                            par_to_node.end(), tree.get_node(subtree_parent->node_id)) +
                                      1);

                for (auto curr : par_to_node) {
                    for (size_t k = 0; k < num_annotations; k++) {
                        if (curr->clade_annotations[k] != "") {
                            new_node->clade_annotations[k] =
                                curr->clade_annotations[k];
                        }
                    }
                    for (auto m : curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            last_subtree_node.push(new_node);
        }
    }
    check_leaves(tree);

    subtree.curr_internal_node = tree.curr_internal_node;

    return subtree;
}
void Mutation_Annotated_Tree::Tree::rotate_for_display(bool reverse) {
    auto dfs = depth_first_expansion();

    std::unordered_map<Node *, int> num_desc;

    for (int i = int(dfs.size()) - 1; i >= 0; i--) {
        auto n = dfs[i];
        int desc = 1;
        for (auto child : n->children) {
            desc += num_desc[child];
        }
        num_desc[n] = desc;
    }

    for (auto n : dfs) {
        if (reverse) {
            tbb::parallel_sort(n->children.begin(), n->children.end(),
                               [&num_desc](Node *n1, Node *n2) {
                                   return num_desc[n1] < num_desc[n2];
                               });
        } else {
            tbb::parallel_sort(n->children.begin(), n->children.end(),
                               [&num_desc](Node *n1, Node *n2) {
                                   return num_desc[n1] > num_desc[n2];
                               });
        }
    }
}
void Mutation_Annotated_Tree::get_random_single_subtree(
    const Mutation_Annotated_Tree::Tree & T, std::vector<Node *> samples,
    std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx,
    bool retain_original_branch_len) {
    // timer.Start();
    std::string preid = "/";
    if (use_tree_idx) {
        preid = "/tree-" + std::to_string(tree_idx) + "-";
    }
    std::set<Node *> leaves_to_keep_set(samples.begin(), samples.end());

    auto all_leaves = T.get_leaves();
    for (size_t i = 0; i < all_leaves.size(); i++) {
        auto l = all_leaves.begin();
        std::advance(l, std::rand() % all_leaves.size());
        leaves_to_keep_set.insert(*l);
        if (leaves_to_keep_set.size() >= subtree_size + samples.size()) {
            break;
        }
    }

    std::vector<Node *> leaves_to_keep(samples.begin(), samples.end());

    auto new_T = get_subtree(T, leaves_to_keep);
    check_leaves(T);
    // Rotate tree for display
    new_T.rotate_for_display();

    // Write subtree to file
    auto subtree_filename = outdir + preid + "single-subtree.nh";
    fprintf(
        stderr,
        "Writing single subtree with %zu randomly added leaves to file %s.\n",
        subtree_size, subtree_filename.c_str());

    std::ofstream subtree_file(subtree_filename.c_str(), std::ofstream::out);
    std::stringstream newick_ss;
    new_T.write_newick_string(newick_ss, new_T.root, true, true,
                              retain_original_branch_len);
    subtree_file << newick_ss.rdbuf();
    subtree_file.close();

    // Write list of mutations on the subtree to file
    auto subtree_mutations_filename =
        outdir + preid + "single-subtree-mutations.txt";
    fprintf(stderr,
            "Writing list of mutations at the nodes of the single subtree to "
            "file %s\n",
            subtree_mutations_filename.c_str());
    FILE *subtree_mutations_file =
        fopen(subtree_mutations_filename.c_str(), "w");

    for (auto n : new_T.depth_first_expansion()) {
        size_t tot_mutations = n->mutations.size();
        fprintf(subtree_mutations_file,
                "%s: ", new_T.get_node_name(n->node_id).c_str());
        for (size_t idx = 0; idx < tot_mutations; idx++) {
            auto m = n->mutations[idx];
            fprintf(subtree_mutations_file, "%s", m.get_string().c_str());
            if (idx + 1 < tot_mutations) {
                fprintf(subtree_mutations_file, ",");
            }
        }
        fprintf(subtree_mutations_file, "\n");
    }

    fclose(subtree_mutations_file);
    check_leaves(T);

    // Expand internal nodes that are condensed
    bool has_condensed = false;
    FILE *subtree_expanded_file = NULL;
    for (auto l : new_T.get_leaves()) {
        if (T.condensed_nodes.find(l->node_id) != T.condensed_nodes.end()) {
            if (!has_condensed) {

                auto subtree_expanded_filename =
                    outdir + preid + "single-subtree-expanded.txt";
                fprintf(stderr,
                        "Subtree has condensed nodes.\nExpanding the condensed "
                        "nodes for the single subtree in file %s\n",
                        subtree_expanded_filename.c_str());
                subtree_expanded_file =
                    fopen(subtree_expanded_filename.c_str(), "w");
                has_condensed = true;
            }
            fprintf(subtree_expanded_file,
                    "%s: ", new_T.get_node_name(l->node_id).c_str());
            auto iter=T.condensed_nodes.find(l->node_id);
            for (auto n : iter->second) {
                fprintf(subtree_expanded_file, "%s ", n.c_str());
            }
            fprintf(subtree_expanded_file, "\n");
        }
    }
    check_leaves(T);
    if (has_condensed) {
        fclose(subtree_expanded_file);
    }
}
void Mutation_Annotated_Tree::get_random_sample_subtrees (const Mutation_Annotated_Tree::Tree& T, std::vector<Node*> samples, std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len) {
    fprintf(stderr, "Computing subtrees for %ld samples. \n\n", samples.size());
    std::string preid = "/";
    if (use_tree_idx) {
        preid = "/tree-" + std::to_string(tree_idx) + "-";
    }
    T.depth_first_expansion();
    // For each final tree, write a subtree of user-specified size around
    // each requested sample in newick format
    // We split the subtree size into two: one half for nearest sequences to
    // the samples and the other half randomly sampled
    //size_t random_subtree_size = print_subtrees_size/2;
    //size_t nearest_subtree_size = print_subtrees_size - random_subtree_size;

    size_t random_subtree_size = subtree_size/5;
    size_t nearest_subtree_size = subtree_size - random_subtree_size;

    //Set a constant random seed
    std::srand(0);

    // Randomly shuffle the leaves for selecting the random subtree
    auto all_leaves = T.get_leaves();
    std::set<Mutation_Annotated_Tree::Node*> random_ordered_leaves;
    for (size_t i=0; i< all_leaves.size(); i++) {
        auto l = all_leaves.begin();
        std::advance(l, std::rand() % all_leaves.size());
        random_ordered_leaves.insert(*l);
        if (random_ordered_leaves.size() >= subtree_size) {
            break;
        }
    }

    // Bool vector to mark which newly placed samples have already been
    // displayed in a subtree (initialized to false)
    std::vector<bool> displayed_samples (samples.size(), false);

    // If the missing sample is not found in the tree, it was not placed
    // because of max_uncertainty. Mark those samples as already
    // displayed.

    int num_subtrees = 0;
    for (size_t i = 0; i < samples.size(); i++) {

        if (displayed_samples[i]) {
            continue;
        }

        Mutation_Annotated_Tree::Node* last_anc = samples[i];
        std::vector<Node*> leaves_to_keep;

        // Keep moving up the tree till a subtree of required size is
        // found
        auto anc=samples[i];
        while (anc) {
            size_t num_leaves = anc->get_num_leaves();
            if (num_leaves < subtree_size) {
                last_anc = anc;
                continue;
            }

            if (num_leaves > subtree_size) {
                struct NodeDist {
                    Mutation_Annotated_Tree::Node* node;
                    uint32_t num_mut;

                    NodeDist(Node* n, uint32_t d) {
                        node = n;
                        num_mut = d;
                    }

                    inline bool operator< (const NodeDist& n) const {
                        return ((*this).num_mut < n.num_mut);
                    }
                };

                for (auto l: T.get_leaves(last_anc)) {
                    leaves_to_keep.emplace_back(l);
                }

                std::vector<NodeDist> node_distances;
                for (auto l: T.get_leaves(anc)) {
                    if (is_anncestor_after_dfs(last_anc, l)) {
                        continue;
                    }

                    uint32_t dist = 0;
                    for (auto a: get_ancestor(l, true)) {
                        if (a == anc) {
                            break;
                        }
                        dist += a->mutations.size();
                    }

                    node_distances.emplace_back(NodeDist(l, dist));
                }

                std::sort(node_distances.begin(), node_distances.end());
                for (auto n: node_distances) {
                    if (leaves_to_keep.size() >= nearest_subtree_size) {
                        break;
                    }
                    leaves_to_keep.emplace_back(n.node);
                }

                if ((nearest_subtree_size < subtree_size) && (nearest_subtree_size < node_distances.size())) {
                    std::vector<NodeDist> remaining_node_distances = {node_distances.begin()+nearest_subtree_size, node_distances.end()};
                    std::shuffle(remaining_node_distances.begin(), remaining_node_distances.end(), std::default_random_engine {});

                    for (auto n: remaining_node_distances) {
                        if (leaves_to_keep.size() == subtree_size) {
                            break;
                        }
                        leaves_to_keep.emplace_back(n.node);
                    }
                }
            } else {
                for (auto l: T.get_leaves(anc)) {
                    if (leaves_to_keep.size() == subtree_size) {
                        break;
                    }
                    leaves_to_keep.emplace_back(l);
                }
            }

            auto new_T = Mutation_Annotated_Tree::get_subtree(T, leaves_to_keep);

            // Rotate tree for display
            new_T.rotate_for_display();

            tbb::parallel_for (tbb::blocked_range<size_t>(i+1, samples.size(), 100),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t j=r.begin(); j<r.end(); ++j) {
                    if (!displayed_samples[j]) {
                        if (new_T.get_node(samples[j]->node_id) != NULL) {
                            displayed_samples[j] = true;
                        }
                    }
                }
            });

            // Write subtree to file
            ++num_subtrees;
            auto subtree_filename = outdir + preid + "subtree-" + std::to_string(num_subtrees) + ".nh";
            fprintf(stderr, "Writing subtree %d to file %s.\n", num_subtrees, subtree_filename.c_str());
            //FILE* subtree_file = fopen(subtree_filename.c_str(), "w");
            //fprintf(subtree_file, "%s\n", newick.c_str());
            //fclose(subtree_file);
            std::ofstream subtree_file(subtree_filename.c_str(), std::ofstream::out);
            std::stringstream newick_ss;
            new_T.write_newick_string(newick_ss, new_T.root, true, true, retain_original_branch_len);
            subtree_file << newick_ss.rdbuf();
            subtree_file.close();


            // Write list of mutations on the subtree to file
            auto subtree_mutations_filename = outdir + preid + "subtree-" + std::to_string(num_subtrees) + "-mutations.txt";

            fprintf(stderr, "Writing list of mutations at the nodes of subtree %d to file %s\n", num_subtrees, subtree_mutations_filename.c_str());
            FILE* subtree_mutations_file = fopen(subtree_mutations_filename.c_str(), "w");

            for (auto n: new_T.depth_first_expansion()) {
                size_t tot_mutations = n->mutations.size();
                fprintf(subtree_mutations_file, "%s: ", new_T.get_node_name(n->node_id).c_str());
                for (size_t idx = 0; idx < tot_mutations; idx++) {
                    auto m = n->mutations[idx];
                    fprintf(subtree_mutations_file, "%s", m.get_string().c_str());
                    if (idx+1 <tot_mutations) {
                        fprintf(subtree_mutations_file, ",");
                    }
                }
                fprintf(subtree_mutations_file, "\n");
            }
            fclose(subtree_mutations_file);

            // Expand internal nodes that are condensed
            bool has_condensed = false;
            FILE* subtree_expanded_file = NULL;
            for (auto l: new_T.get_leaves()) {
                if (T.condensed_nodes.find(l->node_id) != T.condensed_nodes.end()) {
                    if (!has_condensed) {
                        auto subtree_expanded_filename = outdir + preid +  "subtree-" + std::to_string(num_subtrees) + "-expanded.txt";
                        fprintf(stderr, "Subtree %d has condensed nodes.\nExpanding the condensed nodes for subtree %d in file %s\n", num_subtrees, num_subtrees, subtree_expanded_filename.c_str());
                        subtree_expanded_file = fopen(subtree_expanded_filename.c_str(), "w");
                        has_condensed = true;
                    }
                    fprintf(subtree_expanded_file, "%s: ", T.get_node_name(l->node_id).c_str());
                    auto iter=T.condensed_nodes.find(l->node_id);
                    for (auto n: iter->second) {
                        fprintf(subtree_expanded_file, "%s ", n.c_str());
                    }
                    fprintf(subtree_expanded_file, "\n");
                }
            }
            if (has_condensed) {
                fclose(subtree_expanded_file);
            }
            break;
            anc=anc->parent;
        }
    }
}
