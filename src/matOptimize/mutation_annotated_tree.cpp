#include "mutation_annotated_tree.hpp"
#include "src/matOptimize/check_samples.hpp"
#include <algorithm>
#include <atomic>
#include <csignal>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
// Uses one-hot encoding if base is unambiguous
//void check_leaves(const Mutation_Annotated_Tree::Tree& T);
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
    } else {
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

std::vector<Mutation_Annotated_Tree::Node*> depth_first_expansion_no_tree(Mutation_Annotated_Tree::Node* node)  {
    std::vector<Node*> traversal;
    size_t index=0;
    if (node == NULL) {
        return traversal;
    }
    depth_first_expansion_helper(node, traversal,index,0);
    return traversal;
}
std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::depth_first_expansion(Mutation_Annotated_Tree::Node* node) const {
    TIMEIT();
    if (node == NULL) {
        node = root;
    }
    return depth_first_expansion_no_tree(node);
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
void Mutation_Annotated_Tree::Tree::check_leaves() {
    fprintf(stderr,"===================Node ID so far %zu ===================",node_idx);
    for (const auto node: get_leaves()) {
        auto iter=node_names.find(node->node_id);
        if(iter==node_names.end()&&node->parent) {
            fprintf(stderr,"Node ID %zu name not found",node->node_id);
            raise(SIGTRAP);
        }
    }
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
std::vector<Node*> Mutation_Annotated_Tree::Tree::get_leaves(Node* root) const {
    std::vector<Node*> out;
    if (root==nullptr) {
        root=this->root;
    }
    if(root->is_leaf()){
        out.push_back(root);
        return out;
    }
    get_leaves_helper(root, out);
    return out;
}
Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::Tree::copy_tree() {
    Mutation_Annotated_Tree::Tree out;
    out.node_names=node_names;
    out.node_name_to_idx_map=node_name_to_idx_map;
    out.node_idx=node_idx;
    out.max_level=max_level;
    out.condensed_nodes=condensed_nodes;
    out.curr_internal_node=curr_internal_node;
    out.root=new Mutation_Annotated_Tree::Node(*root,nullptr,&out);
    out.all_nodes.resize(all_nodes.size());
    for (auto node : out.depth_first_expansion()) {
        out.all_nodes[node->node_id]=node;
    }
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

        auto anc_node = sample_node;
        // Mutations on the ancestors of added sample
        while (anc_node) {
            auto curr_node_mutations = anc_node->mutations;
            if (curr_node_mutations.size() > 0) {
                curr_node_mutation_string =
                    T->get_node_name_for_log_output(anc_node->node_id) + ":";
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
                T->get_node_name_for_log_output(sample_node->node_id).c_str());
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
static Node* get_subtree_helper(const std::unordered_set<size_t>& nodes_to_place,const Node* main_tree_node,int num_annotations) {
    std::vector<Node*> child_node;
    for (const Node* child : main_tree_node->children) {
        auto ret=get_subtree_helper(nodes_to_place, child, num_annotations);
        if (ret) {
            child_node.push_back(ret);
        }
    }
    auto iter=nodes_to_place.find(main_tree_node->node_id);
    if (iter!=nodes_to_place.end()||child_node.size()>1) {
        auto new_node=new Node(main_tree_node->node_id);
        new_node->clade_annotations.resize(num_annotations);
        new_node->children=std::move(child_node);
        for (auto node : new_node->children) {
            node->parent=new_node;
        }
        return new_node;
    }
    if (child_node.empty()) {
        return nullptr;
    } else {
        return child_node[0];
    }
};
Mutation_Annotated_Tree::Tree
Mutation_Annotated_Tree::get_subtree(const Mutation_Annotated_Tree::Tree &tree,
                                     const std::vector<Node *> &samples,
                                     bool keep_clade_annotations) {
    //Original_State_t ori_state;
    //check_samples(tree.root, ori_state, &tree);
    TIMEIT();
    Tree subtree=tree.get_data_no_nodes();
    std::vector<Mutation_Annotated_Tree::Node*> new_tree_dfs;
    subtree.root=get_subtree_root(tree, samples, new_tree_dfs,keep_clade_annotations);
    for (auto node : new_tree_dfs) {
        subtree.register_node_serial(node);
    }
    return subtree;
}
Mutation_Annotated_Tree::Node*
Mutation_Annotated_Tree::get_subtree_root(const Mutation_Annotated_Tree::Tree &tree,
                                     const std::vector<Node *> &samples,std::vector<Mutation_Annotated_Tree::Node*>& new_tree_dfs,
                                     bool keep_clade_annotations) {
    //Original_State_t ori_state;
    //check_samples(tree.root, ori_state, &tree);
    TIMEIT();
    size_t num_annotations = 0;
    if (keep_clade_annotations) {
        num_annotations = tree.get_num_annotations();
    }
    std::unordered_set<size_t> nodes_to_place;
    nodes_to_place.reserve(samples.size());
    for (const auto samp : samples) {
        nodes_to_place.insert(samp->node_id);
    }
    Node* root=get_subtree_helper(nodes_to_place, tree.root, num_annotations);
    new_tree_dfs=depth_first_expansion_no_tree(root);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,new_tree_dfs.size()),[&new_tree_dfs,&tree,num_annotations](tbb::blocked_range<size_t> range) {
        for(size_t idx=range.begin(); idx<range.end(); idx++) {
            auto new_tree_node=new_tree_dfs[idx];
            const Node* corresponding_main_tree_node=tree.get_node(new_tree_node->node_id);
            if (corresponding_main_tree_node==tree.root) {
                continue;
            }
            const Node* corresponding_main_parent_node=tree.root;
            if (new_tree_node->parent) {
                corresponding_main_parent_node=tree.get_node(new_tree_node->parent->node_id);
            }
            const Node* curr_node=corresponding_main_tree_node;
            std::unordered_map<int, Mutation> mutations;
            while (curr_node!=corresponding_main_parent_node) {
                for (const auto& mut : curr_node->mutations) {
                    auto res=mutations.emplace(mut.get_position(),mut);
                    if (!res.second) {
                        res.first->second.set_par_one_hot(mut.get_par_one_hot());
                    }
                }
                for (size_t k = 0; k < num_annotations; k++) {
                    if (curr_node->clade_annotations[k] != ""&&new_tree_node->clade_annotations[k]=="") {
                        new_tree_node->clade_annotations[k]=curr_node->clade_annotations[k];
                    }
                }
                curr_node=curr_node->parent;
            }
            new_tree_node->mutations.reserve(mutations.size());
            for (const auto& mut : mutations) {
                if (mut.second.get_par_one_hot()&mut.second.get_mut_one_hot()) {
                    continue;
                }
                new_tree_node->mutations.mutations.push_back(mut.second);
            }
            std::sort(new_tree_node->mutations.begin(),new_tree_node->mutations.end());
        }
    });
    //check_samples(tree.root, ori_state, &tree);
    //check_samples(subtree.root, ori_state, &subtree,true);
    return root;
}
static void rotate_each_node(std::vector<Mutation_Annotated_Tree::Node*>& dfs,bool reverse){
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
void Mutation_Annotated_Tree::Tree::rotate_for_display(bool reverse) {
    auto dfs = depth_first_expansion();
    rotate_each_node(dfs, reverse);    
}
std::vector<Node *> Mutation_Annotated_Tree::Tree::rsearch(Node *node, bool include_self) const{
    std::vector<Node *> output;
    if (node == NULL)
        return output;
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

void Mutation_Annotated_Tree::get_random_single_subtree(
    const Mutation_Annotated_Tree::Tree & T, std::vector<Node *> samples,
    std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx,
    bool retain_original_branch_len, std::vector<Node*> anchor_samples) {
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

    // Add "anchor samples" (if any)
    leaves_to_keep_set.insert(anchor_samples.begin(), anchor_samples.end());

    std::vector<Node *> leaves_to_keep(leaves_to_keep_set.begin(), leaves_to_keep_set.end());

    auto new_T = get_subtree(T, leaves_to_keep);
    //check_leaves(T);
    // Rotate tree for display
    new_T.rotate_for_display();

    // Write subtree to file
    auto subtree_filename = outdir + preid + "single-subtree.nh";
    if (anchor_samples.size() > 0) {
        fprintf(stderr,
                "Writing single subtree with %zu randomly added leaves and %zu anchor samples to file %s.\n",
                subtree_size, anchor_samples.size(), subtree_filename.c_str());
    } else {
        fprintf(stderr,
                "Writing single subtree with %zu randomly added leaves to file %s.\n",
                subtree_size, subtree_filename.c_str());
    }

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
                "%s: ", new_T.get_node_name_for_log_output(n->node_id).c_str());
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
    //check_leaves(T);

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
                    "%s: ", new_T.get_node_name_for_log_output(l->node_id).c_str());
            auto iter=T.condensed_nodes.find(l->node_id);
            for (auto n : iter->second) {
                fprintf(subtree_expanded_file, "%s ", n.c_str());
            }
            fprintf(subtree_expanded_file, "\n");
        }
    }
    //check_leaves(T);
    if (has_condensed) {
        fclose(subtree_expanded_file);
    }
}
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
class Least_N{
    //int cur_max;
    public:
    std::vector<NodeDist> nodes;
    size_t n;
    std::vector<Node*> remaining_nodes;
    void add(NodeDist to_add){
        if(!n){
            remaining_nodes.push_back(to_add.node);
        }
        else if(nodes.size()==n){
            if(to_add.num_mut>=nodes[0].num_mut){
                remaining_nodes.push_back(to_add.node);
            }else {
                std::pop_heap(nodes.begin(),nodes.end());
                remaining_nodes.push_back(nodes.back().node);
                nodes.back()=to_add;
                std::push_heap(nodes.begin(),nodes.end());
                /*if((int)nodes[0].num_mut>cur_max){
                    raise(SIGTRAP);
                }
                cur_max=nodes[0].num_mut;*/
            }
        }else {
            nodes.emplace_back(to_add);
            std::make_heap(nodes.begin(),nodes.end());
            //cur_max=nodes[0].num_mut;
        }
    }
};
static void gather_nodes(Node* start_node,Node* exclude_node,Least_N& out,int dist){
    if(start_node==exclude_node){
        return;
    }
    if(start_node->is_leaf()){
        out.add(NodeDist(start_node,dist));
    }
    for (const auto child : start_node->children) {
        gather_nodes(child, exclude_node, out,dist+child->mutations.size());
    }
}
static int set_num_leaves_as_bfs_idx(Node* root){
    int child_cnt=0;
    if(root->is_leaf()){
        child_cnt=1;
    }
    for (auto child : root->children) {
        child_cnt+=set_num_leaves_as_bfs_idx(child);
    }
    root->bfs_index=child_cnt;
    return child_cnt;
}
void Mutation_Annotated_Tree::get_random_sample_subtrees (const Mutation_Annotated_Tree::Tree& T, std::vector<Node*> samples, std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len, std::vector<Node*> anchor_samples) {
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

    std::unordered_map<Node*, size_t> sample_pos;
    for(size_t i=0; i<samples.size();i++){
        sample_pos.emplace(samples[i],i);
    }
    struct to_extract_info{
        std::vector<Node*> leaves_to_keep;
        std::vector<Node*> to_display_samples;
        void add_leave(Node* node,const std::unordered_map<Node*, size_t>& sample_pos,std::vector<bool>& displayed_samples){
            auto iter=sample_pos.find(node);
            if(iter!=sample_pos.end()){
                displayed_samples[iter->second]=true;
                to_display_samples.push_back(node);
            }
            leaves_to_keep.push_back(node);
        }
    };
    struct to_extract_info_back_inserter: public std::iterator< std::output_iterator_tag,
                                                   void, void, void, void >{
        to_extract_info& content;
        const std::unordered_map<Node*, size_t>& sample_pos;
        std::vector<bool>& displayed_samples;
        to_extract_info_back_inserter(
        to_extract_info& content,
        const std::unordered_map<Node*, size_t>& sample_pos,
        std::vector<bool>& displayed_samples
        ):
        content(content),
        sample_pos(sample_pos),
        displayed_samples(displayed_samples){}
        to_extract_info_back_inserter& operator*(){
            return *this;
        }
        to_extract_info_back_inserter& operator++(){
            return *this;
        }
        to_extract_info_back_inserter& operator++(int){
            return *this;
        }
        to_extract_info_back_inserter& operator=(Node* to_add){
            content.add_leave(to_add, sample_pos, displayed_samples);
            return *this;
        }
    };

    std::vector<to_extract_info> trees_to_extract;
    trees_to_extract.reserve(samples.size());
    // Bool vector to mark which newly placed samples have already been
    // displayed in a subtree (initialized to false)
    std::vector<bool> displayed_samples (samples.size(), false);

    // If the missing sample is not found in the tree, it was not placed
    // because of max_uncertainty. Mark those samples as already
    // displayed.
    set_num_leaves_as_bfs_idx(T.root);
    for (size_t i = 0; i < samples.size(); i++) {

        if (displayed_samples[i]) {
            continue;
        }

        Mutation_Annotated_Tree::Node* last_anc = samples[i];
        trees_to_extract.emplace_back();
        auto& curr_tree=trees_to_extract.back();
        curr_tree.leaves_to_keep.reserve(subtree_size + anchor_samples.size());
        subtree_size=std::min(subtree_size,T.root->bfs_index-1);
        // Keep moving up the tree till a subtree of required size is
        // found
        auto anc=samples[i];
        bool is_first=true;
        while (anc) {
            if (!is_first) {
                anc=anc->parent;
                if (!anc) {
                    break;
                }
            }
            is_first=false;
            size_t num_leaves = anc->bfs_index;
            if (num_leaves < subtree_size) {
                last_anc = anc;
                continue;
            }

            if (num_leaves > subtree_size) {
                for (auto l: T.get_leaves(last_anc)) {
                    curr_tree.add_leave(l,sample_pos,displayed_samples);
                }

                Least_N selected;
                selected.n=nearest_subtree_size>curr_tree.leaves_to_keep.size()?nearest_subtree_size-curr_tree.leaves_to_keep.size():0;
                selected.nodes.reserve(selected.n);
                selected.remaining_nodes.reserve(num_leaves-curr_tree.leaves_to_keep.size());
                gather_nodes(anc, last_anc, selected, 0);
                for (auto n: selected.nodes) {
                    curr_tree.add_leave(n.node,sample_pos,displayed_samples);
                }
                to_extract_info_back_inserter inserter{curr_tree,sample_pos,displayed_samples};
                std::sample(selected.remaining_nodes.begin(), selected.remaining_nodes.end(), inserter, subtree_size-curr_tree.leaves_to_keep.size(), std::default_random_engine{});
            } else {
                for (auto l: T.get_leaves(anc)) {
                     if (curr_tree.leaves_to_keep.size() == subtree_size) {
                        break;
                    }
                    curr_tree.add_leave(l,sample_pos,displayed_samples);
                
                }
            }
            break;
        }
        // Add "anchor samples" (if any)
        curr_tree.leaves_to_keep.insert(curr_tree.leaves_to_keep.end(),
                                        anchor_samples.begin(), anchor_samples.end());
    }
    tbb::parallel_for(tbb::blocked_range<size_t>(0,trees_to_extract.size()),[&](tbb::blocked_range<size_t> range){
        for (int num_subtrees=range.begin(); num_subtrees<(int)range.end(); num_subtrees++) {
            std::vector<Node*> new_tree_dfs;
            auto new_T = Mutation_Annotated_Tree::get_subtree_root(T, trees_to_extract[num_subtrees].leaves_to_keep,new_tree_dfs);
            // Rotate tree for display
            rotate_each_node(new_tree_dfs,false);
            
            std::string node_names;
            for(const auto& node:trees_to_extract[num_subtrees].to_display_samples){
                auto cur_name=T.get_node_name(node->node_id);
                node_names+=(cur_name+",");
            }
            // Write subtree to file
            auto subtree_filename = outdir + preid + "subtree-" + std::to_string(num_subtrees+1) + ".nh";
            fprintf(stderr, "Writing subtree %d containg %s to file %s.\n", num_subtrees+1, node_names.c_str(), subtree_filename.c_str());
            //FILE* subtree_file = fopen(subtree_filename.c_str(), "w");
            //fprintf(subtree_file, "%s\n", newick.c_str());
            //fclose(subtree_file);
            std::ofstream subtree_file(subtree_filename.c_str(), std::ofstream::out);
            std::stringstream newick_ss;
            write_newick_string_node(T,newick_ss, new_T, true, true, retain_original_branch_len);
            subtree_file << newick_ss.rdbuf();
            subtree_file.close();


            // Write list of mutations on the subtree to file
            auto subtree_mutations_filename = outdir + preid + "subtree-" + std::to_string(num_subtrees+1) + "-mutations.txt";

            fprintf(stderr, "Writing list of mutations at the nodes of subtree %d to file %s\n", num_subtrees+1, subtree_mutations_filename.c_str());
            FILE* subtree_mutations_file = fopen(subtree_mutations_filename.c_str(), "w");

            for (auto n: new_tree_dfs) {
                size_t tot_mutations = n->mutations.size();
                fprintf(subtree_mutations_file, "%s: ", T.get_node_name_for_log_output(n->node_id).c_str());
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
            if(!T.condensed_nodes.empty()){
                fprintf(stderr, "checking for condensed nodes\n");
            for (auto l: new_tree_dfs) {
                if(l->children.size()){
                    continue;
                }
                if (T.condensed_nodes.find(l->node_id) != T.condensed_nodes.end()) {
                    if (!has_condensed) {
                        auto subtree_expanded_filename = outdir + preid +  "subtree-" + std::to_string(num_subtrees) + "-expanded.txt";
                        fprintf(stderr, "Subtree %d has condensed nodes.\nExpanding the condensed nodes for subtree %d in file %s\n", num_subtrees, num_subtrees, subtree_expanded_filename.c_str());
                        subtree_expanded_file = fopen(subtree_expanded_filename.c_str(), "w");
                        has_condensed = true;
                    }
                    fprintf(subtree_expanded_file, "%s: ", T.get_node_name_for_log_output(l->node_id).c_str());
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
            }
            break;
        }
    });

}