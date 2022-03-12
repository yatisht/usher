#include "mutation_annotated_tree.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <algorithm>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include "usher_graph.hpp"
#include <signal.h>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8
int8_t Mutation_Annotated_Tree::get_nuc_id (char nuc) {
    int8_t ret = 0b1111;
    switch(nuc) {
    case 'a':
    case 'A':
        ret = 0b1;
        break;
    case 'c':
    case 'C':
        ret = 0b10;
        break;
    case 'g':
    case 'G':
        ret = 0b100;
        break;
    case 't':
    case 'T':
        ret = 0b1000;
        break;
    case 'R':
        ret = 0b101;
        break;
    case 'Y':
        ret = 0b1010;
        break;
    case 'S':
        ret = 0b110;
        break;
    case 'W':
        ret = 0b1001;
        break;
    case 'K':
        ret = 0b1100;
        break;
    case 'M':
        ret = 0b11;
        break;
    case 'B':
        ret = 0b1110;
        break;
    case 'D':
        ret = 0b1101;
        break;
    case 'H':
        ret = 0b1011;
        break;
    case 'V':
        ret = 0b111;
    case 'n':
    case 'N':
    default:
        ret = 0b1111;
        break;
    }
    return ret;
}

// Sets bits at positions specified by nuc_vec to 1 in int8
int8_t Mutation_Annotated_Tree::get_nuc_id (std::vector<int8_t> nuc_vec) {
    int8_t ret = 0;
    int8_t one = 1;
    for (auto nuc: nuc_vec) {
        assert((nuc >= 0) && (nuc <=3));
        ret += (one << nuc);
    }
    return ret;
}

// Convert nuc_id back to IUPAC base
char Mutation_Annotated_Tree::get_nuc (int8_t nuc_id) {
    char ret = 'N';
    //assert ((nuc_id >= 1) && (nuc_id <= 15));
    switch(nuc_id) {
    case 1:
        ret = 'A';
        break;
    case 2:
        ret = 'C';
        break;
    case 3:
        ret = 'M';
        break;
    case 4:
        ret = 'G';
        break;
    case 5:
        ret = 'R';
        break;
    case 6:
        ret = 'S';
        break;
    case 7:
        ret = 'V';
        break;
    case 8:
        ret = 'T';
        break;
    case 9:
        ret = 'W';
        break;
    case 10:
        ret = 'Y';
        break;
    case 11:
        ret = 'H';
        break;
    case 12:
        ret = 'K';
        break;
    case 13:
        ret = 'D';
        break;
    case 14:
        ret = 'B';
        break;
    default:
        ret = 'N';
        break;
    }
    return ret;
}

// A:0, C:1, G:2, T:3
int8_t Mutation_Annotated_Tree::get_nt (int8_t nuc_id) {
    int8_t ret = 0;
    switch(nuc_id) {
    case 1:
        ret = 0;
        break;
    case 2:
        ret = 1;
        break;
    case 4:
        ret = 2;
        break;
    case 8:
        ret = 3;
        break;
    default:
        ret = -1;
        break;
    }
    return ret;
}

std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec (char c) {
    switch (c) {
    case 'a':
    case 'A':
        return std::vector<int8_t> {0};
    case 'c':
    case 'C':
        return std::vector<int8_t> {1};
    case 'g':
    case 'G':
        return std::vector<int8_t> {2};
    case 't':
    case 'T':
        return std::vector<int8_t> {3};
    case 'R':
        return std::vector<int8_t> {0,2};
    case 'Y':
        return std::vector<int8_t> {1,3};
    case 'S':
        return std::vector<int8_t> {1,2};
    case 'W':
        return std::vector<int8_t> {0,3};
    case 'K':
        return std::vector<int8_t> {2,3};
    case 'M':
        return std::vector<int8_t> {0,1};
    case 'B':
        return std::vector<int8_t> {1,2,3};
    case 'D':
        return std::vector<int8_t> {0,2,3};
    case 'H':
        return std::vector<int8_t> {0,1,3};
    case 'V':
        return std::vector<int8_t> {0,1,2};
    case 'n':
    case 'N':
        return std::vector<int8_t> {0,1,2,3};
    default:
        return std::vector<int8_t> {0,1,2,3};
    }
}

std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec_from_id (int8_t nuc_id) {
    return get_nuc_vec(get_nuc(nuc_id));
}

// Get newick stringstream for the input subtree rooted at some node (node) in
// the input tree T. Boolean arguments decide whether
// internal node ids and branch lengths are printed. If last boolean argument is
// set, branch lengths from input tree are retained, otherwise, branch length
// for a branch is equal to the number of mutations annotated on that branch
void Mutation_Annotated_Tree::write_newick_string (std::stringstream& ss, const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    TIMEIT();

    std::vector<Node*> traversal = T.depth_first_expansion(node);
    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->branch_length;
        if (!retain_original_branch_len) {
            branch_length = static_cast<float>(n->mutations.size());
        }
        if (curr_level < level) {
            if (!prev_open) {
                ss << ',';
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                ss << '(';
                prev_open = true;
            }
            if (n->is_leaf()) {
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                } else {
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
                }
                prev_open = false;
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        } else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                ss << ')';
                if (print_internal) {
                    ss << node_stack.top();
                }
                if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
                    ss << ':';
                    ss << branch_length_stack.top();
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->is_leaf()) {
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    ss << ',';
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                } else {
                    ss << ',';
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
                }
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        } else {
            prev_open = false;
            if (n->is_leaf()) {
                if (uncondense_leaves && (T.condensed_nodes.find(n->identifier) != T.condensed_nodes.end())) {
                    auto cn = T.condensed_nodes.at(n->identifier);
                    ss << ',';
                    auto cn_size = cn.size();
                    for (size_t idx = 0; idx < cn_size; idx++) {
                        ss << cn[idx];
                        if (idx+1 < cn_size) {
                            ss << ',';
                        }
                    }
                } else {
                    ss << ',';
                    ss << n->identifier;
                }
                if ((print_branch_len) && (branch_length >= 0)) {
                    ss << ':';
                    ss << branch_length;
                }
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        }
        curr_level = level;
    }
    size_t remaining = node_stack.size();
    for (size_t i = 0; i < remaining; i++) {
        ss << ')';
        if (print_internal) {
            ss << node_stack.top();
        }
        if ((print_branch_len) && (branch_length_stack.top() >= 0)) {
            ss << ':';
            ss << branch_length_stack.top();
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    ss << ';';
}

std::string Mutation_Annotated_Tree::get_newick_string (const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    std::stringstream newick_ss;
    write_newick_string(newick_ss, T, node, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
    return newick_ss.str();
}

std::string Mutation_Annotated_Tree::get_newick_string (const Tree& T, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    return get_newick_string(T, T.root, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
}

Mutation_Annotated_Tree::Mutation* Mutation_Annotated_Tree::mutation_from_string(const std::string& mut_string) {
    // Parse string like /[A-Z][0-9]+[A-Z]/ e.g. A23403G; create & return (pointer to) Mutation_Annotated_Tree::Mutation
    // with par_nuc set to ref_nuc (same as VCF parsing code in usher.cpp).
    // If mut_string is not in expected format, print error message and return NULL.
    char ref, alt;
    int position;
    int ret = sscanf(mut_string.c_str(), "%c%d%c", &ref, &position, &alt);
    if (ret != 3 || ref < 'A' || ref > 'Z' || alt < 'A' || alt > 'Z') {
        fprintf(stderr, "mutation_from_string: expected /[A-Z][0-9]+[A-Z/, got '%s'\n", mut_string.c_str());
        return NULL;
    }
    Mutation_Annotated_Tree::Mutation* mut = new Mutation_Annotated_Tree::Mutation();
    mut->ref_nuc = Mutation_Annotated_Tree::get_nuc_id(ref);
    mut->position = position;
    mut->mut_nuc = Mutation_Annotated_Tree::get_nuc_id(alt);
    mut->par_nuc = mut->ref_nuc;
    // Double-check to make sure there aren't additional characters past /[A-Z][0-9]+[A-Z]/:
    if (mut->get_string() != mut_string) {
        fprintf(stderr, "mutation_from_string: unexpected characters at the end of '%s'\n", mut_string.c_str());
        mut = NULL;
    }
    return mut;
}

// Split string into words for a specific delimiter delim
void Mutation_Annotated_Tree::string_split (std::string const& s, char delim, std::vector<std::string>& words) {
    TIMEIT();
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        // if ((end_pos == start_pos) || end_pos >= s.length()) {
        if (end_pos >= s.length()) {
            break;
        }
        words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}

// Split string into words (delimited by space, tabs etc.)
void Mutation_Annotated_Tree::string_split (std::string s, std::vector<std::string>& words) {
    std::string curr = "";
    std::vector<std::string> ret;

    // Used to split string around spaces.
    std::istringstream ss(s);

    std::string word;
    // Traverse through all words
    while (ss >> word) {
        words.push_back(std::move(word));
    };
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::create_tree_from_newick_string (std::string newick_string) {
    TIMEIT();
    Tree T;

    std::vector<std::string> leaves;
    std::vector<size_t> num_open;
    std::vector<size_t> num_close;
    std::vector<std::queue<float>> branch_len (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    string_split(newick_string, ',', s1);

    num_open.reserve(s1.size());
    num_close.reserve(s1.size());

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        bool stop = false;
        bool branch_start = false;
        std::string leaf = "";
        std::string branch = "";
        for (auto c: s) {
            if (c == ':') {
                stop = true;
                branch = "";
                branch_start = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branch_len.size() <= level) {
                    branch_len.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                branch_len[level].push(len);
                level--;
                branch_start = false;
            } else if (!stop) {
                leaf += c;
                branch_start = false;
            } else if (branch_start) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        leaves.push_back(std::move(leaf));
        num_open.push_back(no);
        num_close.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        branch_len[level].push(len);
    }

    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    std::stack<Node*> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = T.new_internal_node_id();
            Node* new_node = NULL;
            if (parent_stack.size() == 0) {
                new_node = T.create_node(nid, branch_len[level].front());
            } else {
                new_node = T.create_node(nid, parent_stack.top(), branch_len[level].front());
            }
            branch_len[level].pop();
            level++;
            parent_stack.push(new_node);
        }
        T.create_node(leaf, parent_stack.top(), branch_len[level].front());
        branch_len[level].pop();
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
            level--;
        }
    }

    if (T.root == NULL) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    return T;
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::create_tree_from_newick (std::string filename) {
    std::ifstream infile(filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the tree file: %s!\n", filename.c_str());
        exit(1);
    }
    std::string newick_string;
    std::getline(infile, newick_string);

    return create_tree_from_newick_string(newick_string);
}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::load_mutation_annotated_tree (std::string filename) {
    TIMEIT();
    Tree tree;

    Parsimony::data data;
#define BIG_SIZE 2000000000l
    boost::iostreams::filtering_istream instream;
    std::ifstream inpfile(filename, std::ios::in | std::ios::binary);
    if (filename.find(".gz\0") != std::string::npos) {
        if (!inpfile) {
            fprintf(stderr, "ERROR: Could not load the mutation-annotated tree object from file: %s!\n", filename.c_str());
            exit(1);
        }
        try {
            instream.push(boost::iostreams::gzip_decompressor());
            instream.push(inpfile);
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        instream.push(inpfile);
    }
    google::protobuf::io::IstreamInputStream stream(&instream);
    google::protobuf::io::CodedInputStream input(&stream);
    //input.SetTotalBytesLimit(BIG_SIZE, BIG_SIZE);
    data.ParseFromCodedStream(&input);
    //check if the pb has a metadata field
    bool hasmeta = (data.metadata_size()>0);
    if (!hasmeta) {
        fprintf(stderr, "WARNING: This pb does not include any metadata. Filling in default values\n");
    }
    tree = create_tree_from_newick_string(data.newick());
    auto dfs = tree.depth_first_expansion();
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, dfs.size()),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t idx = r.begin(); idx < r.end(); idx++) {
            auto node = dfs[idx];
            auto mutation_list = data.node_mutations(idx);
            if (hasmeta) {
                for (int k = 0; k < data.metadata(idx).clade_annotations_size(); k++) {
                    node->clade_annotations.emplace_back(data.metadata(idx).clade_annotations(k));
                }
            }
            for (int k = 0; k < mutation_list.mutation_size(); k++) {
                auto mut = mutation_list.mutation(k);
                Mutation m;
                m.chrom = mut.chromosome();
                m.position = mut.position();
                if (!m.is_masked()) {
                    m.ref_nuc = (1 << mut.ref_nuc());
                    m.par_nuc = (1 << mut.par_nuc());
                    m.is_missing = false;
                    std::vector<int8_t> nuc_vec(mut.mut_nuc_size());
                    for (int n = 0; n < mut.mut_nuc_size(); n++) {
                        nuc_vec[n] = mut.mut_nuc(n);
                    }
                    m.mut_nuc = get_nuc_id(nuc_vec);
                    if (m.mut_nuc != m.par_nuc) {
                        node->add_mutation(m);
                    }
                } else {
                    // Mutation masked
                    m.ref_nuc = 0;
                    m.par_nuc = 0;
                    m.mut_nuc = 0;
                    node->add_mutation(m);
                }
            }
            if (!std::is_sorted(node->mutations.begin(), node->mutations.end())) {
                fprintf(stderr, "WARNING: Mutations not sorted!\n");
                std::sort(node->mutations.begin(), node->mutations.end());
            }
        }
    }, ap);

    size_t num_condensed_nodes = static_cast<size_t>(data.condensed_nodes_size());
    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_condensed_nodes),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t idx = r.begin(); idx < r.end(); idx++) {
            auto cn = data.condensed_nodes(idx);
            tree.condensed_nodes.emplace(std::pair<std::string, std::vector<std::string>>(cn.node_name(), std::vector<std::string>(cn.condensed_leaves_size())));
            for (int k = 0; k < cn.condensed_leaves_size(); k++) {
                tree.condensed_nodes[cn.node_name()][k] = cn.condensed_leaves(k);
                tree.condensed_leaves.emplace(cn.condensed_leaves(k));
            }
        }
    }, ap);

    return tree;
}

void Mutation_Annotated_Tree::save_mutation_annotated_tree (Mutation_Annotated_Tree::Tree tree, std::string filename) {
    TIMEIT();
    Parsimony::data data;
    data.set_newick(get_newick_string(tree, false, true, true));

    auto dfs = tree.depth_first_expansion();

    for (size_t idx = 0; idx < dfs.size(); idx++) {
        auto meta = data.add_metadata();
        for (size_t k = 0; k < dfs[idx]->clade_annotations.size(); k++) {
            meta->add_clade_annotations(dfs[idx]->clade_annotations[k]);
        }
        auto mutation_list = data.add_node_mutations();
        for (auto m: dfs[idx]->mutations) {
            auto mut = mutation_list->add_mutation();
            mut->set_chromosome(m.chrom);
            mut->set_position(m.position);

            if (m.is_masked()) {
                mut->set_ref_nuc(-1);
                mut->set_par_nuc(-1);
            } else {
                int8_t j = get_nt(m.ref_nuc);
                assert (j >= 0);
                mut->set_ref_nuc(j);

                j = get_nt(m.par_nuc);
                assert(j >= 0);
                mut->set_par_nuc(j);

                mut->clear_mut_nuc();
                for (auto nuc: get_nuc_vec_from_id(m.mut_nuc)) {
                    mut->add_mut_nuc(nuc);
                }
            }
        }
    }

    // Add condensed nodes
    for (auto cn: tree.condensed_nodes) {
        auto cn_ptr = data.add_condensed_nodes();
        cn_ptr->set_node_name(cn.first);
        for (auto lid: cn.second) {
            cn_ptr->add_condensed_leaves(lid);
        }
    }

    // Boost library used to stream the contents to the output protobuf file in
    // uncompressed or compressed .gz format
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outbuf;

    if (filename.find(".gz\0") != std::string::npos) {
        try {
            outbuf.push(boost::iostreams::gzip_compressor());
            outbuf.push(outfile);
            std::ostream outstream(&outbuf);
            data.SerializeToOstream(&outstream);
            boost::iostreams::close(outbuf);
            outfile.close();
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        data.SerializeToOstream(&outfile);
        outfile.close();
    }
}

/* === Node === */
bool Mutation_Annotated_Tree::Node::is_leaf () {
    return (children.size() == 0);
}

bool Mutation_Annotated_Tree::Node::is_root() {
    return (parent == NULL);
}

Mutation_Annotated_Tree::Node::Node() {
    level = 0;
    identifier = "";
    parent = NULL;
    branch_length = -1.0;
    clade_annotations.clear();
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, float len) {
    identifier = id;
    parent = NULL;
    level = 1;
    branch_length = len;
    mutations.clear();
}

Mutation_Annotated_Tree::Node::Node (std::string id, Node* p, float len) {
    identifier = id;
    parent = p;
    level = p->level + 1;
    branch_length = len;
    mutations.clear();
}

// Assumes mutations are added in chronological order. If a new mutation occurs
// at the same position, it should either be updated to the new allele or
// removed entirely (in case of reversal mutation)
void Mutation_Annotated_Tree::Node::add_mutation (Mutation mut) {
    auto iter = std::lower_bound(mutations.begin(), mutations.end(), mut);
    // check if mutation at the same position has occured before
    if ((iter != mutations.end()) && (iter->position == mut.position)) {
        // update to new allele
        if (iter->par_nuc != mut.mut_nuc) {
            iter->mut_nuc = mut.mut_nuc;
        }
        //reversal mutation
        else {
            std::vector<Mutation> tmp;
            for (auto m: mutations) {
                if (m.position != iter->position) {
                    tmp.emplace_back(m.copy());
                }
            }
            mutations.clear();
            for (auto m: tmp) {
                mutations.emplace_back(m.copy());
            }
        }
    }
    // new mutation
    else {
        mutations.insert(iter, mut);
    }
}

void Mutation_Annotated_Tree::Node::clear_mutations() {
    mutations.clear();
}

void Mutation_Annotated_Tree::Node::clear_annotations() {
    clade_annotations.clear();
}

/* === Tree === */
size_t Mutation_Annotated_Tree::Tree::get_max_level () const {
    size_t max_level = 0;
    for (auto x: all_nodes) {
        if (x.second->level > max_level) {
            max_level = x.second->level;
        }
    }
    return max_level;
}

size_t Mutation_Annotated_Tree::Tree::get_num_annotations () const {
    size_t ret = 0;
    if (root != NULL) {
        ret = root->clade_annotations.size();
    }
    return ret;
}

void Mutation_Annotated_Tree::Tree::rename_node(std::string old_nid, std::string new_nid) {
    auto n = get_node(old_nid);
    if (n != NULL) {
        if (all_nodes.find(new_nid) != all_nodes.end()) {
            fprintf(stderr, "ERROR: rename_node: node with id '%s' already exists.\n", new_nid.c_str());
            exit(1);
        }
        n->identifier = new_nid;
        all_nodes.erase(old_nid);
        all_nodes[new_nid] = n;
    } else {
        fprintf(stderr, "ERROR: %s not found in the Tree!\n", old_nid.c_str());
        exit(1);
    }
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::get_leaves(std::string nid) {
    std::vector<Node*> leaves;
    if (nid == "") {
        if (root == NULL) {
            return leaves;
        }
        nid = root->identifier;
    }
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves.push_back(curr_node);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves;
}

std::vector<std::string> Mutation_Annotated_Tree::Tree::get_leaves_ids(std::string nid) {
    std::vector<std::string> leaves_ids;
    if (nid == "") {
        if (root == NULL) {
            return leaves_ids;
        }
        nid = root->identifier;
    }
    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        if (curr_node->children.size() == 0)
            leaves_ids.push_back(curr_node->identifier);
        remaining_nodes.pop();
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves_ids;
}

size_t Mutation_Annotated_Tree::Tree::get_num_leaves(Node* node) {
    if (node == NULL) {
        node = root;
    }

    if (node->is_leaf()) {
        return 1;
    }
    size_t num_leaves = 0;
    for (auto c: node->children) {
        num_leaves += get_num_leaves(c);
    }
    return num_leaves;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, float branch_len, size_t num_annotations) {
    all_nodes.clear();
    Node* n = new Node(identifier, branch_len);
    for (size_t k=0; k < num_annotations; k++) {
        n->clade_annotations.emplace_back("");
    }
    root = n;
    all_nodes[identifier] = root;
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, Node* par, float branch_len) {
    if (all_nodes.find(identifier) != all_nodes.end()) {
        fprintf(stderr, "Error: %s already in the tree!\n", identifier.c_str());
        exit(1);
    }
    Node* n = new Node(identifier, par, branch_len);
    size_t num_annotations = get_num_annotations();
    for (size_t k=0; k < num_annotations; k++) {
        n->clade_annotations.emplace_back("");
    }
    all_nodes[identifier] = n;
    par->children.push_back(n);
    return n;
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::create_node (std::string const& identifier, std::string const& parent_id, float branch_len) {
    Node* par = all_nodes[parent_id];
    return create_node(identifier, par, branch_len);
}

Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::Tree::get_node (std::string nid) const {
    if (all_nodes.find(nid) != all_nodes.end()) {
        return all_nodes.at(nid);
    }
    return NULL;

}

bool Mutation_Annotated_Tree::Tree::is_ancestor (std::string anc_id, std::string nid) const {
    Node* node = get_node(nid);
    while (node->parent != NULL) {
        node = node->parent;
        if (node->identifier == anc_id) {
            return true;
        }
    }
    return false;
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::rsearch (const std::string& nid, bool include_self) const {
    std::vector<Node*> ancestors;
    Node* node = get_node(nid);
    if (node==NULL) {
        return ancestors;
    }
    if (include_self) {
        ancestors.reserve(node->level+1);
        ancestors.emplace_back(node);
    } else {
        ancestors.reserve(node->level);
    }
    while (node->parent != NULL) {
        ancestors.emplace_back(node->parent);
        node = node->parent;
    }
    return ancestors;
}

std::string Mutation_Annotated_Tree::Tree::get_clade_assignment (const Node* n, int clade_id, bool include_self) const {
    assert ((size_t)clade_id < get_num_annotations());
    for (auto anc: rsearch(n->identifier, include_self)) {
        if ((int)anc->clade_annotations.size() > clade_id && anc->clade_annotations[clade_id] != "") {
            return anc->clade_annotations[clade_id];
        }
    }
    return "UNDEFINED";
}

void Mutation_Annotated_Tree::Tree::remove_node_helper (std::string nid, bool move_level) {
    auto it = all_nodes.find(nid);
    if (it == all_nodes.end()) {
        fprintf(stderr, "ERROR: Tried to remove node identifier %s but it was not found!\n", nid.c_str());
        exit(1);
    }
    Node* source = it->second;
    Node* curr_parent = source->parent;

    if (curr_parent != NULL) {
        // Remove source from curr_parent
        auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), source);
        assert (iter != curr_parent->children.end());
        curr_parent->children.erase(iter);

        // Remove parent if it no longer has any children
        if (curr_parent->children.size() == 0) {
            if (curr_parent == root) {
                fprintf(stderr, "ERROR: Tree empty!\n");
                exit(1);
            }
            remove_node_helper (curr_parent->identifier, move_level);
        }
        // Move the remaining child one level up if it is the only child of its parent
        else if (move_level && (curr_parent->children.size() == 1)) {
            auto child = curr_parent->children[0];
            if (curr_parent->parent != NULL) {
                for (size_t k=0; k < curr_parent->clade_annotations.size(); k++) {
                    if (child->clade_annotations[k] == "") {
                        child->clade_annotations[k] = curr_parent->clade_annotations[k];
                    }
                }
                child->parent = curr_parent->parent;
                child->level = curr_parent->parent->level + 1;
                child->branch_length += curr_parent->branch_length;

                std::vector<Mutation> tmp;
                for (auto m: child->mutations) {
                    tmp.emplace_back(m);
                }

                //Clear and add back mutations in chrono order
                child->clear_mutations();
                for (auto m: curr_parent->mutations) {
                    child->add_mutation(m);
                }
                for (auto m: tmp) {
                    child->add_mutation(m);
                }

                curr_parent->parent->children.push_back(child);

                iter = std::find(curr_parent->parent->children.begin(), curr_parent->parent->children.end(), curr_parent);
                assert(iter != curr_parent->parent->children.end());
                curr_parent->parent->children.erase(iter);

                // Update levels of source descendants
                std::queue<Node*> remaining_nodes;
                remaining_nodes.push(child);
                while (remaining_nodes.size() > 0) {
                    Node* curr_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    curr_node->level = curr_node->parent->level + 1;
                    for (auto c: curr_node->children) {
                        remaining_nodes.push(c);
                    }
                }
            }

            auto par_it = all_nodes.find(curr_parent->identifier);
            assert (par_it != all_nodes.end());
            all_nodes.erase(par_it);
            delete curr_parent;
        }
    }

    //Remove source and descendants from all_nodes
    std::queue<Node*> desc;
    desc.push(source);
    while (desc.size() > 0) {
        Node* curr_node = desc.front();
        desc.pop();
        for (auto c: curr_node->children) {
            desc.push(c);
        }
        it = all_nodes.find(curr_node->identifier);
        all_nodes.erase(it);
        delete curr_node;
    }
}

void Mutation_Annotated_Tree::Tree::remove_node (std::string nid, bool move_level) {
    TIMEIT();
    remove_node_helper (nid, move_level);
}

void Mutation_Annotated_Tree::Tree::remove_single_child_nodes() {
    auto bfs = breadth_first_expansion();
    for (auto n: bfs) {
        if ((n == root) || (n->children.size() != 1)) {
            continue;
        }

        auto curr_parent = n;
        auto child = n->children[0];
        if (curr_parent->parent != NULL) {
            child->parent = curr_parent->parent;
            child->level = curr_parent->parent->level + 1;
            child->branch_length += curr_parent->branch_length;

            std::vector<Mutation> tmp;
            for (auto m: child->mutations) {
                tmp.emplace_back(m);
            }

            //Clear and add back mutations in chrono order
            child->clear_mutations();
            for (auto m: curr_parent->mutations) {
                child->add_mutation(m);
            }
            for (auto m: tmp) {
                child->add_mutation(m);
            }

            curr_parent->parent->children.push_back(child);

            auto iter = std::find(curr_parent->parent->children.begin(), curr_parent->parent->children.end(), curr_parent);
            assert(iter != curr_parent->parent->children.end());
            curr_parent->parent->children.erase(iter);

            // Update levels of source descendants
            std::queue<Node*> remaining_nodes;
            remaining_nodes.push(child);
            while (remaining_nodes.size() > 0) {
                Node* curr_node = remaining_nodes.front();
                remaining_nodes.pop();
                curr_node->level = curr_node->parent->level + 1;
                for (auto c: curr_node->children) {
                    remaining_nodes.push(c);
                }
            }

            auto par_it = all_nodes.find(curr_parent->identifier);
            assert (par_it != all_nodes.end());
            all_nodes.erase(par_it);
            auto to_delete = curr_parent;
            curr_parent = curr_parent->parent;
            delete to_delete;
        }
    }
}

void Mutation_Annotated_Tree::Tree::move_node (std::string source_id, std::string dest_id, bool move_level) {
    Node* source = all_nodes[source_id];
    Node* destination = all_nodes[dest_id];
    Node* curr_parent = source->parent;

    source->parent = destination;
    source->branch_length = -1.0; // Invalidate source branch length

    destination->children.push_back(source);

    // Remove source from curr_parent
    auto iter = std::find(curr_parent->children.begin(), curr_parent->children.end(), source);
    curr_parent->children.erase(iter);
    if (curr_parent->children.size() == 0) {
        remove_node(curr_parent->identifier, move_level);
    }

    // Update levels of source descendants
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(source);
    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        remaining_nodes.pop();
        curr_node->level = curr_node->parent->level + 1;
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::breadth_first_expansion(std::string nid) {
    std::vector<Node*> traversal;

    if (nid == "") {
        if (root == NULL) {
            return traversal;
        }
        nid = root->identifier;
    }

    Node* node = all_nodes[nid];

    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);

    while (remaining_nodes.size() > 0) {
        Node* curr_node = remaining_nodes.front();
        traversal.push_back(curr_node);
        remaining_nodes.pop();

        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }

    return traversal;
}

void Mutation_Annotated_Tree::Tree::depth_first_expansion_helper(Mutation_Annotated_Tree::Node* node, std::vector<Mutation_Annotated_Tree::Node*>& vec) const {
    node->dfs_idx=vec.size();
    vec.push_back(node);
    for (auto c: node->children) {
        depth_first_expansion_helper(c, vec);
    }
    node->dfs_end_idx=vec.size();
}

std::vector<Mutation_Annotated_Tree::Node*> Mutation_Annotated_Tree::Tree::depth_first_expansion(Mutation_Annotated_Tree::Node* node) const {
    TIMEIT();
    std::vector<Node*> traversal;
    if (node == NULL) {
        node = root;
    }
    if (node == NULL) {
        return traversal;
    }
    depth_first_expansion_helper(node, traversal);
    return traversal;
}

size_t Mutation_Annotated_Tree::Tree::get_parsimony_score() {
    size_t score = 0;
    auto dfs = depth_first_expansion();
    for (auto n: dfs) {
        score += n->mutations.size();
        /*if ((!n->mutations.empty())&&n->mutations[0].is_masked()) {
            score--;
        }*/
    }
    return score;
}

void Mutation_Annotated_Tree::Tree::condense_leaves(std::vector<std::string> missing_samples) {
    if (condensed_nodes.size() > 0) {
        fprintf(stderr, "WARNING: tree contains condensed nodes. Uncondensing fist.\n");
        uncondense_leaves();
    }

    auto tree_leaves = get_leaves_ids();
    for (auto l1_id: tree_leaves) {
        std::vector<Node*> polytomy_nodes;

        auto l1 = get_node(l1_id);
        if (l1 == NULL) {
            continue;
        }
        if (std::find(missing_samples.begin(), missing_samples.end(), l1->identifier) != missing_samples.end()) {
            continue;
        }
        if (l1->mutations.size() > 0) {
            continue;
        }

        for (auto l2: l1->parent->children) {
            if (std::find(missing_samples.begin(), missing_samples.end(), l2->identifier) != missing_samples.end()) {
                continue;
            }
            if (l2->is_leaf() && (get_node(l2->identifier) != NULL) && (l2->mutations.size() == 0)) {
                polytomy_nodes.push_back(l2);
            }
        }
        if (polytomy_nodes.size() > 1) {
            std::string new_node_name = "node_" + std::to_string(1+condensed_nodes.size()) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";

            auto curr_node = get_node(l1->identifier);
            auto new_node = create_node(new_node_name, curr_node->parent, l1->branch_length);

            new_node->clear_mutations();

            condensed_nodes[new_node_name] = std::vector<std::string>(polytomy_nodes.size());

            for (size_t it = 0; it < polytomy_nodes.size(); it++) {
                condensed_nodes[new_node_name][it] = polytomy_nodes[it]->identifier;
                remove_node(polytomy_nodes[it]->identifier, false);
            }
        }
    }
}

void Mutation_Annotated_Tree::Tree::uncondense_leaves() {
    for (auto cn = condensed_nodes.begin(); cn!=condensed_nodes.end(); cn++) {

        auto n = get_node(cn->first);
        auto par = (n->parent != NULL) ? n->parent : n;

        size_t num_samples = cn->second.size();

        if ((num_samples > 1) && (n->mutations.size() > 0)) {
            all_nodes.erase(n->identifier);

            n->identifier = new_internal_node_id();
            all_nodes[n->identifier] = n;

            for (size_t s = 0; s < num_samples; s++) {
                Node* new_n = new Node(cn->second[s], n, -1);
                size_t num_annotations = get_num_annotations();
                for (size_t k=0; k < num_annotations; k++) {
                    new_n->clade_annotations.emplace_back("");
                }
                all_nodes[cn->second[s]] = new_n;

                n->children.push_back(new_n);
            }
        } else if (num_samples > 1) {
            all_nodes.erase(n->identifier);

            n->identifier = cn->second[0];
            all_nodes[n->identifier] = n;

            for (size_t s = 1; s < num_samples; s++) {
                Node* new_n = new Node(cn->second[s], par, n->branch_length);
                size_t num_annotations = get_num_annotations();
                for (size_t k=0; k < num_annotations; k++) {
                    new_n->clade_annotations.emplace_back("");
                }
                all_nodes[cn->second[s]] = new_n;
                par->children.push_back(new_n);
            }
        } else if (num_samples == 1) {
            all_nodes.erase(n->identifier);

            n->identifier = cn->second[0];
            all_nodes[n->identifier] = n;
        }
    }
    condensed_nodes.clear();
    condensed_leaves.clear();
}

void Mutation_Annotated_Tree::Tree::collapse_tree() {
    auto bfs = breadth_first_expansion();

    for (size_t idx = 1; idx < bfs.size(); idx++) {
        auto node = bfs[idx];
        auto mutations = node->mutations;
        if (mutations.size() == 0) {
            auto parent = node->parent;
            auto children = node->children;
            for (auto child: children) {
                move_node(child->identifier, parent->identifier, false);
            }
        }
        //If internal node has one child, the child can be moved up one level
        else if (node->children.size() == 1) {
            auto child = node->children.front();
            auto parent = node->parent;
            for (auto m: mutations) {
                child->add_mutation(m.copy());
            }
            move_node(child->identifier, parent->identifier, false);
        }
    }
}

void Mutation_Annotated_Tree::Tree::rotate_for_display(bool reverse) {
    auto dfs = depth_first_expansion();

    std::unordered_map<Node*, int> num_desc;

    for (int i=int(dfs.size())-1; i>=0; i--) {
        auto n = dfs[i];
        int desc = 1;
        for (auto child: n->children) {
            desc += num_desc[child];
        }
        num_desc[n] = desc;
    }

    for (auto n: dfs) {
        if (reverse) {
            tbb::parallel_sort(n->children.begin(), n->children.end(),
            [&num_desc](Node* n1, Node* n2) {
                return num_desc[n1] < num_desc[n2];
            });
        } else {
            tbb::parallel_sort(n->children.begin(), n->children.end(),
            [&num_desc](Node* n1, Node* n2) {
                return num_desc[n1] > num_desc[n2];
            });
        }
    }
}

void Mutation_Annotated_Tree::Tree::rotate_for_consistency() {
    auto dfs = depth_first_expansion();

    std::unordered_map<Node*, int> num_desc;
    std::unordered_map<Node*, std::string> smallest_desc;

    for (int i=int(dfs.size())-1; i>=0; i--) {
        auto n = dfs[i];
        int desc = 1;
        for (auto child: n->children) {
            desc += num_desc[child];
        }
        num_desc[n] = desc;

        if (n->is_leaf()) {
            smallest_desc[n] = n->identifier;
        } else {
            std::string smallest = smallest_desc[n->children[0]];
            for (auto child: n->children) {
                if (smallest_desc[child] < smallest) {
                    smallest = smallest_desc[child];
                }
            }
            smallest_desc[n] = smallest;
        }
    }

    for (auto n: dfs) {
        tbb::parallel_sort(n->children.begin(), n->children.end(),
        [&](Node* n1, Node* n2) {
            return ((smallest_desc[n1] < smallest_desc[n2]) ||
                    ((smallest_desc[n1] == smallest_desc[n2]) &&
                     (num_desc[n1] < num_desc[n2])));
        });
    }
}


Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_tree_copy(const Mutation_Annotated_Tree::Tree& tree, const std::string& identifier) {
    TIMEIT();
    auto root = tree.root;
    if (identifier != "") {
        root = tree.get_node(identifier);
    }

    Tree copy = create_tree_from_newick_string (get_newick_string(tree, root, true, true));

    std::vector<Node*> dfs1;
    std::vector<Node*> dfs2;

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, 2),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t k=r.begin(); k<r.end(); ++k) {
            if (k==0) {
                dfs1 = tree.depth_first_expansion(root);
            } else {
                dfs2 = copy.depth_first_expansion();
            }
        }
    }, ap);


    tbb::parallel_for(tbb::blocked_range<size_t>(0, dfs1.size()),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t k=r.begin(); k<r.end(); ++k) {
            auto n1 = dfs1[k];
            auto n2 = dfs2[k];
            n2->clade_annotations.resize(n1->clade_annotations.size());
            for (size_t i=0; i<n1->clade_annotations.size(); i++) {
                n2->clade_annotations[i] = n1->clade_annotations[i];
            }
            for (auto m: n1->mutations) {
                Mutation m2 = m.copy();
                n2->add_mutation(m2);
            }
        }
    }, ap);

    size_t num_condensed_nodes = static_cast<size_t>(tree.condensed_nodes.size());
    tbb::parallel_for( tbb::blocked_range<size_t>(0, num_condensed_nodes),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t idx = r.begin(); idx < r.end(); idx++) {
            auto cn = tree.condensed_nodes.begin();
            std::advance(cn, idx);
            copy.condensed_nodes.insert(std::pair<std::string, std::vector<std::string>>(cn->first, std::vector<std::string>(cn->second.size())));
            for (size_t k = 0; k < cn->second.size(); k++) {
                copy.condensed_nodes[cn->first][k] = cn->second[k];
                copy.condensed_leaves.insert(cn->second[k]);
            }
        }
    }, ap);

    return copy;
}

// Get the last common ancestor of two node identifiers. Return NULL if does not
// exist
Mutation_Annotated_Tree::Node* Mutation_Annotated_Tree::LCA (const Mutation_Annotated_Tree::Tree& tree, const std::string& nid1, const std::string& nid2) {
    TIMEIT();

    if ((tree.get_node(nid1) == NULL) || (tree.get_node(nid2) == NULL)) {
        return NULL;
    }

    auto n2_ancestors = tree.rsearch(nid2, true);

    for (auto anc1: tree.rsearch(nid1, true)) {
        for (auto anc2: n2_ancestors) {
            if (anc1 == anc2) {
                return anc1;
            }
        }
    }

    return NULL;
}

// Extract the subtree consisting of the specified set of samples. This routine
// maintains the internal node names of the input tree. Mutations are copied
// from the tree such that the path of mutations from root to the sample is
// same as the original tree.
Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::get_subtree (const Mutation_Annotated_Tree::Tree& tree, const std::vector<std::string>& samples, bool keep_clade_annotations) {
    TIMEIT();
    Tree subtree;

    // Set of leaf and internal nodes corresponding to the subtree
    tbb::concurrent_unordered_set<Node*> subtree_nodes;
    // Maintain a set of all ancestors of a sample for each sample
    std::vector<tbb::concurrent_unordered_set<Node*>> all_ancestors(samples.size());

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t k=r.begin(); k<r.end(); ++k) {
            subtree_nodes.insert(tree.get_node(samples[k]));
            for (auto anc: tree.rsearch(samples[k], true)) {
                all_ancestors[k].insert(anc);
            }
        }
    }, ap);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, samples.size()),
    [&](tbb::blocked_range<size_t> r) {
        for (size_t i=r.begin(); i<r.end(); ++i) {
            for (size_t j=i+1; j<samples.size(); ++j) {
                for (auto anc: tree.rsearch(samples[i], true)) {
                    if (all_ancestors[j].find(anc) != all_ancestors[j].end()) {
                        subtree_nodes.insert(anc);
                        break;
                    }
                }
            }
        }
    }, ap);

    auto dfs = tree.depth_first_expansion();
    size_t num_annotations = 0;
    if (keep_clade_annotations) {
        num_annotations = tree.get_num_annotations();
    }

    std::stack<Node*> last_subtree_node;
    for (auto n: dfs) {
        // If the node is in subtree_nodes, it should be added to the subtree
        if (subtree_nodes.find(n) != subtree_nodes.end()) {
            Node* subtree_parent = NULL;
            if (last_subtree_node.size() > 0) {
                while (!tree.is_ancestor(last_subtree_node.top()->identifier, n->identifier)) {
                    last_subtree_node.pop();
                }
                subtree_parent = last_subtree_node.top();
            }
            // Add as root of the subtree
            if (subtree_parent == NULL) {
                // for root node, need to size the annotations vector
                Node* new_node = subtree.create_node(n->identifier, -1.0, num_annotations);
                // need to assign any clade annotations which would belong to that root as well
                for (size_t k = 0; k < num_annotations; k++) {
                    if (n->clade_annotations[k] != "") {
                        new_node->clade_annotations[k] = n->clade_annotations[k];
                    }
                }
                std::vector<Node*> root_to_node = tree.rsearch(n->identifier, true);
                std::reverse(root_to_node.begin(), root_to_node.end());
                //root_to_node.emplace_back(n);

                for (auto curr: root_to_node) {
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            // Add to the parent identified
            else {
                Node* new_node = subtree.create_node(n->identifier, subtree_parent->identifier, num_annotations);

                auto par_to_node = tree.rsearch(n->identifier, true);
                std::reverse(par_to_node.begin(), par_to_node.end());
                par_to_node.erase(par_to_node.begin(), std::find(par_to_node.begin(), par_to_node.end(), subtree_parent)+1);


                for (auto curr: par_to_node) {
                    for (size_t k = 0; k < num_annotations; k++) {
                        if (curr->clade_annotations[k] != "") {
                            new_node->clade_annotations[k] = curr->clade_annotations[k];
                        }
                    }
                    for (auto m: curr->mutations) {
                        new_node->add_mutation(m);
                    }
                }
            }
            last_subtree_node.push(n);
        }
    }

    subtree.curr_internal_node = tree.curr_internal_node;

    return subtree;
}

void Mutation_Annotated_Tree::clear_tree(Mutation_Annotated_Tree::Tree& T) {
    for (auto n: T.depth_first_expansion()) {
        delete(n);
    }
}

void Mutation_Annotated_Tree::get_random_single_subtree (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len) {
    //timer.Start();
    std::string preid = "/";
    if (use_tree_idx) {
        preid = "/tree-" + std::to_string(tree_idx) + "-";
    }
    std::set<Mutation_Annotated_Tree::Node*> leaves_to_keep_set;
    for (auto s: samples) {
        leaves_to_keep_set.insert(T->get_node(s));
    }

    auto all_leaves = T->get_leaves();
    for (size_t i=0; i< all_leaves.size(); i++) {
        auto l = all_leaves.begin();
        std::advance(l, std::rand() % all_leaves.size());
        leaves_to_keep_set.insert(*l);
        if (leaves_to_keep_set.size() >= subtree_size + samples.size()) {
            break;
        }
    }

    std::vector<std::string> leaves_to_keep;
    for (auto l: leaves_to_keep_set) {
        leaves_to_keep.emplace_back(l->identifier);
    }

    auto new_T = Mutation_Annotated_Tree::get_subtree(*T, leaves_to_keep);

    // Rotate tree for display
    new_T.rotate_for_display();

    // Write subtree to file
    auto subtree_filename = outdir + preid + "single-subtree.nh";
    fprintf(stderr, "Writing single subtree with %zu randomly added leaves to file %s.\n", subtree_size, subtree_filename.c_str());

    std::ofstream subtree_file(subtree_filename.c_str(), std::ofstream::out);
    std::stringstream newick_ss;
    write_newick_string(newick_ss, new_T, new_T.root, true, true, retain_original_branch_len);
    subtree_file << newick_ss.rdbuf();
    subtree_file.close();


    // Write list of mutations on the subtree to file
    auto subtree_mutations_filename = outdir + preid + "single-subtree-mutations.txt";
    fprintf(stderr, "Writing list of mutations at the nodes of the single subtree to file %s\n", subtree_mutations_filename.c_str());
    FILE* subtree_mutations_file = fopen(subtree_mutations_filename.c_str(), "w");

    for (auto n: new_T.depth_first_expansion()) {
        size_t tot_mutations = n->mutations.size();
        fprintf(subtree_mutations_file, "%s: ", n->identifier.c_str());
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
        if (T->condensed_nodes.find(l->identifier) != T->condensed_nodes.end()) {
            if (!has_condensed) {

                auto subtree_expanded_filename = outdir + preid + "single-subtree-expanded.txt";
                fprintf(stderr, "Subtree has condensed nodes.\nExpanding the condensed nodes for the single subtree in file %s\n", subtree_expanded_filename.c_str());
                subtree_expanded_file = fopen(subtree_expanded_filename.c_str(), "w");
                has_condensed = true;
            }
            fprintf(subtree_expanded_file, "%s: ", l->identifier.c_str());
            for (auto n: T->condensed_nodes[l->identifier]) {
                fprintf(subtree_expanded_file, "%s ", n.c_str());
            }
            fprintf(subtree_expanded_file, "\n");
        }
    }
    if (has_condensed) {
        fclose(subtree_expanded_file);
    }
}

void Mutation_Annotated_Tree::get_random_sample_subtrees (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len) {
    fprintf(stderr, "Computing subtrees for %ld samples. \n\n", samples.size());
    std::string preid = "/";
    if (use_tree_idx) {
        preid = "/tree-" + std::to_string(tree_idx) + "-";
    }
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
    auto all_leaves = T->get_leaves();
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
    for (size_t ms_idx = 0; ms_idx < samples.size(); ms_idx++) {
        if (T->get_node(samples[ms_idx]) == NULL) {
            displayed_samples[ms_idx] = true;
        }
    }

    int num_subtrees = 0;
    for (size_t i = 0; i < samples.size(); i++) {

        if (displayed_samples[i]) {
            continue;
        }

        Mutation_Annotated_Tree::Node* last_anc = T->get_node(samples[i]);
        std::vector<std::string> leaves_to_keep;

        // Keep moving up the tree till a subtree of required size is
        // found
        for (auto anc: T->rsearch(samples[i], true)) {
            size_t num_leaves = T->get_num_leaves(anc);
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

                for (auto l: T->get_leaves(last_anc->identifier)) {
                    leaves_to_keep.emplace_back(l->identifier);
                }

                std::vector<NodeDist> node_distances;
                for (auto l: T->get_leaves(anc->identifier)) {
                    if (T->is_ancestor(last_anc->identifier, l->identifier)) {
                        continue;
                    }

                    uint32_t dist = 0;
                    for (auto a: T->rsearch(l->identifier, true)) {
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
                    leaves_to_keep.emplace_back(n.node->identifier);
                }

                if ((nearest_subtree_size < subtree_size) && (nearest_subtree_size < node_distances.size())) {
                    std::vector<NodeDist> remaining_node_distances = {node_distances.begin()+nearest_subtree_size, node_distances.end()};
                    std::shuffle(remaining_node_distances.begin(), remaining_node_distances.end(), std::default_random_engine {});

                    for (auto n: remaining_node_distances) {
                        if (leaves_to_keep.size() == subtree_size) {
                            break;
                        }
                        leaves_to_keep.emplace_back(n.node->identifier);
                    }
                }
            } else {
                for (auto l: T->get_leaves(anc->identifier)) {
                    if (leaves_to_keep.size() == subtree_size) {
                        break;
                    }
                    leaves_to_keep.emplace_back(l->identifier);
                }
            }

            auto new_T = Mutation_Annotated_Tree::get_subtree(*T, leaves_to_keep);

            // Rotate tree for display
            new_T.rotate_for_display();

            tbb::parallel_for (tbb::blocked_range<size_t>(i+1, samples.size(), 100),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t j=r.begin(); j<r.end(); ++j) {
                    if (!displayed_samples[j]) {
                        if (new_T.get_node(samples[j]) != NULL) {
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
            write_newick_string(newick_ss, new_T, new_T.root, true, true, retain_original_branch_len);
            subtree_file << newick_ss.rdbuf();
            subtree_file.close();


            // Write list of mutations on the subtree to file
            auto subtree_mutations_filename = outdir + preid + "subtree-" + std::to_string(num_subtrees) + "-mutations.txt";

            fprintf(stderr, "Writing list of mutations at the nodes of subtree %d to file %s\n", num_subtrees, subtree_mutations_filename.c_str());
            FILE* subtree_mutations_file = fopen(subtree_mutations_filename.c_str(), "w");

            for (auto n: new_T.depth_first_expansion()) {
                size_t tot_mutations = n->mutations.size();
                fprintf(subtree_mutations_file, "%s: ", n->identifier.c_str());
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
                if (T->condensed_nodes.find(l->identifier) != T->condensed_nodes.end()) {
                    if (!has_condensed) {
                        auto subtree_expanded_filename = outdir + preid +  "subtree-" + std::to_string(num_subtrees) + "-expanded.txt";
                        fprintf(stderr, "Subtree %d has condensed nodes.\nExpanding the condensed nodes for subtree %d in file %s\n", num_subtrees, num_subtrees, subtree_expanded_filename.c_str());
                        subtree_expanded_file = fopen(subtree_expanded_filename.c_str(), "w");
                        has_condensed = true;
                    }
                    fprintf(subtree_expanded_file, "%s: ", l->identifier.c_str());
                    for (auto n: T->condensed_nodes[l->identifier]) {
                        fprintf(subtree_expanded_file, "%s ", n.c_str());
                    }
                    fprintf(subtree_expanded_file, "\n");
                }
            }
            if (has_condensed) {
                fclose(subtree_expanded_file);
            }
            break;
        }
    }
}

void Mutation_Annotated_Tree::get_sample_mutation_paths (Mutation_Annotated_Tree::Tree* T, std::vector<std::string> samples, std::string mutation_paths_filename) {
    FILE* mutation_paths_file = fopen(mutation_paths_filename.c_str(), "w");

    for (size_t s=0; s<samples.size(); s++) {
        auto sample = samples[s];
        auto sample_node = T->get_node(sample);

        // If the missing sample is not found in the tree, it was not placed
        // because of max_uncertainty.
        if (T->get_node(sample) == NULL) {
            continue;
        }

        // Stack for last-in first-out ordering
        std::stack<std::string> mutation_stack;
        std::string curr_node_mutation_string;

        // Mutations on the added sample
        auto curr_node_mutations = sample_node->mutations;
        if (curr_node_mutations.size() > 0) {
            curr_node_mutation_string = sample + ":";
            size_t num_mutations = curr_node_mutations.size();
            for (size_t k = 0; k < num_mutations; k++) {
                curr_node_mutation_string += curr_node_mutations[k].get_string();
                if (k < num_mutations-1) {
                    curr_node_mutation_string += ',';
                } else {
                    curr_node_mutation_string += ' ';
                }
            }
            mutation_stack.push(curr_node_mutation_string);
        }

        // Mutations on the ancestors of added sample
        for (auto anc_node: T->rsearch(sample)) {
            curr_node_mutations = anc_node->mutations;
            if (curr_node_mutations.size() > 0) {
                curr_node_mutation_string = anc_node->identifier + ":";
                size_t num_mutations = curr_node_mutations.size();
                for (size_t k = 0; k < num_mutations; k++) {
                    curr_node_mutation_string += curr_node_mutations[k].get_string();
                    if (k < num_mutations-1) {
                        curr_node_mutation_string += ',';
                    } else {
                        curr_node_mutation_string += ' ';
                    }
                }
                mutation_stack.push(curr_node_mutation_string);
            }
        }

        fprintf(mutation_paths_file, "%s\t", sample.c_str());
        while (mutation_stack.size()) {
            fprintf(mutation_paths_file, "%s", mutation_stack.top().c_str());
            mutation_stack.pop();
        }
        fprintf(mutation_paths_file, "\n");
    }
    fclose(mutation_paths_file);
}

void Mutation_Annotated_Tree::read_vcf(Mutation_Annotated_Tree::Tree* T, std::string &vcf_filename, std::vector<Missing_Sample>& missing_samples, bool create_new_mat) {
    if (create_new_mat) {
        // If called with a tree file that needs to create a MAT

        // Vector used to store all tree nodes in breadth-first search (BFS) order
        std::vector<Node*> bfs;
        // Map the node identifier string to index in the BFS traversal
        std::unordered_map<std::string, size_t> bfs_idx;

        // Variables below used to store the different fields of the input VCF file
        bool header_found = false;
        std::vector<std::string> variant_ids;

        // timer object to be used to measure runtimes of individual stages
        Timer timer;

        // Breadth-first expansion to populate bfs and bfs_idx
        bfs = T->breadth_first_expansion();
        for (size_t idx = 0; idx < bfs.size(); idx++) {
            bfs_idx[bfs[idx]->identifier] = idx;
        }

        fprintf(stderr, "Loading VCF file.\n");
        timer.Start();

        // Boost library used to stream the contents of the input VCF file in
        // uncompressed or compressed .gz format
        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", vcf_filename.c_str());
            exit(1);
        }
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }

        fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

        fprintf(stderr, "Computing parsimonious assignments for input variants.\n");
        timer.Start();

        // A TBB flow graph containing a single source_node (reader) connected
        // to several mappers. The source_node sequentially reads in the different
        // lines of the input VCF file and constructs a mapper task for each
        // VCF line. Each mapper task takes a mapper_input as input, which stores
        // the alternate alleles, ambiguous bases and missing data (Ns) for
        // different tree samples at the corresponding VCF line/position. The
        // mappers use Fitch-Sankoff algorithm to assign mutations at different
        // branches of the tree and update the mutation-annotated tree (T)
        // accordingly.
        tbb::flow::graph mapper_graph;

        tbb::flow::function_node<mapper_input, int> mapper(mapper_graph, tbb::flow::unlimited, mapper_body());
        tbb::flow::source_node <mapper_input> reader (mapper_graph,
        [&] (mapper_input &inp) -> bool {

            //check if reached end-of-file
            int curr_char = instream.peek();
            if(curr_char == EOF)
                return false;

            std::string s;
            std::getline(instream, s);
            std::vector<std::string> words;
            string_split(s, words);
            inp.variant_pos = -1;

            // Header is found when "POS" is the second word in the line
            if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                    // Sample names start from the 10th word in the header
                    for (size_t j=9; j < words.size(); j++) {
                        variant_ids.emplace_back(words[j]);
                        // If sample name not in tree, add it to missing_samples
                        if (bfs_idx.find(words[j]) == bfs_idx.end()) {
                            missing_samples.emplace_back(Missing_Sample(words[j]));
                        }
                    }
                    header_found = true;
                }
            } else if (header_found) {
                if (words.size() != 9+variant_ids.size()) {
                    fprintf(stderr, "ERROR! Incorrect VCF format.\n");
                    exit(1);
                }
                std::vector<std::string> alleles;
                alleles.clear();
                inp.variant_pos = std::stoi(words[1]);
                string_split(words[4], ',', alleles);
                // T will be modified by the mapper with mutation
                // annotations
                inp.T = T;
                inp.chrom = words[0];
                inp.bfs = &bfs;
                inp.bfs_idx = &bfs_idx;
                inp.variant_ids = &variant_ids;
                inp.missing_samples = &missing_samples;
                // Ref nuc id uses one-hot encoding (A:0b1, C:0b10, G:0b100,
                // T:0b1000)
                inp.ref_nuc = get_nuc_id(words[3][0]);
                assert((inp.ref_nuc & (inp.ref_nuc-1)) == 0); //check if it is power of 2
                inp.variants.clear();
                for (size_t j=9; j < words.size(); j++) {
                    if (isdigit(words[j][0])) {
                        int allele_id = std::stoi(words[j]);
                        if (allele_id > 0) {
                            std::string allele = alleles[allele_id-1];
                            inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id(allele[0])));
                        }
                    } else {
                        inp.variants.emplace_back(std::make_tuple(j-9, get_nuc_id('N')));
                    }
                }
            }
            return true;
        }, true );
        tbb::flow::make_edge(reader, mapper);
        mapper_graph.wait_for_all();
    } else {
        // Read vcf with existing mat

        fprintf(stderr, "Loading VCF file\n");
        // Boost library used to stream the contents of the input VCF file in
        // uncompressed or compressed .gz format
        std::ifstream infile(vcf_filename, std::ios_base::in | std::ios_base::binary);
        if (!infile) {
            fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", vcf_filename.c_str());
            exit(1);
        }
        boost::iostreams::filtering_istream instream;
        try {
            if (vcf_filename.find(".gz\0") != std::string::npos) {
                instream.push(boost::iostreams::gzip_decompressor());
            }
            instream.push(infile);
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }

        bool header_found = false;
        std::vector<std::string> variant_ids;
        std::vector<size_t> missing_idx;
        std::string s;
        // This while loop reads the VCF file line by line and populates
        // missing_samples and missing_sample_mutations based on the names and
        // variants of missing samples. If a sample name in the VCF is already
        // found in the tree, it gets ignored with a warning message
        while (instream.peek() != EOF) {
            std::getline(instream, s);
            std::vector<std::string> words;
            string_split(s, words);
            if ((not header_found) && (words.size() > 1)) {
                if (words[1] == "POS") {
                    for (size_t j=9; j < words.size(); j++) {
                        variant_ids.emplace_back(words[j]);
                        if ((T->get_node(words[j]) == NULL) && (T->condensed_leaves.find(words[j]) == T->condensed_leaves.end())) {
                            missing_samples.emplace_back(Missing_Sample(words[j]));
                            missing_idx.emplace_back(j);
                        } else {
                            fprintf(stderr, "WARNING: Ignoring sample %s as it is already in the tree.\n", words[j].c_str());
                        }
                    }
                    header_found = true;
                }
            } else if (header_found) {
                if (words.size() != 9+variant_ids.size()) {
                    fprintf(stderr, "ERROR! Incorrect VCF format. Expected %zu columns but got %zu.\n", 9+variant_ids.size(), words.size());
                    exit(1);
                }
                std::vector<std::string> alleles;
                alleles.clear();
                string_split(words[4], ',', alleles);
                for (size_t k = 0; k < missing_idx.size(); k++) {
                    size_t j = missing_idx[k];
                    auto iter = missing_samples.begin();
                    std::advance(iter, k);
                    if (iter != missing_samples.end()) {
                        Mutation m;
                        m.chrom = words[0];
                        m.position = std::stoi(words[1]);
                        m.ref_nuc = get_nuc_id(words[3][0]);
                        assert((m.ref_nuc & (m.ref_nuc-1)) == 0); //check if it is power of 2
                        m.par_nuc = m.ref_nuc;
                        // Alleles such as '.' should be treated as missing
                        // data. if the word is numeric, it is an index to one
                        // of the alleles
                        if (isdigit(words[j][0])) {
                            int allele_id = std::stoi(words[j]);
                            if (allele_id > 0) {
                                std::string allele = alleles[allele_id-1];
                                if (allele[0] == 'N') {
                                    m.is_missing = true;
                                    m.mut_nuc = get_nuc_id('N');
                                } else {
                                    auto nuc = get_nuc_id(allele[0]);
                                    if (nuc == get_nuc_id('N')) {
                                        m.is_missing = true;
                                    } else {
                                        m.is_missing = false;
                                    }
                                    m.mut_nuc = nuc;
                                }
                                (*iter).mutations.emplace_back(m);
                            }
                        } else {
                            m.is_missing = true;
                            m.mut_nuc = get_nuc_id('N');
                            (*iter).mutations.emplace_back(m);
                        }
                        if ((m.mut_nuc & (m.mut_nuc-1)) !=0) {
                            (*iter).num_ambiguous++;
                        }
                    }
                }
            }
        }
    }
}
