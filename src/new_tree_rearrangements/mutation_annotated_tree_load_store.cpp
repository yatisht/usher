#include "mutation_annotated_tree.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <cstdio>
#include <fcntl.h>
#include <iostream>
#include <boost/iostreams/filter/gzip.hpp>
#include "parsimony.pb.h"
#include <stack>
#include <fstream>
#include <sstream>
#include "tbb/parallel_for.h"
#include <queue>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <google/protobuf/io/coded_stream.h>
#include <string>
std::vector<int8_t> Mutation_Annotated_Tree::get_nuc_vec_from_id (int8_t nuc_id) {
    return get_nuc_vec(get_nuc(nuc_id));
}

static void print_node(
    std::stringstream &ss, Mutation_Annotated_Tree::Node *n,
    bool print_branch_len, bool retain_original_branch_len,
    bool uncondense_leaves,
    const tbb::concurrent_unordered_map<std::string, std::vector<std::string>>
        &condensed_nodes,bool print_internal) {
    if (n->is_leaf() && uncondense_leaves &&
        (condensed_nodes.find(n->identifier) != condensed_nodes.end())) {
        auto cn = condensed_nodes.at(n->identifier);
        auto cn_size = cn.size();
        for (size_t idx = 0; idx < cn_size; idx++) {
            ss << cn[idx];
            if (idx + 1 < cn_size) {
                ss << ',';
            }
        }
    } else if(print_internal||n->is_leaf()) {
        ss << n->identifier;
    }
    float branch_length=retain_original_branch_len?n->branch_length:n->mutations.size();
    if ((print_branch_len) && (branch_length >= 0)) {
        ss << ':';
        ss << branch_length;
    }
}
// Get newick stringstream for the input subtree rooted at some node (node) in 
// the input tree T. Boolean arguments decide whether
// internal node ids and branch lengths are printed. If last boolean argument is
// set, branch lengths from input tree are retained, otherwise, branch length
// for a branch is equal to the number of mutations annotated on that branch 
void Mutation_Annotated_Tree::write_newick_string (std::stringstream& ss, const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    TIMEIT();

    struct stack_content{Node* this_node;size_t child_idx;};
    std::stack<stack_content> node_stack;
    node_stack.push({T.root,0});

    while (!node_stack.empty()) {
        //puts(ss.str().c_str());
        auto& stack_top=node_stack.top();
        bool have_children=!stack_top.this_node->children.empty();
        auto& child_idx=stack_top.child_idx;
        bool after_last=have_children&&child_idx==stack_top.this_node->children.size();
        if(after_last){
            ss<<')';
            print_node(ss,stack_top.this_node , print_branch_len,  retain_original_branch_len,  uncondense_leaves, T.condensed_nodes,print_internal);
            node_stack.pop();
            continue;
        }

        if(have_children){
            if (child_idx==0) {
                ss<<'(';
            }
            node_stack.push({stack_top.this_node->children[child_idx],0});
            child_idx++;
            if ((child_idx!=1)&&(child_idx<=stack_top.this_node->children.size())) {
                ss<<',';
            }
            continue;
        }else {
            print_node(ss,stack_top.this_node , print_branch_len,  retain_original_branch_len,  uncondense_leaves, T.condensed_nodes,print_internal);
            node_stack.pop();
            continue;
        }
    }
    ss<<';';
}

std::string Mutation_Annotated_Tree::get_newick_string (const Mutation_Annotated_Tree::Tree& T, Mutation_Annotated_Tree::Node* node, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    std::stringstream newick_ss;
    write_newick_string(newick_ss, T, node, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
    return newick_ss.str();
}

std::string Mutation_Annotated_Tree::get_newick_string (const Tree& T, bool print_internal, bool print_branch_len, bool retain_original_branch_len, bool uncondense_leaves) {
    return get_newick_string(T, T.root, print_internal, print_branch_len, retain_original_branch_len, uncondense_leaves);
}

// Split string into words for a specific delimiter delim
void Mutation_Annotated_Tree::string_split (std::string const& s, char delim, std::vector<std::string>& words) {
    TIMEIT();
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if ((end_pos == start_pos) || end_pos >= s.length()) {
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

static void populate_node(const char*& in,const char* end, Mutation_Annotated_Tree::Node* this_node, Mutation_Annotated_Tree::Tree& T,bool use_internal_node_label){
    auto start_iter=in;
    errno=0;
    char *endptr;
    float loaded=strtof(in,&endptr);
    if(errno==0){
    if (*endptr==':'||*endptr==','||*endptr==')'||*endptr==';') {
        assert(!this_node->is_leaf());
        this_node->branch_length=loaded;
        this_node->identifier=std::to_string(++T.curr_internal_node);
        T.all_nodes.emplace(this_node->identifier,this_node);
        in=endptr-1;
        return;
    }
    }
    while (in!=end&&*in!=':'&&*in!=','&&*in!=')') {
        in++;
    }
    if (use_internal_node_label||this_node->is_leaf()) {
        assert(end>start_iter);
        this_node->identifier=std::string(start_iter,in-start_iter);
    }else {
        this_node->identifier=std::to_string(++T.curr_internal_node);
    }
    T.all_nodes.emplace(this_node->identifier,this_node);
    if (*in==':') {
        in++;
    }
    errno=0;
    loaded=strtof(in,&endptr);
    if(errno==0){
this_node->branch_length=loaded;
    }
    in=endptr-1;
    
}
static void populate_node(Mutation_Annotated_Tree::Node* to_populate,const char*& start,Mutation_Annotated_Tree::Tree* tree){
    const char* name=start;
    while (*start!=':') {
        start++;
    }
    start++;
    to_populate->identifier=std::string(name,start-1);
    char* end;
    to_populate->branch_length=strtod(start,&end);
    start=end;
    tree->all_nodes.emplace(to_populate->identifier,to_populate);
}
static void load_with_internal_label(const char* start,const char* end,Mutation_Annotated_Tree::Tree* tree){
    assert(*start=='(');
    Mutation_Annotated_Tree::Node* curr_node=new Mutation_Annotated_Tree::Node();
    tree->root=curr_node;
    start++;
    while (start!=end) {
        if (*start=='(') {
            Mutation_Annotated_Tree::Node* par_node=curr_node;
            curr_node=new Mutation_Annotated_Tree::Node();
            curr_node->parent=par_node;
            par_node->children.push_back(curr_node);
            start++;
        }else if (*start==',') {
            start++;
        }else if (*start==';') {
            start++;
        }
        else if (*start==')') {
            start++;
            populate_node(curr_node,start,tree);
            curr_node=curr_node->parent;
        }else  {
            Mutation_Annotated_Tree::Node* new_node=new Mutation_Annotated_Tree::Node();
            new_node->parent=curr_node;
            curr_node->children.push_back(new_node);
            populate_node(new_node,start,tree);
        }
    }
}
void Mutation_Annotated_Tree::Tree::load_from_newick(const std::string& newick_string,bool use_internal_node_label){
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
            }
            else if (c == '(') {
                no++;
                level++;
                if (branch_len.size() <= level) {
                  branch_len.resize(level*2);
                }
            }
            else if (c == ')') {
                stop = true;
                nc++;
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                branch_len[level].push(len);
                level--;
                branch_start = false;
            }
            else if (!stop) {
                leaf += c;
                branch_start = false;
            }
            else if (branch_start) {
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

    curr_internal_node = 0;
    std::stack<Node*> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = std::to_string(++curr_internal_node);
            Node* new_node = NULL;
            if (parent_stack.size() == 0) {
                new_node = create_node(nid, branch_len[level].front());
            }
            else {
                new_node = create_node(nid, parent_stack.top(), branch_len[level].front());
            }
            branch_len[level].pop();
            level++;
            parent_stack.push(new_node);
        }
        create_node(leaf, parent_stack.top(), branch_len[level].front());
        branch_len[level].pop();
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
            level--;
        }
    }

    if (root == NULL) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

}

Mutation_Annotated_Tree::Tree Mutation_Annotated_Tree::create_tree_from_newick_string (const std::string newick_string) {
    TIMEIT();
    Tree T;
    T.load_from_newick(newick_string);
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

    struct stat stat_buf;
    stat(filename.c_str(),&stat_buf);
    size_t file_size=stat_buf.st_size;
    auto fd=open(filename.c_str(), O_RDONLY);
    uint8_t* maped_file=(uint8_t*)mmap(nullptr, file_size, PROT_READ, MAP_SHARED,fd , 0);
    close(fd);
    google::protobuf::io::CodedInputStream input(maped_file,file_size);
    input.SetTotalBytesLimit(file_size*4, file_size*4);
    data.ParseFromCodedStream(&input);
    munmap(maped_file, file_size);
    //check if the pb has a metadata field
    bool hasmeta = (data.metadata_size()>0);
    if (!hasmeta) {
        fprintf(stderr, "WARNING: This pb does not include any metadata. Filling in default values\n");
    }
    tree.load_from_newick(data.newick());
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
                  char mut_one_hot=1<<mut.mut_nuc(0);
                  char all_major_alleles=mut_one_hot;
                  if (mut.position()>0) {
                     for (int n = 1; n < mut.mut_nuc_size(); n++) {
                        all_major_alleles|= (1<<mut.mut_nuc(n));
                     }
                  }
                    Mutation m(mut.chromosome(),mut.position(),nuc_one_hot(mut_one_hot),two_bit_to_one_hot(mut.par_nuc()),all_major_alleles,two_bit_to_one_hot(mut.ref_nuc()));
                    node->add_mutation(m);
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
               }
            }
    }, ap);

    return tree;
}

void Mutation_Annotated_Tree::save_mutation_annotated_tree (const Mutation_Annotated_Tree::Tree& tree, std::string filename) {
    TIMEIT();
    Parsimony::data data;
    auto dfs = tree.depth_first_expansion();


    data.set_newick(get_newick_string(tree, true, true));
    for (size_t idx = 0; idx < dfs.size(); idx++) {
        auto meta = data.add_metadata();
        for (size_t k = 0; k < dfs[idx]->clade_annotations.size(); k++) {
            meta->add_clade_annotations(dfs[idx]->clade_annotations[k]);
        }
        auto mutation_list = data.add_node_mutations();
        for (auto m: dfs[idx]->mutations) {

            auto mut = mutation_list->add_mutation();
            mut->set_chromosome(m.get_chromosome());
            mut->set_position(m.get_position());
            
            if (m.is_masked()) {
                mut->set_ref_nuc(-1);
                mut->set_par_nuc(-1);
            }
            else {
                int8_t j = one_hot_to_two_bit(m.get_ref_one_hot()) ;
                assert (j >= 0);
                mut->set_ref_nuc(j);

                j = one_hot_to_two_bit(m.get_par_one_hot()) ;
                assert(j >= 0);
                mut->set_par_nuc(j);

                mut->clear_mut_nuc();
                mut->add_mut_nuc(one_hot_to_two_bit(m.get_mut_one_hot()));
                /*if (dfs[idx]->is_leaf()) {
                    nuc_one_hot other_mut=m.get_all_major_allele()&(~m.get_mut_one_hot());
                    if (other_mut) {
                        for (int i=0; i<4; i++) {
                            if ((1<<i)&other_mut) {
                                mut->add_mut_nuc(i);
                            }
                        }
                    }
                }*/
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
        }
        catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    }
    else {
        data.SerializeToOstream(&outfile);
        outfile.close();
    }
}

Mutation_Annotated_Tree::Mutation::Mutation(const std::string &chromosome,
                                            int position, nuc_one_hot mut,
                                            nuc_one_hot par, nuc_one_hot tie,
                                            nuc_one_hot ref)
    : position(position), par_mut_nuc((par << 4) | (mut)),
      boundary1_all_major_allele(tie) {
    auto ins_result = chromosome_map.emplace(chromosome, chromosome_map.size());
    if (ins_result.second) {
        std::lock_guard<std::mutex> lk(ref_lock);
        chromosomes.push_back(chromosome);
    }
    chrom_idx = ins_result.first->second;
    if (ref) {
        std::lock_guard<std::mutex> lk(ref_lock);
        refs.resize(std::max((int)refs.size(), position + 1));
        assert((refs[position].is_invalid()) || refs[position] == ref);
        refs[position] = ref;
    }
};
