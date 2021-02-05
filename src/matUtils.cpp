#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
//#include "usher_graph.hpp"
#include "usher_mapper.cpp"//it shouldn't need this, but if I just include the hpp it insists that the mapper2_body function isn't defined. 

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

po::variables_map check_options(int argc, char** argv) {
    // Check command line options and return variable map.
    po::options_description desc{"Options"};
    desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("output-mat,o", po::value<std::string>()->default_value(""),
         "Use to output a full processed mutation-annotated tree file.")
        ("restricted-samples,s", po::value<std::string>()->default_value(""), 
         "Sample names to restrict. Use to perform masking") //this is now optional, as the Utils may be doing things other than masking. Should still perform the same given the same commands as previous.
        ("find-epps,e", po::bool_switch(),
        "Use to calculate and store the number of equally parsimonious placements for all nodes")
        ("write-vcf,v", po::value<string>()->default_value(""),
         "Output VCF file ")
        ("no-genotypes,n", po::bool_switch(),
        "Do not include sample genotype columns in VCF output. Used only with the vcf option")
        ("write-tree,t", po::value<string>()->default_value(""),
         "Use to write a newick tree to the indicated file.")
        ("help,h", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(p)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

/*
As a general principle, I intend to write this utility such that each function is modular and self-contained and any or all of them can be called based on command line usage.
Each relevant function will both take and return a MAT object
The main() function will only contain MAT and option read in, a series of if() then function calls, and saving the tree at the end
*/

MAT::Tree restrictSamples (std::string samples_filename, MAT::Tree T) {
    // Load restricted sampl0e names from the input file and add it to the set
    std::ifstream infile(samples_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the restricted samples file: %s!\n", samples_filename.c_str());
        exit(1);
    }    
    std::unordered_set<std::string> restricted_samples;
    std::string sample;
    while (std::getline(infile, sample)) {
        if (T.get_node(sample) == NULL) {
            fprintf(stderr, "ERROR: Sample %s missing in input MAT!\n", sample.c_str());
            exit(1);
        }
        restricted_samples.insert(std::move(sample));
    }

    // Set of nodes rooted at restricted samples
    std::unordered_set<MAT::Node*> restricted_roots;
    std::unordered_map<std::string, bool> visited;
    for (auto s: restricted_samples) {
        visited[s] = false;
    }
    for (auto cn: T.breadth_first_expansion()) { 
        auto s = cn->identifier;
        if (restricted_samples.find(s) == restricted_samples.end()) {
            continue;
        }
        if (visited[s]) {
            continue;
        }
        auto curr_node = T.get_node(s);
        for (auto n: T.rsearch(s)) {
            bool found_unrestricted = false;
            for (auto l: T.get_leaves_ids(n->identifier)) {
                if (restricted_samples.find(l)  == restricted_samples.end()) {
                    found_unrestricted = true;
                    break;
                }
            }
            if (!found_unrestricted) {
                for (auto l: T.get_leaves_ids(n->identifier)) {
                    visited[l] = true;
                }
                curr_node = n;
                break;
            }
        }
        restricted_roots.insert(curr_node);
    }

    fprintf(stderr, "Restricted roots size: %zu\n\n", restricted_roots.size());

    // Map to store number of occurences of a mutation in the tree
    std::unordered_map<std::string, int> mutations_counts;
    for (auto n: T.depth_first_expansion()) {
        for (auto mut: n->mutations) {
            if (mut.is_masked()) {
                continue;
            }
            auto mut_string = mut.get_string();
            if (mutations_counts.find(mut_string) == mutations_counts.end()) {
                mutations_counts[mut_string] = 1;
            }
            else {
                mutations_counts[mut_string] += 1;
            }
        }
    }
    
    // Reduce mutation counts for mutations in subtrees rooted at 
    // restricted_roots. Mutations specific to restricted samples 
    // will now be set to 0. 
    for (auto r: restricted_roots) {
        //fprintf(stdout, "At restricted root %s\n", r->identifier.c_str());
        for (auto n: T.depth_first_expansion(r)) {
            for (auto mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                mutations_counts[mut_string] -= 1;
            }
        }
    }

    for (auto r: restricted_roots) {
        for (auto n: T.depth_first_expansion(r)) {
            for (auto& mut: n->mutations) {
                if (mut.is_masked()) {
                    continue;
                }
                auto mut_string = mut.get_string();
                if (mutations_counts[mut_string] == 0) {
                    fprintf(stderr, "Masking mutation %s at node %s\n", mut_string.c_str(), n->identifier.c_str());
                    mut.position = -1;
                    mut.ref_nuc = 0;
                    mut.par_nuc = 0;
                    mut.mut_nuc = 0;
                }
            }
        }
    }
    return T;
}

MAT::Tree findEPPs (MAT::Tree Tobj) {
    TIMEIT()
    //all comments on function JDM

    MAT::Tree TCopy = MAT::get_tree_copy(Tobj);
    MAT::Tree* T = &TCopy; //wants the pointer for mapping. Going to be popping from and adding to the tree copy for now.

    auto fdfs = Tobj.depth_first_expansion(); //the full tree expanded for outer loop iteration.
    for (size_t s=0; s<fdfs.size(); s++){ //this loop is not a parallel for because its going to contain a parallel for
        //get the node object.
        auto node = fdfs[s];
        if (node->is_leaf()) {
            //retrieve the full set of mutations associated with this Node object from root to it
            //to do this, get the full set of ancestral nodes and their mutations
            //code copied from the usher mapper.
            std::vector<int> anc_positions; //tracking positions is required to account for backmutation/overwriting along the path
            std::vector<MAT::Mutation> ancestral_mutations;
            //first load in the current mutations
            for (auto m: node->mutations){
                if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                    ancestral_mutations.emplace_back(m);
                    if (!m.is_masked()) {
                        anc_positions.emplace_back(m.position);
                    }            
                }
            }
            //then load in ancestral mutations
            for (auto n: Tobj.rsearch(node->identifier)) {
                for (auto m: n->mutations) {
                    if (m.is_masked() || (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                        ancestral_mutations.emplace_back(m);
                        if (!m.is_masked()) {
                            anc_positions.emplace_back(m.position);
                        }
                    }
                }
            }
            //now that we have the mutation set, pop the current node from the Tree
            T->remove_node(node->identifier, true); //this should modify in-place. pop it from the copy
            //the ancestral_mutations vector, plus the mutations assigned to this specific node, constitute the "missing_sample" equivalents for calling the mapper

            auto dfs = T->depth_first_expansion();
            size_t total_nodes = dfs.size();

            // Stores the excess mutations to place the sample at each
            // node of the tree in DFS order. When placement is as a
            // child, it only contains parsimony-increasing mutations in
            // the sample. When placement is as a sibling, it contains 
            // parsimony-increasing mutations as well as the mutations
            // on the placed node in common with the new sample. Note
            // guaranteed to be corrrect only for optimal nodes since
            // the mapper can terminate the search early for non-optimal
            // nodes
            std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
            // Stores the imputed mutations for ambiguous bases in the
            // sampled in order to place the sample at each node of the 
            // tree in DFS order. Again, guaranteed to be corrrect only 
            // for pasrimony-optimal nodes 
            std::vector<std::vector<MAT::Mutation>> node_imputed_mutations(total_nodes);

            // Stores the parsimony score to place the sample at each
            // node of the tree in DFS order.
            std::vector<int> node_set_difference;
            size_t best_node_num_leaves = 0;
            // The maximum number of mutations is bound by the number
            // of mutations in the missing sample (place at root)
            //int best_set_difference = 1e9;
            // TODO: currently number of root mutations is also added to
            // this value since it forces placement as child but this
            // could be changed later 
            int best_set_difference = ancestral_mutations.size() + T->root->mutations.size() + 1;

            size_t best_j = 0;
            size_t num_best = 1;
            bool best_node_has_unique = false;
            MAT::Node* best_node = T->root;

            std::vector<bool> node_has_unique(total_nodes, false);
            std::vector<size_t> best_j_vec;
            best_j_vec.emplace_back(0);

            // Parallel for loop to search for most parsimonious
            // placements. Real action happens within mapper2_body
            auto grain_size = 400; 
            tbb::parallel_for( tbb::blocked_range<size_t>(0, total_nodes, grain_size),
                    [&](tbb::blocked_range<size_t> r) {
                    for (size_t k=r.begin(); k<r.end(); ++k){
                        mapper2_input inp;
                        inp.T = T;
                        inp.node = dfs[k];
                        inp.missing_sample_mutations = &ancestral_mutations;
                        inp.excess_mutations = &node_excess_mutations[k];
                        inp.imputed_mutations = &node_imputed_mutations[k];
                        inp.best_node_num_leaves = &best_node_num_leaves;
                        inp.best_set_difference = &best_set_difference;
                        inp.best_node = &best_node;
                        inp.best_j =  &best_j;
                        inp.num_best = &num_best;
                        inp.j = k;
                        inp.has_unique = &best_node_has_unique;
                        inp.best_j_vec = &best_j_vec;
                        inp.node_has_unique = &(node_has_unique);

                        mapper2_body(inp, false);
                    }       
                    }); 
            //BACK TO MY CODE (JDM)
            //put the node back.
            T->add_node(node, node->parent); //simplest option? assuming the node object doesn't get deleted along with the tree vector attribute

            //give the original tree node, which hasn't moved, the metadata.
            auto cnode = Tobj.get_node(node->identifier);
            cnode->epps = num_best;
            fprintf(stderr, "Node metadata updated- ID %s ", node->identifier.c_str());
            fprintf(stderr, "EPPs %ld\n", num_best);
        }
    }
    return Tobj; //return the actual object.
}

//CODE FOR VCF CONVERSION

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes) {
    auto vcf_filepath = outdir + "/" + vcf_filename;
    FILE *vcf_file = fopen(vcf_filepath.c_str(), "w");
    vector<Mutation_Annotated_Tree::Node*> dfs = T.depth_first_expansion();
    write_vcf_header(vcf_file, dfs, !no_genotypes);
    write_vcf_rows(vcf_file, T, dfs, !no_genotypes);
    fclose(vcf_file);
}

void write_vcf_header(FILE *vcf_file, vector<Mutation_Annotated_Tree::Node*> &dfs,
                      bool print_genotypes) {
    // Write minimal VCF header with sample names in same order that genotypes
    // will be printed out (DFS).
    fprintf(vcf_file, "##fileformat=VCFv4.2\n");
    fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if (print_genotypes) {
        fprintf(vcf_file, "\tFORMAT");
        for (auto node: dfs) {
            if (node->is_leaf()) {
                fprintf(vcf_file, "\t%s", node->identifier.c_str());
            }
        }
    }
    fputc('\n', vcf_file);
}

uint count_leaves(vector<Mutation_Annotated_Tree::Node*> &dfs) {
    // Return the number of leaf nodes in dfs
    uint count = 0;
    for (auto node: dfs) {
        if (node->is_leaf()) {
            count++;
        }
    }
    return count;
}

int8_t *new_gt_array(int size, int8_t ref) {
    // Allocate and return an array of int8_t (encoded nucleotide values) initialized to ref.
    int8_t *gt_array = new int8_t[size];
    for (int i = 0;  i < size;  i++) {
        gt_array[i] = ref;
    }
    return gt_array;
}

uint r_add_genotypes(MAT::Node *node,
                     unordered_map<string, vector<int8_t *>> &chrom_pos_genotypes, 
                     unordered_map<string, vector<int8_t>> &chrom_pos_ref,
                     uint leaf_count, uint leaf_ix, vector<struct MAT::Mutation *> &mut_stack) {
    // Traverse tree, adding leaf/sample genotypes for mutations annotated on path from root to node
    // to chrom_pos_genotypes (and reference allele to chrom_pos_ref).
    for (auto &mut: node->mutations) {
      if (mut.is_masked()) {
          continue;
      }
      mut_stack.push_back(&mut);
    }
    if (node->is_leaf()) {
        // Store genotypes in this leaf's column for all mutations on the path from root to leaf
        for (auto mut: mut_stack) {
            string chrom = mut->chrom;
            uint pos = (uint)mut->position;
            if (chrom.empty()) {
              fprintf(stderr, "mut->chrom is empty string at node '%s', position %u\n",
                      node->identifier.c_str(), pos);
            }
            if (chrom_pos_genotypes.find(chrom) == chrom_pos_genotypes.end()) {
                // First variant on chrom: initialize a vector mapping position to genotype.
                // Assume a genome size similar to SARS-CoV-2, resize if necessary.
                uint initSize = 30000;
                chrom_pos_genotypes[chrom] = vector<int8_t *>(initSize);
                chrom_pos_ref[chrom] = vector<int8_t>(initSize);
            }
            if (pos >= chrom_pos_genotypes[chrom].size()) {
                // chrom has larger positions than we assumed; allocate a larger vector.
                uint newSize = chrom_pos_genotypes[chrom].size() * 2;
                chrom_pos_genotypes[chrom].resize(newSize);
                chrom_pos_ref[chrom].resize(newSize);
            }
            if (! chrom_pos_genotypes[chrom][pos]) {
                // First variant reported at this position; allocate genotype array and
                // store reference allele (which is stored in par_nuc not ref_nuc).
                chrom_pos_genotypes[chrom][pos] = new_gt_array(leaf_count, mut->par_nuc);
                chrom_pos_ref[chrom][pos] = mut->par_nuc;
            }
            // Store the allele/genotype for this chrom / pos / sample.
            chrom_pos_genotypes[chrom][pos][leaf_ix] = mut->mut_nuc;
        }
        leaf_ix++;
    }
    for (auto child: node->children) {
        leaf_ix = r_add_genotypes(child, chrom_pos_genotypes, chrom_pos_ref, leaf_count, leaf_ix,
                                  mut_stack);
    }
    for (auto mut: node->mutations) {
        mut_stack.pop_back();
    }
    return leaf_ix;
}

unordered_map<int8_t, uint>count_alleles(int8_t *gt_array, uint gtCount)  {
    // Tally up the count of each allele (both ref and alts) from sample genotypes.
    unordered_map<int8_t, uint> allele_counts;
    for (uint i = 0;  i < gtCount;  i++) {
        int8_t allele = gt_array[i];
        if (allele_counts.find(allele) == allele_counts.end()) {
            allele_counts.insert({allele, 1});
        } else {
            allele_counts[allele]++;
        }
    }
    return allele_counts;
}

bool cmp_allele_count_desc(const pair<int8_t, uint>& a, const pair<int8_t, uint>& b) {
    // Compare counts of two alleles, for sorting in descending order.
    return a.second > b.second;
}

map<int8_t, uint>make_alts(unordered_map<int8_t, uint> &allele_counts, int8_t ref) {
    // Map alternate alleles, ordered by count (highest first), to counts.
    vector<pair<int8_t, uint>> pairs;
    for (auto &itr : allele_counts) {
        if (itr.first != ref) {
            pairs.push_back(itr);
        }
    }
    sort(pairs.begin(), pairs.end(), cmp_allele_count_desc);
    map<int8_t, uint> alts;
    for (auto &itr : pairs) {
      alts.insert(itr);
    }
    return alts;
}

string make_id(int8_t ref, uint pos, map<int8_t, uint> &alts) {
    // Return a C string comma-sep list of the form <ref><pos><alt1>[,<ref><pos><alt2>[,...]].
    string id;
    for (auto &itr : alts) {
        if (! id.empty()) {
            id += ",";
        }
        id += MAT::get_nuc(ref) + to_string(pos) + MAT::get_nuc(itr.first);
    }
    return id;
}

string make_alt_str(map<int8_t, uint> &alts) {
    // Return a C string comma-sep list of alternate alleles.
    string alt_str;
    for (auto &itr : alts) {
        if (! alt_str.empty()) {
          alt_str += ",";
        }
        alt_str += MAT::get_nuc(itr.first);
    }
    return alt_str;
}

string make_info(map<int8_t, uint> &alts, uint leaf_count) {
    // Return a C string VCF INFO value with AC (comma-sep list of alternate allele counts)
    // and AN (total genotype count).
    string alt_count_str;
    for (auto &itr : alts) {
        if (! alt_count_str.empty()) {
            alt_count_str += ",";
        }
        alt_count_str += to_string(itr.second);
    }
    string info = "AC=" + alt_count_str + ";AN=" + to_string(leaf_count);
    return info;
}

int *make_allele_codes(int8_t ref, map<int8_t, uint> &alts) {
    // Return an array that maps binary-encoded nucleotide to VCF genotype encoding:
    // 0 for reference allele, 1 for first alternate allele, and so on.
    int *al_codes = new int[256];
    for (int i = 0;  i < 256;  i++) {
        al_codes[i] = 0;
    }
    al_codes[(uint8_t)ref] = 0;
    int altIx = 1;
    for (auto &itr : alts) {
        al_codes[itr.first] = altIx++;
    }
    return al_codes;
}

void write_vcf_rows(FILE *vcf_file, MAT::Tree T, vector<MAT::Node*> &dfs, bool print_genotypes) {
    // Fill in a matrix of genomic positions and sample genotypes in the same order as the
    // sample names in the header, compute allele counts, and output VCF rows.
    uint leaf_count = count_leaves(dfs);
    // The int8_t here is mutation_annotated_tree.hpp's binary encoding of IUPAC nucleotide bases.
    unordered_map<string, vector<int8_t *>> chrom_pos_genotypes;
    unordered_map<string, vector<int8_t>> chrom_pos_ref;
    vector<struct MAT::Mutation *> mut_stack;
    r_add_genotypes(T.root, chrom_pos_genotypes, chrom_pos_ref, leaf_count, 0, mut_stack);
    // Write row of VCF for each variant in chrom_pos_genotypes[chrom]
    for (auto itr = chrom_pos_genotypes.begin();  itr != chrom_pos_genotypes.end();  ++itr) {
        string chrom = itr->first;
        vector<int8_t *> pos_genotypes = itr->second;
        for (uint pos = 0;  pos < pos_genotypes.size();  pos++) {
            int8_t *gt_array = pos_genotypes[pos];
            if (gt_array) {
              int8_t ref = chrom_pos_ref[chrom][pos];
              unordered_map<int8_t, uint>allele_counts = count_alleles(gt_array, leaf_count);
              map<int8_t, uint>alts = make_alts(allele_counts, ref);
              string id = make_id(ref, pos, alts);
              string alt_str = make_alt_str(alts);
              string info = make_info(alts, leaf_count);
              fprintf(vcf_file, "%s\t%d\t%s\t%c\t%s\t.\t.\t%s",
                      chrom.c_str(), pos, id .c_str(), MAT::get_nuc(ref), alt_str.c_str(),
                      info.c_str());
              if (print_genotypes) {
                  int *allele_codes = make_allele_codes(ref, alts);
                  fprintf(vcf_file, "\tGT");
                  for (uint i = 0;  i < leaf_count;  i++) {
                      int8_t allele = gt_array[i];
                      fprintf(vcf_file, "\t%d", allele_codes[allele]);
                  }
              }
              fputc('\n', vcf_file);
            }
        }
    }
}

int main(int argc, char** argv) {

    // Command line options
    po::variables_map vm = check_options(argc, argv);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    std::string samples_filename = vm["restricted-samples"].as<std::string>();
    std::string tree_filename = vm["write-tree"].as<string>();
    std::string vcf_filename = vm["write-vcf"].as<string>();

    bool fepps = vm["find-epps"].as<bool>();
    bool no_genotypes = vm['no-genotypes'].as<bool>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    // If a restricted samples file was provided, perform masking procedure
    if (samples_filename != "") {
        T = restrictSamples(samples_filename, T);
    }
    // If the argument to calculate equally parsimonious placements was used, perform this operation
    if (fepps) {
        fprintf(stderr, "Attempting to calculate EPPs\n");
        T = findEPPs(T);
    }
    //if a vcf filename was given, write a vcf to it
    if (vcf_filename != "") {
        make_vcf(T, vcf_filename, no_genotypes);
    }
    //if a newick tree filename was given, write a tree to it
    if (tree_filename != "") {
        auto tree_filepath = outdir + "/" + tree_filename;
        FILE *tree_file = fopen(tree_filepath.c_str(), "w");
        fprintf(tree_file, "%s\n",
            MAT::get_newick_string(T, true, true, true).c_str());
        fclose(tree_file);        
    }
    // Store final MAT to output file if indicated
    if (output_mat_filename != "") {
        MAT::save_mutation_annotated_tree(T, output_mat_filename);
    }

    return 0;
}

