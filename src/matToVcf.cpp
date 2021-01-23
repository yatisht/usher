#define matToVCF
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp> 
#include <boost/filesystem.hpp>
#include "mutation_annotated_tree.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;
using namespace std;

po::variables_map check_options(int argc, char** argv) {
    // Check command line options and return variable map.
    po::options_description desc{"Options"};
    desc.add_options()
        ("mat,i", po::value<string>()->required(),
         "Mutation-annotated tree file to convert to VCF [REQUIRED]")
        ("vcf,v", po::value<string>()->required(),
         "Output VCF file [REQUIRED]")
        ("tree,t", po::value<string>()->default_value(""),
         "Output tree file")
        ("outdir,d", po::value<string>()->default_value("."),
         "Output directory to dump output and log files [DEFAULT uses current directory]")
        ("no-genotypes,n", "Do not include sample genotype columns in VCF output")
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
    catch(exception &e){
        cerr << desc << endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
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
    string mat_filename = vm["mat"].as<string>();
    string vcf_filename = vm["vcf"].as<string>();
    string tree_filename = vm["tree"].as<string>();
    string outdir = vm["outdir"].as<string>();
    bool no_genotypes = vm.count("no-genotypes");

    // Create outdir if it does not exist
    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory %s.\n", outdir.c_str());
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();

    // Load the tree and uncondense if it has condensed nodes.
    auto T = MAT::load_mutation_annotated_tree(mat_filename);
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }

    // Save VCF
    auto vcf_filepath = outdir + "/" + vcf_filename;
    FILE *vcf_file = fopen(vcf_filepath.c_str(), "w");
    vector<Mutation_Annotated_Tree::Node*> dfs = T.depth_first_expansion();
    write_vcf_header(vcf_file, dfs, !no_genotypes);
    write_vcf_rows(vcf_file, T, dfs, !no_genotypes);
    fclose(vcf_file);

    // Save Newick if specified
    if (tree_filename != "") {
        auto tree_filepath = outdir + "/" + tree_filename;
        FILE *tree_file = fopen(tree_filepath.c_str(), "w");
        fprintf(tree_file, "%s\n",
                MAT::get_newick_string(T, true, true, true).c_str());
        fclose(tree_file);
    }
}
