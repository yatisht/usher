#include "convert.hpp"

po::variables_map parse_convert_command(po::parsed_options parsed) {

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description conv_desc("convert options");
    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("write-vcf,v", po::value<std::string>()->default_value(""),
         "Output VCF file ")
        ("no-genotypes,n", po::bool_switch(),
        "Do not include sample genotype columns in VCF output. Used only with the vcf option")
        ("write-tree,t", po::value<std::string>()->default_value(""),
         "Use to write a newick tree to the indicated file.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void convert_main(po::parsed_options parsed) {
    //the convert subcommand converts the protobuf into another file format, returning some other file type that represents an equivalent structure
    po::variables_map vm = parse_convert_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string tree_filename = vm["write-tree"].as<std::string>();
    std::string vcf_filename = vm["write-vcf"].as<std::string>();
    bool no_genotypes = vm["no-genotypes"].as<bool>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    //if a vcf filename was given, write a vcf to it
    if (vcf_filename != "") {
        fprintf(stderr, "Generating VCF\n");
        make_vcf(T, vcf_filename, no_genotypes);
    }
    //if a newick tree filename was given, write a tree to it
    if (tree_filename != "") {
        fprintf(stderr, "Generating Newick file\n");
        auto tree_filepath = "./" + tree_filename; //for simplicity, write it to this directory
        FILE *tree_file = fopen(tree_filepath.c_str(), "w");
        fprintf(tree_file, "%s\n",
            MAT::get_newick_string(T, true, true, true).c_str());
        fclose(tree_file);        
    }
}

void write_vcf_header(FILE *vcf_file, std::vector<Mutation_Annotated_Tree::Node*> &dfs,
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

uint count_leaves(std::vector<Mutation_Annotated_Tree::Node*> &dfs) {
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
                     std::unordered_map<std::string, std::vector<int8_t *>> &chrom_pos_genotypes, 
                     std::unordered_map<std::string, std::vector<int8_t>> &chrom_pos_ref,
                     uint leaf_count, uint leaf_ix, std::vector<struct MAT::Mutation *> &mut_stack) {
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
            std::string chrom = mut->chrom;
            uint pos = (uint)mut->position;
            if (chrom.empty()) {
              fprintf(stderr, "mut->chrom is empty std::string at node '%s', position %u\n",
                      node->identifier.c_str(), pos);
            }
            if (chrom_pos_genotypes.find(chrom) == chrom_pos_genotypes.end()) {
                // First variant on chrom: initialize a vector mapping position to genotype.
                // Assume a genome size similar to SARS-CoV-2, resize if necessary.
                uint initSize = 30000;
                chrom_pos_genotypes[chrom] = std::vector<int8_t *>(initSize);
                chrom_pos_ref[chrom] = std::vector<int8_t>(initSize);
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

std::unordered_map<int8_t, uint>count_alleles(int8_t *gt_array, uint gtCount)  {
    // Tally up the count of each allele (both ref and alts) from sample genotypes.
    std::unordered_map<int8_t, uint> allele_counts;
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

bool cmp_allele_count_desc(const std::pair<int8_t, uint>& a, const std::pair<int8_t, uint>& b) {
    // Compare counts of two alleles, for sorting in descending order.
    return a.second > b.second;
}

std::map<int8_t, uint>make_alts(std::unordered_map<int8_t, uint> &allele_counts, int8_t ref) {
    // Map alternate alleles, ordered by count (highest first), to counts.
    std::vector<std::pair<int8_t, uint>> pairs;
    for (auto &itr : allele_counts) {
        if (itr.first != ref) {
            pairs.push_back(itr);
        }
    }
    std::sort(pairs.begin(), pairs.end(), cmp_allele_count_desc);
    std::map<int8_t, uint> alts;
    for (auto &itr : pairs) {
      alts.insert(itr);
    }
    return alts;
}

std::string make_id(int8_t ref, uint pos, std::map<int8_t, uint> &alts) {
    // Return a C std::string comma-sep list of the form <ref><pos><alt1>[,<ref><pos><alt2>[,...]].
    std::string id;
    for (auto &itr : alts) {
        if (! id.empty()) {
            id += ",";
        }
        id += MAT::get_nuc(ref) + std::to_string(pos) + MAT::get_nuc(itr.first);
    }
    return id;
}

std::string make_alt_str(std::map<int8_t, uint> &alts) {
    // Return a C std::string comma-sep list of alternate alleles.
    std::string alt_str;
    for (auto &itr : alts) {
        if (! alt_str.empty()) {
          alt_str += ",";
        }
        alt_str += MAT::get_nuc(itr.first);
    }
    return alt_str;
}

std::string make_info(std::map<int8_t, uint> &alts, uint leaf_count) {
    // Return a C std::string VCF INFO value with AC (comma-sep list of alternate allele counts)
    // and AN (total genotype count).
    std::string alt_count_str;
    for (auto &itr : alts) {
        if (! alt_count_str.empty()) {
            alt_count_str += ",";
        }
        alt_count_str += std::to_string(itr.second);
    }
    std::string info = "AC=" + alt_count_str + ";AN=" + std::to_string(leaf_count);
    return info;
}

int *make_allele_codes(int8_t ref, std::map<int8_t, uint> &alts) {
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

void write_vcf_rows(FILE *vcf_file, MAT::Tree T, std::vector<MAT::Node*> &dfs, bool print_genotypes) {
    // Fill in a matrix of genomic positions and sample genotypes in the same order as the
    // sample names in the header, compute allele counts, and output VCF rows.
    uint leaf_count = count_leaves(dfs);
    // The int8_t here is mutation_annotated_tree.hpp's binary encoding of IUPAC nucleotide bases.
    std::unordered_map<std::string, std::vector<int8_t *>> chrom_pos_genotypes;
    std::unordered_map<std::string, std::vector<int8_t>> chrom_pos_ref;
    std::vector<struct MAT::Mutation *> mut_stack;
    r_add_genotypes(T.root, chrom_pos_genotypes, chrom_pos_ref, leaf_count, 0, mut_stack);
    // Write row of VCF for each variant in chrom_pos_genotypes[chrom]
    for (auto itr = chrom_pos_genotypes.begin();  itr != chrom_pos_genotypes.end();  ++itr) {
        std::string chrom = itr->first;
        std::vector<int8_t *> pos_genotypes = itr->second;
        for (uint pos = 0;  pos < pos_genotypes.size();  pos++) {
            int8_t *gt_array = pos_genotypes[pos];
            if (gt_array) {
              int8_t ref = chrom_pos_ref[chrom][pos];
              std::unordered_map<int8_t, uint>allele_counts = count_alleles(gt_array, leaf_count);
              std::map<int8_t, uint>alts = make_alts(allele_counts, ref);
              std::string id = make_id(ref, pos, alts);
              std::string alt_str = make_alt_str(alts);
              std::string info = make_info(alts, leaf_count);
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

void make_vcf (MAT::Tree T, std::string vcf_filename, bool no_genotypes) {
    auto vcf_filepath = "./" + vcf_filename; //taking away a touch of functionality for simplicity. just write it to this directory.
    FILE *vcf_file = fopen(vcf_filepath.c_str(), "w");
    std::vector<Mutation_Annotated_Tree::Node*> dfs = T.depth_first_expansion();
    write_vcf_header(vcf_file, dfs, !no_genotypes);
    write_vcf_rows(vcf_file, T, dfs, !no_genotypes);
    fclose(vcf_file);
}

