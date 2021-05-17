#include "convert.hpp"
#include "nlohmann_json.hpp"

using json = nlohmann::json;

void write_vcf_header(std::ostream& vcf_file, std::vector<Mutation_Annotated_Tree::Node*> &dfs,
                      bool print_genotypes) {
    // Write minimal VCF header with sample names in same order that genotypes
    // will be printed out (DFS).
    //fprintf(vcf_file, "##fileformat=VCFv4.2\n");
    vcf_file << "##fileformat=VCFv4.2\n";
    //fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if (print_genotypes) {
        //fprintf(vcf_file, "\tFORMAT");
        vcf_file << "\tFORMAT";
        for (auto node: dfs) {
            if (node->is_leaf()) {
                //fprintf(vcf_file, "\t%s", node->identifier.c_str());
                vcf_file << boost::format("\t%s") % node->identifier.c_str();
            }
        }
    }
    vcf_file << "\n";
    //fputc('\n', vcf_file);
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

void write_vcf_rows(std::ostream& vcf_file, MAT::Tree T, std::vector<MAT::Node*> &dfs, bool print_genotypes) {
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
              if (alts.size() == 0) {
                  fprintf(stderr, "WARNING: no-alternative site encountered in vcf output; skipping\n");
                  continue;
              }
              std::string id = make_id(ref, pos, alts);
              std::string alt_str = make_alt_str(alts);
              std::string info = make_info(alts, leaf_count);
              //fprintf(vcf_file, "%s\t%d\t%s\t%c\t%s\t.\t.\t%s",
                    //   chrom.c_str(), pos, id .c_str(), MAT::get_nuc(ref), alt_str.c_str(),
                    //   info.c_str());
              vcf_file << boost::format("%s\t%d\t%s\t%c\t%s\t.\t.\t%s")
                      % chrom.c_str() % pos % id .c_str() % MAT::get_nuc(ref) % alt_str.c_str() % info.c_str();
              if (print_genotypes) {
                  int *allele_codes = make_allele_codes(ref, alts);
                  //fprintf(vcf_file, "\tGT");
                  vcf_file << "\tGT";
                  for (uint i = 0;  i < leaf_count;  i++) {
                      int8_t allele = gt_array[i];
                      vcf_file << boost::format("\t%d") % allele_codes[allele];
                      //fprintf(vcf_file, "\t%d", allele_codes[allele]);
                  }
              }
              //fputc('\n', vcf_file);
              vcf_file << '\n';
            }
        }
    }
}

void make_vcf (MAT::Tree T, std::string vcf_filepath, bool no_genotypes) {
    try {
        std::ofstream outfile(vcf_filepath, std::ios::out | std::ios::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
        if (vcf_filepath.find(".gz\0") != std::string::npos) {
            outbuf.push(boost::iostreams::gzip_compressor());
        } 
        outbuf.push(outfile);
        std::ostream vcf_file(&outbuf);
        //FILE *vcf_file = fopen(vcf_filepath.c_str(), "w");
        std::vector<Mutation_Annotated_Tree::Node*> dfs = T.depth_first_expansion();
        write_vcf_header(vcf_file, dfs, !no_genotypes);
        write_vcf_rows(vcf_file, T, dfs, !no_genotypes);
        boost::iostreams::close(outbuf);
        outfile.close();
        //fclose(vcf_file);
    } catch (const boost::iostreams::gzip_error& e) {
        std::cout << e.what() << '\n';
    }
}


/// JSON functions below

std::string write_mutations(MAT::Node *N) { // writes muts as a list, e.g. "A23403G,G1440A,G23403A,G2891A" for "nuc mutations" under labels
    std::string muts = ": \"" ;
    for (unsigned int m = 0 ; m < N->mutations.size() ; m ++ ){
        auto mut_string = N->mutations[m].get_string();
        muts += mut_string ;
        if ( m < (N->mutations.size()-1) ){
            muts +=  "," ; // if not last, add comma
        }
    }
    return muts ;
}

//the below is more formal JSON parsing code. 
//the goal is to load a nextstrain JSON into a MAT structure
//which is compatible with various downstream tools.

void create_node_from_json(MAT::Tree* T, json nodeinfo, MAT::Node* parent = NULL, size_t counter = 0, size_t* warning_counter = 0) {
    //this function is recursive.
    //generate and save a node from the parent first
    //then for each child, run this function on them, with this node as the parent.
    if (nodeinfo.contains("branch_attrs")) {    
        std::string nid;
        if (nodeinfo.contains("name")) {
            nid = nodeinfo["name"];
        } else {
            nid = std::to_string(counter);
            counter++;
        }
        MAT::Node* n;
        if (parent != NULL) {
            n = T->create_node(nid, parent);    
        } else {
            n = T->create_node(nid, 0.0, 1);
        }
        auto battrs = nodeinfo.at("branch_attrs");
        if (battrs["mutations"].contains("nuc")) {
            std::vector<std::string> mutations;
            for (auto n: battrs["mutations"]["nuc"]) {
                mutations.push_back(n);
            }
            float blen = mutations.size();
            n->branch_length = blen;
            for (auto m: mutations) {
                MAT::Mutation mut;
                mut.chrom = "NC_045512"; //hardcoded for sars-cov-2, in line with the jsons.
                //the json encodes ambiguous bases as - sometimes, it seems.
                int8_t nucid;
                if (static_cast<char>(m[0]) == '-') {
                    //the MAT .pb format does NOT support ambiguous parent bases.
                    //skip these. 
                    //record the number of entries skipped for warning printouts.
                    (*warning_counter)++;
                    continue;
                } else {
                    nucid = MAT::get_nuc_id(m[0]);
                    mut.par_nuc = nucid; 
                    mut.ref_nuc = nucid; //JSON does not track the original reference vs the parent. We're going to treat the parent as reference.
                }
                mut.position = std::stoi(m.substr(1, m.size()-1));
                if (static_cast<char>(m[m.size()-1]) == '-') {
                    //skip these as well. causes issues with usher to have ambiguous bases on internal nodes.
                    (*warning_counter)++;
                    continue;
                } else {
                    mut.mut_nuc = MAT::get_nuc_id(m[m.size()-1]);
                }
                n->add_mutation(mut);
            }
        }
        if (battrs.contains("labels")) {
            if (battrs["labels"].contains("clade")) {
                n->clade_annotations.push_back(nodeinfo["branch_attrs"]["labels"]["clade"]);
            }
        }
        if (nodeinfo.contains("children")) {
            for (auto& cl: nodeinfo["children"].items()) {
                create_node_from_json(T, cl.value(), n, counter, warning_counter);
            }
        }
    } else {
        //there are sometimes "nodes" in the json
        //which do not represent actual samples, or anything, and do not have 
        //branch_attrs. in these cases, we want to continue with their children
        //treating the parent of this pseudo-node as their parent.
        if (nodeinfo.contains("children")) {
            for (auto& cl: nodeinfo["children"].items()) {
                create_node_from_json(T, cl.value(), parent, counter, warning_counter);
            }
        }
    }
}

MAT::Tree load_mat_from_json(std::string json_filename) {
    MAT::Tree T;
    std::ifstream json_in(json_filename);
    json j;
    json_in >> j;
    size_t wc = 0;
    create_node_from_json(&T, j["tree"], NULL, 0, &wc);
    if (wc > 0) {
        fprintf(stderr, "WARNING: %ld mutations are removed for ambiguity\n", wc);
    }
    return T;
}

json get_json_entry(MAT::Node* n, std::map<std::string,std::map<std::string,std::string>>* catmeta, size_t div = 0, bool use_clade_zero = false, bool use_clade_one = false) {
    //each node has 3 constituent attributes
    //node_attrs, branch_attrs, and children. If its a leaf,
    //it also has a simple name attribute.
    //branch_attrs contains mutation information.
    //node_attrs contains clade information.
    //children contains nested information about the child nodes.
    json sj;
    std::vector<std::string> mutids;
    std::string muts;
    for (auto m: n->mutations) {
        mutids.push_back(m.get_string());
        muts.append(m.get_string());
        if (m.get_string() != n->mutations.back().get_string()) {
            muts.append(",");
        }
    }
    std::map<std::string,std::vector<std::string>> nmap {{"nuc", mutids},};
    //mutation information is encoded twice in the output we're using. 
    //don't ask me why. 
    std::map<std::string,std::string> mutv {{"nuc mutations", muts}};
    sj["branch_attrs"] = {{"labels",mutv}, {"mutations",nmap}};
    //note: the below is pretty much sars-cov-2 specific. but so is all json-related things.
    //need to declare maps to get nlohmann to interpret these as key pairs.
    div += mutids.size();
    std::string c1av = "";
    std::string c2av = "";
    if (n->clade_annotations.size() >= 1) {
        c1av = n->clade_annotations[0];
    }
    if (n->clade_annotations.size() >= 2) {
        //json output supports two simultaneous clade annotations
        //being nextstrain and pangolin, since this is very specific
        //to sars-cov-2 phylogenomics. Additional fields are ignored
        //at this point.
        c1av = n->clade_annotations[1];
    }
    std::map<std::string,std::string> c1a {{"value",c1av}};
    std::map<std::string,std::string> c2a {{"value",c2av}};
    std::string country = n->identifier.substr(0, n->identifier.find("/"));
    std::string date = n->identifier.substr( n->identifier.find_last_of("|")+1, n->identifier.size() );
    std::map<std::string,std::string> com {{"value",country}};
    std::map<std::string,std::string> dam {{"value",date}};
    if ((n->is_leaf()) && (country.length() != n->identifier.size()) && (date.length() != n->identifier.size()) ) {
        sj["node_attrs"] = { {"country",com}, {"date",dam} ,{"div", div}, {"MAT_Clade_0", c1a}, {"MAT_Clade_1", c2a} };
    } else {
<<<<<<< HEAD
=======
        // sj["node_attrs"] = {{"div", div}, {"MAT_Clade_0", c1a}, {"MAT_Clade_1", c2a} };
>>>>>>> upstream/master
        sj["node_attrs"]["div"] = div;
        if (use_clade_zero) {
            sj["node_attrs"]["MAT_Clade_0"] = c1a;
        }
        if (use_clade_one) {
            sj["node_attrs"]["MAT_Clade_1"] = c2a;
        }
    }
    for (const auto& cmi: *catmeta) {
        if (cmi.second.find(n->identifier) != cmi.second.end()) {
            //store the metadata on both the branch and the tip for now. 
            sj["branch_attrs"]["labels"][cmi.first] = cmi.second.at(n->identifier);
            sj["node_attrs"][cmi.first]["value"] = cmi.second.at(n->identifier);
        }
    }
    sj["name"] = n->identifier;
    std::vector<json> child_json;
    for (auto cn: n->children) {
        json cj = get_json_entry(cn, catmeta, div, use_clade_zero, use_clade_one);
        child_json.push_back(cj);
        sj["children"] = child_json;
    }
    return sj;
}

void write_json_from_mat(MAT::Tree* T, std::string output_filename, std::map<std::string,std::map<std::string,std::string>>* catmeta) {
    json nj;
    std::string desc = "JSON generated by matUtils. If you have metadata you wish to display, you can now drag on a CSV/TSV file and it will be added into this view, [see here](https://docs.nextstrain.org/projects/auspice/en/latest/advanced-functionality/drag-drop-csv-tsv.html) for more info.";
    std::map<std::string,std::string> lm = {{"branch_label", "none"}};
    nj = {
        {"version","v2"},
        {"meta", {
            {"title","mutation_annotated_tree"},
            {"filters",json::array({"country"})},
            {"panels",json::array({"tree"})},
            {"colorings",{ {{"key","country"},{"title","Country"},{"type","categorical"}} }},
            {"display_defaults",lm},
            {"description",desc}
        }},
        {"tree",{
            {"name","wrapper"},
        }}
    };
    //add metadata to the header colorings if any exist
    if (catmeta->size()>0) {
        for (const auto& cmi: *catmeta) {
            std::map<std::string,std::string> mmap {{"key",cmi.first},{"title",cmi.first},{"type","categorical"}};
            nj["meta"]["colorings"].push_back(mmap);
        }
    }
    //check whether each of the mat clade annotation fields are used by any sample. 
    bool uses_clade_0 = false;
    bool uses_clade_1 = false;
    for (auto n: T->depth_first_expansion()) {
        if (n->clade_annotations.size() >= 1) {
            if (n->clade_annotations[0] != "") {
                uses_clade_0 = true;
            }
            if (n->clade_annotations.size() >= 2) {
                if (n->clade_annotations[1] != "") {
                    uses_clade_1 = true;
                }
            }
        }
        if ((uses_clade_0) && (uses_clade_1)) {
            break;
        }
    }
    if (uses_clade_0) {
        std::map<std::string,std::string> c1map {{"key","MAT_Clade_0"},{"title","MAT_Clade_1"},{"type","categorical"}};
        nj["meta"]["colorings"].push_back(c1map);
    }
    if (uses_clade_1) {
        std::map<std::string,std::string> c2map {{"key","MAT_Clade_1"},{"title","MAT_Clade_2"},{"type","categorical"}};
        nj["meta"]["colorings"].push_back(c2map);
    }
    auto treestuff = get_json_entry(T->root, catmeta, 0, uses_clade_0, uses_clade_1);
    nj["tree"]["children"] = json::array({treestuff});
    std::ofstream out(output_filename);
    // out << std::setw(4) << nj << std::endl;
    out << nj << std::endl;
    out.close();
}
