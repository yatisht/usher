#include "introduce.hpp"

po::variables_map parse_introduce_command(po::parsed_options parsed) {

    // uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    // std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";

    po::variables_map vm;
    po::options_description filt_desc("introduce options");
    filt_desc.add_options()
        ("input-mat,i", po::value<std::string>()->required(),
         "Input mutation-annotated tree file [REQUIRED]")
        ("population-samples,s", po::value<std::string>()->required(), 
         "Names of samples from the population of interest [REQUIRED].") 
        ("additional-info,a", po::bool_switch(),
        "Set to print clade membership and full mutation paths at inferred introduction points.")
        ("clade-regions,c", po::value<std::string>()->default_value(""),
        "Set to optionally record, for each clade root in the tree, the support for that clade root being IN each region in the input, as a tsv with the indicated name.")
        ("output,o", po::value<std::string>()->required(),
        "Name of the file to save the introduction information to.")
        // ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try{
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << filt_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

std::map<std::string, std::vector<std::string>> read_two_column (std::string sample_filename) {
    std::map<std::string, std::vector<std::string>> amap;
    //this variant on a sample reader can optionally parse a second column of values
    //because this is handling multiple vectors, it takes vector pointers instead of defining its own and returning them
    std::ifstream infile(sample_filename);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the file: %s!\n", sample_filename.c_str());
        exit(1);
    }    
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<std::string> words;
        MAT::string_split(line, words);

        std::string sname = words[0];
        std::string rname;
        if (words.size() == 1) {
            //default region for a single-column file is "default".
            rname = "default";
            //carriage return handling
            if (sname[sname.size()-1] == '\r') {
                sname = sname.substr(0,sname.size()-1);
            }
        } else if (words.size() == 2) {
            rname = words[1];
            if (rname[rname.size()-1] == '\r') {
                rname = rname.substr(0,rname.size()-1);
            }
        } else {
            fprintf(stderr, "ERROR: Too many columns in file %s- check format\n", sample_filename.c_str());
            exit(1);
        }
        amap[rname].push_back(sname);
    }
    infile.close();
    return amap;
}

void record_clade_regions(MAT::Tree* T, std::map<std::string, std::map<std::string, float>> region_assignments, std::string filename) {
    //record a tsv with a column for each annotated region (single will have label default)
    //and a row for each clade label. The contents are the assignment support for that specific clade root being IN the indicated region
    std::ofstream of;
    of.open(filename);
    //write the header line
    //save the regions into an explicit vector just to make sure we don't get scrambled (iteration over hash maps...)
    of << "clade\t";
    std::vector<std::string> region_order;
    for (auto r: region_assignments) {
        of << r.first << "\t";
        region_order.push_back(r.first);
    }
    of << "\n";

    for (auto n: T->depth_first_expansion()) {
        for (auto ca: n->clade_annotations){
            if (ca.size() > 0) {
                //we have a clade root here
                std::stringstream rowstr;
                rowstr << ca << "\t";
                for (auto r: region_order) {
                    auto search = region_assignments[r].find(n->identifier);
                    //this can be anything from 0 to 1.
                    rowstr << search->second << "\t";
                }
                of << rowstr.str() << "\n";
            }
        }
    }
}

std::map<std::string, float> get_assignments(MAT::Tree* T, std::unordered_set<std::string> sample_set) {
    /*
    This function applies a heuristic series of steps to label internal nodes as in or out of a geographic area
    based on their relationship to the samples in the input list. The rules are:
    1. If the node is a leaf, it is IN if it is in samples, otherwise OUT
    2. If all descendents of the node are IN, its IN
    3. If all descendents of the node are OUT, its OUT
    4. Nodes are assigned the state which yields the minimum distance to the nearest leaf of that type divided by the total number of leaves of that type descended from it
    This last step is the most complex one, and is a weighted minimum which takes into account both the number of descendents and the total distance to those descendents
    It's worth noting that if a sample is an identical child to the internal node, that state will always win
    5. On a tie, the node is assigned to the state of its parent

    Introductions are identified as locations where assignments in an rsearch from an IN sample shift from IN to OUT.
    */

    std::map<std::string, float> assignments;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        if (n->is_leaf()) {
            //rule 1
            if (sample_set.find(n->identifier) != sample_set.end()) {
                assignments[n->identifier] = 1; 
            } else {
                assignments[n->identifier] = 0;
            }
        } else {
            auto leaves = T->get_leaves_ids(n->identifier);
            //to apply rules 2-3, we need to check the state of each leaf
            std::unordered_set<std::string> in_leaves;
            std::unordered_set<std::string> out_leaves;
            for (auto l: leaves) {
                if (sample_set.find(l) != sample_set.end()) {
                    in_leaves.insert(l);
                } else {
                    out_leaves.insert(l);
                }
            }
            if (out_leaves.size() == 0) {
                //rule 2
                assignments[n->identifier] = 1;
            } else if (in_leaves.size() == 0) {
                //rule 3
                assignments[n->identifier] = 0;
            } else {
                //rule 4...
                //the best way to do this is to keep in mind that the nearest descendent leaf 
                //is going to be the next leaf encountered in DFS order
                //so we're going to iterate over DFS until we encounter a leaf of each type and record those distances
                size_t min_to_in = 0;
                size_t min_to_out = 0;
                for (auto d: T->depth_first_expansion(n)) {
                    if ((min_to_in > 0) & (min_to_out > 0)){
                        //if both of these are assigned, we're done. break out of this loop and move on
                        break;
                    }
                    if (d->is_leaf()){
                        if ((sample_set.find(d->identifier) != sample_set.end()) & (min_to_in == 0)) {
                            //found the nearest IN.
                            //rsearch back from it so that we're getting the direct path to the original node
                            //maybe there's a better way to do this, but I want to make sure I'm not including side branches
                            //on my iteration
                            size_t total_traveled = 0;
                            for (auto a: T->rsearch(d->identifier, true)) {
                                total_traveled += a->mutations.size();
                                if (a->identifier == n->identifier) {
                                    //back to the original query, break here
                                    min_to_in = total_traveled;
                                    break;
                                }
                            }
                        } else if ((sample_set.find(d->identifier) == sample_set.end()) & (min_to_out == 0)) {
                            //found the nearest out
                            //same routine as above
                            size_t total_traveled = 0;
                            for (auto a: T->rsearch(d->identifier, true)) {
                                total_traveled += a->mutations.size();
                                if (a->identifier == n->identifier) {
                                    //back to the original query, break here
                                    min_to_out = total_traveled;
                                    break;
                                }
                            }
                        } 
                    }
                }
                //update to rule 4- 1 is IN, 0 is OUT, but we want to have in-between numbers to represent relative confidence.
                //we calculate the balance by computing C=1/(1+((OUT_MD/OUT_LEAVES)/(IN_MD/IN_LEAVES)))
                //C is near 0 when OUT is large, C is near 1 when IN is large, C is 0.5 when they are the same
                //now we complete rule 4 by checking the balance.
                //I am getting lazier about tiebreaking with the new method but... oh well.
                if (min_to_in == 0) {
                    //this calculation is unnecessary in these cases.
                    //tiebreaker for both being 0 is IN with this ordering.
                    //identical IN sample, its IN.
                    assignments[n->identifier] = 1;
                } else if (min_to_out == 0) {
                    assignments[n->identifier] = 0;
                } else {
                    // fprintf(stderr, "DEBUG: min %ld, mout %ld, ols %ld, ils %ld\n", min_to_in, min_to_out, out_leaves.size(), in_leaves.size());                    
                    //unnecessary variable declarations because I was hitting floating point exceptions.
                    //float vor = (min_to_out/out_leaves.size());
                    //float vir = (min_to_in/in_leaves.size());
                    float vor = (static_cast<float>(min_to_out) / static_cast<float>(out_leaves.size()));
                    float vir = (static_cast<float>(min_to_in) / static_cast<float>(in_leaves.size()));
                    // fprintf(stderr, "DEBUG: or %f, ir %f\n", vor, vir);
                    float r = (vor/vir);
                    // fprintf(stderr, "DEBUG: ratio is %f, ", r);
                    float c = (1/(1+r));
                    // fprintf(stderr, "confidence is %f\n", c);
                    if (isnan(c)) {
                        fprintf(stderr, "DEBUG: nan out. min %ld, mout %ld, ols %ld, ils %ld,", min_to_in, min_to_out, out_leaves.size(), in_leaves.size());                    
                        fprintf(stderr, " vor %f, vir %f, r %f\n", vor, vir, r);
                    }
                    assignments[n->identifier] = c;
                }
            }
        }
    }
    return assignments;
}

std::vector<std::string> find_introductions(MAT::Tree* T, std::map<std::string, std::vector<std::string>> sample_regions, bool add_info, std::string clade_output) {
    //for every region, independently assign IN/OUT states
    //and save these assignments into a map of maps
    //so we can check membership of introduction points in each of the other groups
    //this allows us to look for migrant flow between regions
    std::map<std::string, std::map<std::string, float>> region_assignments;
    //TODO: This could be parallel for a significant speedup when dozens or hundreds of regions are being passed in
    //I also suspect I could use pointers for the assignment maps to make this more memory efficient
    for (auto ms: sample_regions) {
        std::string region = ms.first;
        std::vector<std::string> samples = ms.second;
        fprintf(stderr, "Processing region %s with %ld samples\n", region.c_str(), samples.size());
        std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
        auto assignments = get_assignments(T, sample_set);
        region_assignments[region] = assignments;
    }
    //if requested, record the clade output
    if (clade_output.size() > 0) {
        fprintf(stderr, "Clade root region support requested; recording...\n");
        record_clade_regions(T, region_assignments, clade_output);
    }
    //there are some significant runtime issues when checking assignment states repeatedly for different groups
    //so I'm going to generate a series of sets of nodes which are 1 for any given assignment
    //these will be used when looking for the origin of an introduction, as that only
    //cares about whether its 1 for a given region, not about the context.
    
    //this structure holds the ID of every node which is 1 in at least one region
    //and all corresponding regions it is 1 for. 
    std::map<std::string, std::vector<float>> region_cons;
    std::map<std::string, std::vector<std::string>> region_ins;
    for (auto ra: region_assignments) {
        for (auto ass: ra.second) {
            if (ass.second > 0.5) {
                region_ins[ass.first].push_back(ra.first);
                region_cons[ass.first].push_back(ass.second);
            }
        }
    }
    // fprintf(stderr, "DEBUG: %ld nodes are IN across all regions\n", region_ins.size());

    fprintf(stderr, "Regions processed; identifying introductions.\n");



    //now that we have all assignments sorted out, pass over it again
    //looking for introductions.
    std::vector<std::string> outstrs;
    if (region_assignments.size() == 1) {
        if (add_info) {
            outstrs.push_back("sample\tintroduction_node\tintro_confidence\tparent_confidence\tdistance\tclades\tmutation_path\n");
        } else {
            outstrs.push_back("sample\tintroduction_node\tintro_confidence\tparent_confidence\tdistance\n");
        }
    } else {
        //add a column for the putative origin of the sample introduction
        //if we're using multiple regions.
        if (add_info) {
            outstrs.push_back("sample\tintroduction_node\tintro_confidence\tparent_confidence\tdistance\tregion\torigins\torigins_confidence\tclades\tmutation_path\n");            
        } else {
            outstrs.push_back("sample\tintroduction_node\tintro_confidence\tparent_confidence\tdistance\tregion\torigins\torigins_confidence\n");
        }
    }
    for (auto ra: region_assignments) {
        std::string region = ra.first;
        auto assignments = ra.second;
        std::vector<std::string> samples = sample_regions[region];
        //its time to identify introductions. we're going to do this by iterating
        //over all of the leaves. For each sample, rsearch back until it hits a 0 assignment
        //then record the last encountered 1 assignment as the point of introduction
        //int total_processed = 0;
        for (auto s: samples) {
            //timer.Start();
            //total_processed++;
            //everything in this vector is going to be 1 (IN) this region
            std::string last_encountered = s;
            float last_anc_state = 1;
            size_t traversed = 0;
            for (auto a: T->rsearch(s,true)) {
                float anc_state;
                if (a->is_root()) {
                    //if we get back to the root, the root is necessarily the point of introduction for this sample
                    last_encountered = a->identifier;
                    anc_state = 0;
                } else {
                    //every node should be in assignments at this point.
                    anc_state = assignments.find(a->identifier)->second;
                }
                if (anc_state < 0.5) {
                    //check whether this 0 node is 1 in any other region
                    //record each region where this is true
                    //(in the single region case, its never true, but its only like two operations to check anyways)
                    std::string origins;
                    std::stringstream origins_cons;
                    //can't assign region of origin if introduction point is root (no information about parent)
                    if ((region_assignments.size() > 1) & (!a->is_root())) {
                        auto assign_search = region_ins.find(a->identifier);
                        if (assign_search != region_ins.end()) {
                            for (auto r: assign_search->second) {
                                if (origins.size() == 0) {
                                    origins += r;
                                } else {
                                    origins += "," + r;
                                }
                            }
                            for (auto v: region_cons.find(a->identifier)->second) {
                                origins_cons << v << ",";
                            }
                        } else {
                            origins = "indeterminate";
                            origins_cons << 0.0;
                        }
                    }
                    if (origins.size() == 0) {
                        //if we didn't find anything which has the pre-introduction node at 1, we don't know where it came from
                        origins = "indeterminate";
                        origins_cons << 0.0;
                    }
                    //collect additional information if requested.
                    std::string intro_clades = "";
                    std::string intro_mut_path = "";
                    if (add_info) {
                        //both of these require an rsearch back from the point of origin to the tree root
                        //mutation path is going to be in reverse for simplicity, so using < to indicate direction
                        for (auto a: T->rsearch(a->identifier, true)) {
                            //collect mutations
                            std::string mutstr;
                            for (auto m: a->mutations) {
                                if (mutstr.size() == 0) {
                                    mutstr += m.get_string();
                                } else {
                                    mutstr += "," + m.get_string();
                                }
                            }
                            intro_mut_path += mutstr + "<";
                            //check for any clade identifiers, record comma delineated
                            for (auto ann: a->clade_annotations) {
                                if (ann.size() > 0) {
                                    if (intro_clades.size() == 0) {
                                        intro_clades += ann;
                                    } else {
                                        intro_clades += "," + ann;
                                    }    
                                }
                            }
                        }
                        if (intro_clades.size() == 0) {
                            intro_clades = "none";
                        }
                    }
                    std::stringstream ostr;
                    if (region_assignments.size() == 1) {
                        ostr << s << "\t" << last_encountered << "\t" << last_anc_state << "\t" << anc_state << "\t" << traversed;
                        if (add_info) {
                            ostr << "\t" << intro_clades << "\t" << intro_mut_path << "\n";
                        } else {
                            ostr << "\n";
                        }
                    } else {
                        ostr << s << "\t" << last_encountered << "\t" << last_anc_state << "\t" << anc_state << "\t" << traversed << "\t" << region << "\t" << origins << "\t" << origins_cons.str();
                        if (add_info) {
                            ostr << "\t" << intro_clades << "\t" << intro_mut_path << "\n";
                        } else {
                            ostr << "\n";
                        }
                    }
                    outstrs.push_back(ostr.str());
                    break;
                } else {
                    last_encountered = a->identifier;
                    last_anc_state = anc_state;
                    traversed += a->mutations.size();
                }
            }
        //fprintf(stderr, "Found introduction for sample %s, time taken %ld msec, %d processed\n", s.c_str(), timer.Stop(), total_processed);
        }
    }
    return outstrs;
}

void introduce_main(po::parsed_options parsed) {
    po::variables_map vm = parse_introduce_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples_filename = vm["population-samples"].as<std::string>();
    std::string clade_regions = vm["clade-regions"].as<std::string>();
    bool add_info = vm["additional-info"].as<bool>();
    std::string output_file = vm["output"].as<std::string>();
    // int32_t num_threads = vm["threads"].as<uint32_t>();

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
      T.uncondense_leaves();
    }
    auto region_map = read_two_column(samples_filename);
    auto outstrings = find_introductions(&T, region_map, add_info, clade_regions);

    std::ofstream of;
    of.open(output_file);
    for (std::string o: outstrings) {
        of << o;
    }
    of.close();
}
