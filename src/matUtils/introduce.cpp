#include "introduce.hpp"
#include "select.hpp"

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
     "Set to calculate additional phylogenetic trait association statistics for whole regions and individual introductions. WARNING: Adds significantly to runtime.")
    ("clade-regions,c", po::value<std::string>()->default_value(""),
     "Set to optionally record, for each clade root in the tree, the support for that clade root being IN each region in the input, as a tsv with the indicated name.")
    ("date-metadata,M", po::value<std::string>()->default_value(""),
     "Pass a TSV or CSV containing a 'date' column to use for date information. If not used, date will be inferred from the sample name where possible.")
    ("full-output,o", po::value<std::string>()->required(),
     "Name of the file to save by-sample introduction information to.")
    ("origin-confidence,C", po::value<float>()->default_value(0.5),
     "Set the threshold for recording of putative origins of introductions. Default is 0.5")
    ("evaluate-metadata,E", po::bool_switch(),
     "Set to assign each leaf a confidence value based on ancestor distance and confidence.")
    ("dump-assignments,D", po::value<std::string>()->default_value(""),
     "Indicate a directory to which two-column text files containing node assignment values should be dumped for downstream processing.")
    ("latest-date,l", po::value<std::string>()->default_value("1500/1/1"),
     "Use to filter to clusters which have samples after the indicated date.")
    ("cluster-output,u", po::value<std::string>()->default_value(""),
    "Write by-cluster summary output to the indicated file.")
    // ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(filt_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
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

float get_association_index(MAT::Tree* T, std::map<std::string, float> assignments, bool permute, MAT::Node* subroot) {
    /*
    The association index was introduced by Wang et al 2005 for the estimation of phylogeny and trait correlation. Parker et al 2008 has a good summary.
    It's an index that is small for strong correlation and large for weak correlation, with non-integer values.
    AI = sum(for all internal nodes) (1-tips_with_trait) / (2 ^ (total_tips - 1))
    This can be calculated for a full tree or for an introduction-specific subtree.

    My implementation is an efficient one that relies on dynamic programming and useful ordering of node traversal.
    This searches over a reverse breadth first order. For each internal node, check the children.
    Count the number of direct leaf children which are in/out and get the records for the counts for in/out
    from the internal_tracker map. Add these up to get the total number of in/out leaves associated with a node.
    This avoids repeatedly traversing the tree to identify leaves for each internal node, since reverse breadth-first search order
    means that all children of a node will necessarily be traversed over before the node itself is.

    The permute boolean, when true, sets a mode where in/out traits are randomly assigned on a uniform distribution based on
    the baseline frequency of the trait across the tree instead of actually checking membership. Used to evaluate the results with
    some statistical grounding, though full repeated permutations are too computationally costly.
    */
    //timer.Start();
    float total_ai = 0.0;
    std::map<std::string,std::pair<size_t,size_t>> internal_tracker;
    std::vector<MAT::Node*> bfs;
    if (subroot != NULL) {
        bfs = T->breadth_first_expansion(subroot->identifier);
    } else {
        bfs = T->breadth_first_expansion();
    }

    srand(time(nullptr));
    size_t leaf_count = 0;
    size_t permuted_inc = 0;
    size_t sample_count = 0;
    if (permute) {
        for (auto b: bfs) {
            if (b->is_leaf()) {
                leaf_count++;
                auto search = assignments.find(b->identifier);
                if (search != assignments.end()) {
                    if (search->second > 0.5) {
                        sample_count++;
                    }
                }
            }
        }
    }

    std::reverse(bfs.begin(), bfs.end());
    for (auto n: bfs) {
        if (!n->is_leaf()) {
            size_t in_c = 0;
            size_t out_c = 0;
            for (auto c: n->children) {
                if (c->is_leaf()) {
                    if (permute) {
                        int random = rand()%leaf_count;
                        if (random <= static_cast<float>(sample_count)) {
                            in_c++;
                            permuted_inc++;
                        } else {
                            out_c++;
                        }
                    } else {
                        auto search = assignments.find(c->identifier);
                        if (search != assignments.end()) {
                            if (search->second > 0.5) {
                                in_c++;
                            } else {
                                out_c++;
                            }
                        }
                    }
                } else {
                    //check to see if we've recorded the internal child- we definitely should have in reverse BFS
                    auto search = internal_tracker.find(c->identifier);
                    if (search != internal_tracker.end()) {
                        in_c += search->second.first;
                        out_c += search->second.second;
                    } else {
                        //this SHOULD not happen, logically...
                        fprintf(stderr, "ERROR: AI calculation encountered mystery internal child node\n");
                        exit(1);
                    }
                }
            }
            internal_tracker[n->identifier] = std::make_pair(in_c,out_c);
            size_t total_leaves = in_c + out_c;
            total_ai += ((1 - std::max(in_c,out_c)/total_leaves) / (pow(2, (total_leaves-1))));
        }
    }
    return total_ai;
}

size_t get_monophyletic_cladesize(MAT::Tree* T, std::map<std::string, float> assignments, MAT::Node* subroot) {
    /*
    The monophyletic clade statistic was introduced by Salemi et al 2005. Parker et al 2008 has a good summary.
    MC is bigger for strong correlations, bounded 1 to N where N is the number of samples in the subtree.
    This is a simple qualifier which just searches across the subtree and identifies the largest clade which entirely and only contains IN samples.

    My implementation accomplishes this by looking across the depth-first expansion of leaves only
    and identifying the longest contiguous string of IN sample identifiers, which reflects the largest
    clade subtree which is entirely IN.
    */
    size_t biggest = 0;
    std::vector<std::string> acls;
    //depth-first search order is required for this implementation.
    for (auto n: T->depth_first_expansion(subroot)) {
        if (n->is_leaf()) {
            acls.push_back(n->identifier);
        }
    }
    size_t current = 0;
    for (auto l: acls) {
        auto search = assignments.find(l);
        if (search != assignments.end()) {
            if (search->second >= 0.5) {
                current++;
            } else {
                if (current > biggest) {
                    biggest = current;
                }
                current = 0;
            }
        }
    }
    if (current > biggest) {
        biggest = current;
    }
    return biggest;
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
        for (auto ca: n->clade_annotations) {
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

std::map<std::string, float> get_assignments(MAT::Tree* T, std::unordered_set<std::string> sample_set, bool eval_uncertainty) {
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
                    if ((min_to_in > 0) & (min_to_out > 0)) {
                        //if both of these are assigned, we're done. break out of this loop and move on
                        break;
                    }
                    if (d->is_leaf()) {
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
                    float r = (vir/vor);
                    float c = (1/(1+r));
                    if (isnan(c)) {
                        fprintf(stderr, "ERROR: Invalid introduction assignment calculation. Debug information follows.\n");
                        fprintf(stderr, "min %ld, mout %ld, ols %ld, ils %ld,", min_to_in, min_to_out, out_leaves.size(), in_leaves.size());
                        fprintf(stderr, " vor %f, vir %f, r %f\n", vor, vir, r);
                        exit(1);
                    }
                    assignments[n->identifier] = c;
                }
            }
        }
    }
    if (eval_uncertainty) {
        timer.Start();
        fprintf(stderr, "Leaf label uncertainty estimate requested; calculating...\n");
        //update the assignments for specific leaves against the rest of the dataset
        for (auto l: T->get_leaves()) {
            float total_conf = 0.0;
            float max_conf = 0.0;
            float traversed = static_cast<float>(l->mutations.size());
            for (auto anc: T->rsearch(l->identifier, false)) {
                float acv = assignments.find(anc->identifier)->second;
                total_conf += (acv / ((1+traversed) * (1+traversed)));
                max_conf += (1 / ((1+traversed) * (1+traversed)));
                traversed += static_cast<float>(anc->mutations.size());
            }
            assignments[l->identifier] = total_conf / max_conf;
        }
        fprintf(stderr, "All leaves processed in %ld msec.\n", timer.Stop());
    }
    return assignments;
}

std::pair<boost::gregorian::date,boost::gregorian::date> get_nearest_date(MAT::Tree* T, MAT::Node* n, std::set<std::string>* in_samples, std::map<std::string, std::string> datemeta) {
    boost::gregorian::date earliest = boost::gregorian::day_clock::universal_day();
    boost::gregorian::date latest = boost::gregorian::date(1500,1,1);
    for (auto l: T->get_leaves_ids(n->identifier)) {
        if (in_samples->find(l) != in_samples->end()) {
            if (datemeta.size() > 0) {
                if (datemeta.find(l) != datemeta.end()) {
                    boost::gregorian::date leafdate;
                    try {
                        leafdate = boost::gregorian::from_string(datemeta.find(l)->second);
                    } catch (boost::bad_lexical_cast &e) {
                        fprintf(stderr, "WARNING: Malformed date %s provided in date file for sample %s; ignoring sample date\n", datemeta.find(l)->second.c_str(), l.c_str());
                        continue;
                    }
                    if (leafdate < earliest) {
                        earliest = leafdate;
                    }
                    if (leafdate > latest) {
                        latest = leafdate;
                    }
                }
            } else {
                std::string datend = l.substr(l.rfind("|")+1, std::string::npos);
                if (datend.size() > 0) {
                    boost::gregorian::date leafdate;
                    try {
                        if (datend.size() == 8) {
                            leafdate = boost::gregorian::from_string("20"+datend);
                        } else if (datend.size() == 10) {
                            leafdate = boost::gregorian::from_string(datend);
                        } else {
                            continue;
                        }
                    } catch (const std::out_of_range& oor) {
                        continue;
                    }
                    if (leafdate < earliest) {
                        earliest = leafdate;
                    }
                    if (leafdate > latest) {
                        latest = leafdate;
                    }
                }
            }
        }
    }
    if ((earliest == boost::gregorian::day_clock::universal_day()) && (latest == boost::gregorian::date(1500,1,1))) {
        return std::pair<boost::gregorian::date,boost::gregorian::date> ();
    }
    return std::pair<boost::gregorian::date,boost::gregorian::date> (earliest,latest);
}

std::vector<std::string> find_introductions(MAT::Tree* T, std::map<std::string, std::vector<std::string>> sample_regions, bool add_info, std::string clade_output, float min_origin_confidence, std::string bycluster, std::string dump_assignments, bool eval_uncertainty, std::string latest_date = "1700/1/1", std::map<std::string, std::string> datemeta = {}) {
    //for every region, independently assign IN/OUT states
    //and save these assignments into a map of maps
    //so we can check membership of introduction points in each of the other groups
    //this allows us to look for migrant flow between regions
    std::map<std::string, std::map<std::string, float>> region_assignments;
    boost::gregorian::date recency_filter;
    std::vector<std::string> bycluster_output;
    try {
        recency_filter = boost::gregorian::from_string(latest_date);
    } catch (const std::out_of_range& oor) {
        fprintf(stderr, "ERROR: Minimum date argument (-l) could not be parsed. Check that it is formatted year-month-day and try again.\n");
        exit(1);
    }
    //TODO: This could be parallel for a significant speedup when dozens or hundreds of regions are being passed in
    //I also suspect I could use pointers for the assignment maps to make this more memory efficient
    for (auto ms: sample_regions) {
        std::string region = ms.first;
        std::vector<std::string> samples = ms.second;
        fprintf(stderr, "Processing region %s with %ld total samples\n", region.c_str(), samples.size());
        std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
        auto assignments = get_assignments(T, sample_set, eval_uncertainty);
        if (add_info) {
            size_t global_mc = get_monophyletic_cladesize(T, assignments);
            float global_ai = get_association_index(T, assignments);
            fprintf(stderr, "Region largest monophyletic clade: %ld, regional association index: %f\n", global_mc, global_ai);
            std::vector<float> permvec;
            for (int i = 0; i < 100; i++) {
                float perm = get_association_index(T, assignments, true);
                permvec.push_back(perm);
            }
            std::sort(permvec.begin(), permvec.end());
            fprintf(stderr, "Real value %f. Quantiles of random expected AI for this sample size: %f, %f, %f, %f, %f\n", global_ai, permvec[5],permvec[25],permvec[50],permvec[75],permvec[95]);
        }

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
            if (ass.second > min_origin_confidence) {
                region_ins[ass.first].push_back(ra.first);
                region_cons[ass.first].push_back(ass.second);
            }
        }
    }
    fprintf(stderr, "Regions processed; identifying introductions.\n");
    //now that we have all assignments sorted out, pass over it again
    //looking for introductions.
    std::vector<std::string> outstrs;
    std::string header = "sample\tintroduction_node\tintroduction_rank\tgrowth_score\tearliest_date\tlatest_date\tcluster_size\tcluster_span\tintro_confidence\tparent_confidence\tdistance\torigin_gap";
    if (region_assignments.size() > 1) {
        header += "\tregion\torigins\torigins_confidence";
    }
    header += "\tclades\tmutation_path";
    if (eval_uncertainty) {
        header += "\tmeta_uncertainty";
    }
    if (add_info) {
        header += "\tmonophyl_size\tassoc_index\n";
    } else {
        header += "\n";
    }
    outstrs.emplace_back(header);
    for (auto ra: region_assignments) {
        std::string region = ra.first;
        auto assignments = ra.second;
        std::vector<std::string> samples = sample_regions[region];
        std::set<std::string> sampleset (samples.begin(), samples.end());
        std::map<std::string, size_t> recorded_mc;
        std::map<std::string, float> recorded_ai;
        std::map<std::string,std::map<std::string,std::string>> clusters;
        std::map<std::string,std::string> clustermeta;
        //its time to identify introductions. we're going to do this by iterating
        //over all of the leaves. For each sample, rsearch back until it hits a 0 assignment
        //then record the last encountered 1 assignment as the point of introduction
        int total_processed = 0;
        for (auto s: samples) {
            //timer.Start();
            //total_processed++;
            //everything in this vector is going to be considered 1 (IN) this region
            std::string last_encountered = s;
            MAT::Node* last_node = NULL;
            float last_anc_state = 1;
            auto node = T->get_node(s);
            if (node == NULL) {
                fprintf(stderr, "WARNING: query sample %s not found in tree. continuing\n", s.c_str());
                continue;
            }
            size_t traversed = node->mutations.size();
            for (auto a: T->rsearch(s,false)) {
                float anc_state;
                if (a->is_root()) {
                    //if we get back to the root, the root is necessarily the point of introduction for this sample
                    last_encountered = a->identifier;
                    anc_state = 0;
                } else {
                    //every node should be in assignments at this point.
                    // fprintf(stderr, "DEBUG: checking ancestor %s\n", a->identifier.c_str());
                    anc_state = assignments.find(a->identifier)->second;
                }
                if (anc_state < min_origin_confidence) {
                    //check whether this 0 node is 1 in any other region
                    //record each region where this is true
                    //(in the single region case, its never true, but its only like two operations to check anyways)
                    //revise this to assign origin as the highest confidence IN for ANY region, no matter how lacking in confidence it is
                    std::string origins;
                    std::stringstream origins_cons;
                    //can't assign region of origin if introduction point is root (no information about parent)
                    if ((region_assignments.size() > 1) & (!a->is_root())) {
                        auto assign_search = region_ins.find(a->identifier);
                        // float highest_conf = 0.0;
                        // std::string highest_conf_origin = "indeterminate";
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
                    //get the mutations and clade information
                    //both of these require an rsearch back from the point of origin to the tree root
                    //mutation path is going to be in reverse for simplicity, so using < to indicate direction
                    for (MAT::Node* as: T->rsearch(a->identifier, true)) {
                        //collect mutations
                        std::string mutstr;
                        for (auto m: as->mutations) {
                            if (mutstr.size() == 0) {
                                mutstr += m.get_string();
                            } else {
                                mutstr += "," + m.get_string();
                            }
                        }
                        intro_mut_path += mutstr + "<";
                        //check for any clade identifiers, record comma delineated
                        for (auto ann: as->clade_annotations) {
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
                    //calculate the trait phylogeny association metrics
                    //only if additional info is requested because of how it affects total runtime.
                    size_t mc = 0;
                    float ai = 0.0;
                    if (add_info) {
                        auto mc_s = recorded_mc.find(a->identifier);
                        if (mc_s != recorded_mc.end()) {
                            mc = mc_s->second;
                        } else {
                            mc = get_monophyletic_cladesize(T, assignments, last_node);
                            recorded_mc[a->identifier] = mc;
                        }
                        auto ai_s = recorded_ai.find(a->identifier);
                        if (ai_s != recorded_ai.end()) {
                            ai = ai_s->second;
                        } else {
                            ai = get_association_index(T, assignments, false, last_node);
                            recorded_ai[a->identifier] = ai;
                        }
                    }
                    //ostr stores by-sample and cluster-specific information both together for the full output
                    std::stringstream ostr;
                    //mcl stores just cluster-specific information for a separate output file.
                    std::stringstream mcl;
                    if (region_assignments.size() == 1) {
                        ostr << "\t" << last_anc_state << "\t" << anc_state << "\t" << traversed << "\t" << a->mutations.size() << "\t" << intro_clades << "\t" << intro_mut_path;
                        mcl << last_anc_state << "\t" << anc_state << "\t" << a->mutations.size() << "\t" << intro_clades << "\t" << intro_mut_path;
                        if (eval_uncertainty) {
                            ostr << "\t" << assignments.find(s)->second;
                        }
                        if (add_info) {
                            ostr << "\t" << mc << "\t" << ai << "\n";
                            mcl << "\t" << mc << "\t" << ai << "\n";
                        } else {
                            ostr << "\n";
                            mcl << "\n";
                        }
                    } else {
                        ostr << "\t" << last_anc_state << "\t" << anc_state << "\t" << traversed << "\t" << a->mutations.size() << "\t" << region << "\t" << origins << "\t" << origins_cons.str() << "\t" << intro_clades << "\t" << intro_mut_path;
                        mcl << last_anc_state << "\t" << anc_state << "\t" << a->mutations.size() << "\t" << region << "\t" << origins << "\t" << origins_cons.str() << "\t" << intro_clades << "\t" << intro_mut_path;
                        if (eval_uncertainty) {
                            ostr << "\t" << assignments.find(s)->second;
                        }
                        if (add_info) {
                            ostr << "\t" << mc << "\t" << ai << "\n";
                            mcl << "\t" << mc << "\t" << ai << "\n";
                        } else {
                            ostr << "\n";
                            mcl << "\n";
                        }
                    }
                    clusters[last_encountered][s] = ostr.str();
                    //keying the meta for clusters on each sample is redundant and gross, but I can't think of a better way to pass information from 
                    clustermeta[last_encountered] = mcl.str();
                    total_processed++;
                    break;
                } else {
                    last_encountered = a->identifier;
                    last_node = a;
                    last_anc_state = anc_state;
                    traversed += a->mutations.size();
                }
            }
        }
        std::vector<float> growthv;
        std::map<float,std::vector<std::string>> cgm;
        std::map<std::string, std::string> date_tracker;
        for (auto cs: clusters) {
            std::string ldatestr;
            MAT::Node* nn = T->get_node(cs.first);
            std::pair<boost::gregorian::date,boost::gregorian::date> ldates = get_nearest_date(T, nn, &sampleset, datemeta);
            boost::gregorian::days diff(0);
            if ((ldates.first.is_not_a_date()) || (ldates.second.is_not_a_date())) {
                fprintf(stderr, "WARNING: Cluster %s has no valid dates included among samples\n", cs.first.c_str());
                ldatestr = "no-valid-date\tno-valid-date";
            } else {
                if (recency_filter > ldates.second) {
                    continue;
                }
                ldatestr = boost::gregorian::to_simple_string(ldates.first) + "\t" + boost::gregorian::to_simple_string(ldates.second);
                diff = (ldates.second - ldates.first);
            }
            date_tracker[cs.first] = ldatestr;
            float gv;
            gv = static_cast<float>(cs.second.size()) / static_cast<float>((int)(diff.days()/7)+1);
            growthv.emplace_back(gv);
            cgm[gv].emplace_back(cs.first);
        }
        //sort by default goes from smallest to largest
        //I want to rank by largest to smallest, so this is reversed.
        //tiebreaker ordering has to do with the order of samples encountered.
        std::sort(growthv.begin(), growthv.end());
        std::reverse(growthv.begin(), growthv.end());
        auto rm = std::unique(growthv.begin(), growthv.end());
        growthv.erase(rm, growthv.end());
        size_t rankr = 0;
        for (size_t i = 0; i < growthv.size(); i++) {
            float gv = growthv[i];
            for (auto cid: cgm[gv]) {
                //calculate the total span of branch lengths across this cluster.
                std::vector<std::string> cs;
                for (auto ss: clusters[cid]) {
                    cs.push_back(ss.first);
                }
                size_t span = 0;
                if (cs.size() > 1) {
                    std::set<std::string> ancm;
                    for (auto s: cs) {
                        for (auto a: T->rsearch(s, true)) {
                            if (a->identifier == cid) {
                                break;
                            } else if (ancm.find(a->identifier) == ancm.end()) {
                                span += a->mutations.size();
                                ancm.insert(a->identifier);
                            } else {
                                break;
                            }
                        }
                    }
                } else {
                    span = T->get_node(cs[0])->mutations.size();
                }
                //yes, I'm iterating over this multiple times. Nothing is ever easy.
                //new- store a single line describing cluster-unique attributes and counts. 
                std::stringstream clo;
                clo << cid << "\t" << clusters[cid].size() << "\t" << date_tracker[cid] << "\t" << gv << "\t" << span << "\t" << clustermeta[cid];
                bycluster_output.push_back(clo.str());
                rankr++;
                for (auto ss: clusters[cid]) {
                    std::stringstream cout;
                    //in order, first seven columns are
                    //sample id, cluster id, cluster rank, cluster growth score, earliest date, latest date, cluster size
                    //then the rest are the by-sample information (path, distance of this specific sample, yadda yadda)
                    if (date_tracker.find(cid) == date_tracker.end()) {
                        continue;
                    }
                    cout << ss.first << "\t" << cid << "\t" << rankr << "\t" << gv << "\t" << date_tracker[cid] << "\t" << clusters[cid].size() << "\t" << span << ss.second;
                    outstrs.push_back(cout.str());
                }
            }
        }
        fprintf(stderr, "Region %s complete, %d samples processed.\n", region.c_str(), total_processed);
    }
    if (dump_assignments != "") {
        boost::filesystem::path path(dump_assignments);
        if (!boost::filesystem::exists(path)) {
            fprintf(stderr, "Creating output directory to dump region assignments.\n\n");
            boost::filesystem::create_directory(dump_assignments);
        }
        for (auto ra: region_assignments) {
            std::ofstream rof(dump_assignments + "/" + ra.first + "_assignments.tsv");
            rof << "sample\tconfidence_continuous\n";
            for (auto ass: ra.second) {
                //only save nodes with non-zero confidence values for the sake of file size.
                if (ass.second > 0) {
                    rof << ass.first << "\t" << ass.second << "\n";
                }
            }
            rof.close();
        }
    }
    if (bycluster != "") {
        std::ofstream cof(bycluster);
        cof << "cluster_id\tsample_count\tearliest_date\tlatest_date\tgrowth_score\tspan\tintro_confidence\tparent_confidence\torigin_gap\t";
        if (add_info) {
            cof << "monophyletic_cladesize\tassociation_index\t";
        }
        if (region_assignments.size() == 1) {
            cof << "clade\tmutation_path\n";
        } else {
            cof << "region\tinferred_origin\tinferred_origin_confidence\tclade\tmutation_path\n";
        }
        for (auto os: bycluster_output) {
            cof << os;
        }
        cof.close();
    }
    return outstrs;
}

void introduce_main(po::parsed_options parsed) {
    po::variables_map vm = parse_introduce_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string samples_filename = vm["population-samples"].as<std::string>();
    std::string clade_regions = vm["clade-regions"].as<std::string>();
    std::string metafile = vm["date-metadata"].as<std::string>();
    std::string latest_date = vm["latest-date"].as<std::string>();
    bool add_info = vm["additional-info"].as<bool>();
    std::string output_file = vm["full-output"].as<std::string>();
    std::string dump_assignments = vm["dump-assignments"].as<std::string>();
    float moconf = vm["origin-confidence"].as<float>();
    bool leafconf = vm["evaluate-metadata"].as<bool>();
    std::string clusterout = vm["cluster-output"].as<std::string>();
    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    //T here is the actual object.
    if (T.condensed_nodes.size() > 0) {
        T.uncondense_leaves();
    }
    auto region_map = read_two_column(samples_filename);
    std::map<std::string,std::string> datemeta = {};
    if (metafile.size() > 0) {
        std::set<std::string> samples;
        for (auto kv: region_map) {
            samples.insert(kv.second.begin(), kv.second.end());
        }
        auto metad = read_metafile(metafile, samples);
        if (metad.find("date") != metad.end()) {
            datemeta = metad["date"];
        } else {
            fprintf(stderr, "ERROR: Metadata file does not contain required column 'date'; exiting\n");
            exit(1);
        }
    }
    auto outstrings = find_introductions(&T, region_map, add_info, clade_regions, moconf, clusterout, dump_assignments, leafconf, latest_date, datemeta);

    std::ofstream of;
    of.open(output_file);
    for (std::string o: outstrings) {
        of << o;
    }
    of.close();
}
