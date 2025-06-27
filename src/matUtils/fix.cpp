#include "fix.hpp"

// Default: move node if it has at least 1 descendent
int default_min_descendent_count = 1;
int default_iterations = 1;

po::variables_map parse_fix_command(po::parsed_options parsed) {

    po::variables_map vm;
    po::options_description filt_desc("fix options");
    filt_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-mat,o", po::value<std::string>()->required(),
     "Path to output fixed mutation-annotated tree file [REQUIRED]")
    ("iterations,n", po::value<int>()->default_value(default_iterations),
     "Number of iterations to run")
    ("min-descendent-count,c", po::value<int>()->default_value(default_min_descendent_count),
     "Minimum number of descendents required to move a node")
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

static int fix_grandparent_reversions_r(MAT::Tree *T, MAT::Node *node, MAT::Node *ggp_node,
                                        MAT::Node *gp_node, MAT::Node *p_node,
                                        int min_descendent_count) {
    // Recursively scan the tree for cases of grandparent-reversion, i.e. N > A > B > revA;
    // when a case like that is found, move the revA node to be a child of N having mutation B.
    int descendent_count = 0;
    // First recursively descend to children, looking each one up by identifier in case it has
    // been removed by the time we get to it:
    std::vector<std::string> child_ids;
    for (auto child: node->children) {
        child_ids.push_back(child->identifier);
    }
    for (auto child_id: child_ids) {
        MAT::Node *child = T->get_node(child_id);
        if (child != NULL) {
            descendent_count += fix_grandparent_reversions_r(T, child, gp_node, p_node, node,
                                                             min_descendent_count);
        }
    }
    // Now determine whether this node has only one mutation that reverts its grandparent's only
    // mutation; if so (and if parent has only one mut, othw parsimony score would increase),
    // then move this node to be a child of its great-grandparent, with only its parent's mut.
    if (ggp_node != NULL && node->mutations.size() == 1 && gp_node->mutations.size() == 1 &&
        p_node->mutations.size() == 1) {
        MAT::Mutation node_mut = node->mutations.front();
        MAT::Mutation gp_mut = gp_node->mutations.front();
        if (node_mut.position == gp_mut.position &&
            node_mut.chrom == gp_mut.chrom &&
            node_mut.mut_nuc == gp_mut.par_nuc &&
            node_mut.par_nuc == gp_mut.mut_nuc &&
            descendent_count >= min_descendent_count) {
            fprintf(stderr, "Node %s mutation %s reverts grandparent %s's %s%s, moving to %s with %s (%d descendents)\n",
                    node->identifier.c_str(), node_mut.get_string().c_str(),
                    gp_node->identifier.c_str(),
                    (gp_mut.mut_nuc == gp_mut.ref_nuc ? "reversion " : ""),
                    gp_mut.get_string().c_str(),
                    ggp_node->identifier.c_str(),
                    p_node->mutations.front().get_string().c_str(), descendent_count);
            node->mutations.clear();
            node->mutations = std::vector<MAT::Mutation>(p_node->mutations);
            T->move_node(node->identifier, ggp_node->identifier);
        }
    }
    return descendent_count + 1;
}

static void fix_grandparent_reversions(MAT::Tree *T, int iterations, int min_descendent_count) {
    // Find and fix nodes that reverse their grandparent's mutation.
    // For each case of N > A > B > revA, move revA to N > B (if N > B already exists then move all
    // children of revA to N > B).
    int i;
    for (i = 0;  i < iterations;  i++) {
        fix_grandparent_reversions_r(T, T->root, NULL, NULL, NULL, min_descendent_count);
    }
}

void fix_main(po::parsed_options parsed) {
    po::variables_map vm = parse_fix_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string output_mat_filename = vm["output-mat"].as<std::string>();
    int iterations = vm["iterations"].as<int>();
    int min_descendent_count = vm["min-descendent-count"].as<int>();

    // Load the input MAT
    timer.Start();
    fprintf(stderr, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    // No need to uncondense the tree

    timer.Start();
    fprintf(stderr, "Fixing grandparent-reversions\n");
    fix_grandparent_reversions(&T, iterations, min_descendent_count);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    timer.Start();
    fprintf(stderr, "Saving Final Tree\n");
    MAT::save_mutation_annotated_tree(T, output_mat_filename);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}
