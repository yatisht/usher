#include "filter.hpp"

/*
Functions in this module take a MAT and a set of samples and return the a smaller MAT
representing those samples in an optimal manner.
*/
MAT::Tree filter_master(const MAT::Tree& T, std::vector<std::string> sample_names, bool prune) {
    MAT::Tree subtree;
    if (prune) {
        subtree = prune_leaves(T, sample_names);
    } else if (sample_names.size() < 10000) {
        //for retaining only a subtree, get_subtree is the most effective method
        //for subtress up to about 10 thousand samples in size; after that, pruning
        //all other nodes becomes more efficient, because get_subtree scales with 
        //the size of the input sample set, while prune takes a similar time for 
        //any sample set size, while scaling on total tree size
        subtree = get_sample_subtree(T, sample_names);
    } else {
        subtree = get_sample_prune(T, sample_names);
    }
    return subtree;
}

MAT::Tree prune_leaves (const MAT::Tree& T, std::vector<std::string> sample_names) {

    timer.Start();
    fprintf(stderr, "Pruning %zu samples.\n", sample_names.size());
    auto subtree = MAT::get_tree_copy(T);
    for (auto s: sample_names) {
        if (subtree.get_node(s) == NULL) {
            fprintf(stderr, "ERROR: Sample %s not found in the tree!\n", s.c_str()); 
        }
        else {
            assert (subtree.get_node(s)->is_leaf());
            subtree.remove_node(s, true);
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_subtree (const MAT::Tree& T, std::vector<std::string> sample_names) {

    timer.Start();
    fprintf(stderr, "Extracting subtree of %zu samples.\n", sample_names.size());
    auto subtree = MAT::get_subtree(T, sample_names);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_prune (const MAT::Tree& T, std::vector<std::string> sample_names) {
    
    timer.Start();
    fprintf(stderr, "Large sample input; building subtree by pruning all but %zu samples.\n", sample_names.size());
    //its going to be much quicker to check membership in a set then a vector O(logN)
    //which is important when we're looking at large sample inputs
    //so convert the sample names vector to an unordered set
    std::unordered_set<std::string> setnames(sample_names.begin(),sample_names.end());

    auto subtree = MAT::get_tree_copy(T);
    auto dfs = T.depth_first_expansion();
    for (auto s: dfs) {
        //only call the remover on leaf nodes (can't be deleting the root...)
        if (s->is_leaf()) {
            //if the node is NOT in the set, remove it
            if (setnames.find(s->identifier) == setnames.end()) {
                subtree.remove_node(s->identifier, true);
            }
        }
    }    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    return subtree;
}
