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
                //BUG NOTE: I'm setting the move to false here to patch over this problem
                //if its set to true, then sometimes it will fail to save properly- overwriting root?
                //leave this as not using the move for now, but needs to be revisited
                subtree.remove_node(s->identifier, false);
            }
        }
    }    
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    return subtree;
}

void resolve_polytomy(MAT::Tree* T, std::vector<MAT::Node*> polytomy_nodes) {
    assert (polytomy_nodes.size() > 2);
    std::vector<MAT::Node*> new_parents;
    auto og_parent = polytomy_nodes[0]->parent;
    new_parents.push_back(og_parent);
    //generate a number of new internal nodes to hang these samples on
    //the number is the total number of polytomy nodes - 2 (remember, minimum 3 nodes to resolve)
    for (size_t i=0; i<polytomy_nodes.size()-2; i++) {
        //create a unique identifier indicating that this is a polytomy resolving node
        const std::string nid = og_parent->identifier + "_polytomy_" + std::to_string(i);
        //generate the actual node; its parent is going to be the last entry of new_parents
        //for the first one, that's the original parent, for each after, its the last generated
        auto nnode = T->create_node(nid, new_parents.back(), 0);
        //I don't think these are initialized with mutations, but its called in the condenser, so...
        nnode->clear_mutations();
        new_parents.push_back(nnode);
    }
    //now, assign each of the polytomy samples to one of these new parents
    //the first polytomy node can be ignored and stay a child to its current parent
    //the rest are going to be moved to the parent which matches their index
    //except for the last one, which will go to the last parent in line (index - 1),
    //as the last parent can support two children.
    for (size_t i=1; i<polytomy_nodes.size()-1; i++) {
        //the node mover routine resets the branch length and we don't want to do that
        //save it, move, then reassign the branch length
        float length_to_keep = polytomy_nodes[i]->branch_length;
        T->move_node(polytomy_nodes[i]->identifier, new_parents[i]->identifier, false);
        polytomy_nodes[i]->branch_length = length_to_keep;
    }
    float length_to_keep = polytomy_nodes.back()->branch_length;
    T->move_node(polytomy_nodes.back()->identifier, new_parents.back()->identifier, false);
    polytomy_nodes.back()->branch_length = length_to_keep;
}

MAT::Tree resolve_all_polytomies(MAT::Tree T) {
    //go through the tree, identify every uncondensed polytomy
    //then resolve the polytomy positionally and save the results
    //this will conflict with condensing the tree but that's fine.
    for (auto l: T.get_leaves()) {
        //for every leaf, check the leaf's parent's children
        //parallel to the implementation for node condensing in the class definition
        //adding lots of safety checks to minimize strange behavior, though leaves we care about should always pass these
        //same checks as the node condenser routine.
        if (l == NULL || l->mutations.size() > 0) {
            continue;
        }
        std::vector<MAT::Node*> polytomy_nodes;
        for (auto l2: l->parent->children) {
            if (l2->is_leaf() && (T.get_node(l2->identifier) != NULL) && (l2->mutations.size() == 0)) {
                polytomy_nodes.push_back(l2);
            }
        }
        //greater than 2 because a polytomy of two nodes is already resolved (unlike condenser, which can condense a polytomy of 2)
        if (polytomy_nodes.size() > 2) {
            //pass by reference, modify in place
            resolve_polytomy(&T, polytomy_nodes);
        }
    }
    return T;
}
