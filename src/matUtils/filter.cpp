#include "filter.hpp"

/*
Functions in this module take a MAT and a set of samples and return the a smaller MAT
representing those samples in an optimal manner.
*/
MAT::Tree filter_master(const MAT::Tree& T, std::vector<std::string> sample_names, bool prune, bool keep_clade_annotations) {
    MAT::Tree subtree;
    assert (sample_names.size() > 0);
    if (prune) {
        subtree = prune_leaves(T, sample_names);
    } else if (sample_names.size() < 10000) {
        //for retaining only a subtree, get_subtree is the most effective method
        //for subtress up to about 10 thousand samples in size; after that, pruning
        //all other nodes becomes more efficient, because get_subtree scales with
        //the size of the input sample set, while prune takes a similar time for
        //any sample set size, while scaling on total tree size
        subtree = get_sample_subtree(T, sample_names, keep_clade_annotations);
    } else {
        subtree = get_sample_prune(T, sample_names, keep_clade_annotations);
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
        } else {
            assert (subtree.get_node(s)->is_leaf());
            subtree.remove_node(s, true);
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_subtree (const MAT::Tree& T, std::vector<std::string> sample_names, bool keep_clade_annotations) {

    timer.Start();
    fprintf(stderr, "Extracting subtree of %zu samples.\n", sample_names.size());
    auto subtree = MAT::get_subtree(T, sample_names, keep_clade_annotations);
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());

    return subtree;
}

MAT::Tree get_sample_prune (const MAT::Tree& T, std::vector<std::string> sample_names, bool keep_clade_annotations) {

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
    if (!keep_clade_annotations) {
        dfs = subtree.depth_first_expansion();
        for (auto s: dfs) {
            s->clear_annotations();
        }
    }
    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
    return subtree;
}

void resolve_polytomy(MAT::Tree* T, std::vector<MAT::Node*> polytomy_nodes) {
    assert (polytomy_nodes.size() > 2);
    // fprintf(stderr, "DEBUG: Resolving polytomy with %ld children\n", polytomy_nodes.size());
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
        T->move_node(polytomy_nodes[i]->identifier, new_parents[i]->identifier, false);
    }
    T->move_node(polytomy_nodes.back()->identifier, new_parents.back()->identifier, false);
}

MAT::Tree resolve_all_polytomies(MAT::Tree T) {
    //go through the tree, identify every uncondensed polytomy
    //then resolve the polytomy positionally and save the results
    //this will conflict with condensing the tree but that's fine.
    //step one is go through the tree and delete all single children nodes.
    T.collapse_tree();
    for (auto n: T.depth_first_expansion()) {
        //greater than 2 because a polytomy of two nodes is already resolved (unlike condenser, which can condense a polytomy of 2)
        if (n->children.size() > 2) {
            //pass by reference, modify in place
            resolve_polytomy(&T, n->children);
        }
    }
    return T;
}

void reroot_tree(MAT::Tree* T, std::string rnid) {
    //to reroot, we traverse from the root down to the new root via an rsearch
    //moving the parent to be the child of the next node in each up to the final node.
    auto norder = T->rsearch(rnid,true);
    if (norder.size() == 0) {
        fprintf(stderr, "ERROR: New root selection not found in tree. Exiting\n");
        exit(1);
    }
    fprintf(stderr, "Moving root to node %s over a distance of %ld connections.\n", rnid.c_str(), norder.size());
    std::reverse(norder.begin(), norder.end()); //reverse so the root is first.
    for (size_t i = 0; i < norder.size() - 1; i++) {

        MAT::Node* source = norder[i];
        MAT::Node* destination = norder[i+1];
        source->parent = destination;

        destination->parent = NULL;
        T->root = destination;

        auto iter = std::find(source->children.begin(), source->children.end(), destination);
        assert (iter != source->children.end());
        source->children.erase(iter);
        destination->children.push_back(source);
        assert (source->parent->identifier == destination->identifier);
        assert (destination->parent == NULL);
        //update the level values for all descendents
        destination->level = 0;
        std::queue<MAT::Node*> remaining_nodes;
        for (auto c: destination->children) {
            remaining_nodes.push(c);
        }
        while (remaining_nodes.size() > 0) {
            MAT::Node* curr_node = remaining_nodes.front();
            remaining_nodes.pop();
            curr_node->level = curr_node->parent->level + 1;
            for (auto c: curr_node->children) {
                remaining_nodes.push(c);
            }
        }
    }
    assert (T->get_node(rnid)->is_root());
}