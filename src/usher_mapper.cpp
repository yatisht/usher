#include "usher_graph.hpp"

tbb::mutex data_lock;
tbb::reader_writer_lock rd_wr_lock;

int mapper_body::operator()(mapper_input input) {
    TIMEIT();

    // Valid variants should have positions >= 0
    if (input.variant_pos >= 0) {
        data_lock.lock();
        fprintf(stderr, "At variant site %i\n", input.variant_pos);
        data_lock.unlock();

        size_t num_nodes = input.bfs->size();

        // State and score vector for the Fitch-Sankoff algorithm. State uses
        // last 4 bits of int8 to mark which of the four bases: A,C,G,T is
        // parsimony-assigned. scores vector maintains a vector of 4 values, one
        // for each possible base (A,C,G,T)
        std::vector<int8_t> states(num_nodes);
        std::vector<std::vector<int>> scores(num_nodes);

        // Initialize scores and states to 0 values
        for (size_t i=0; i<num_nodes; i++) {
            scores[i].resize(4);
            for (int j=0; j<4; j++) {
                scores[i][j] = 0;
            }
            states[i] = 0;
        }

        // For leaf nodes in the tree, start by intializing scores for
        // non-reference bases to a large value (num_nodes) here
        for (auto l: input.T->get_leaves()) {
            size_t node_idx = (*input.bfs_idx)[l->identifier];
            int8_t ref_nuc_id = MAT::get_nt(input.ref_nuc);
            assert (ref_nuc_id >= 0);
            for (int j=0; j<4; j++) {
                if (j != ref_nuc_id) {
                    scores[node_idx][j] = (int) num_nodes;
                }
            }
        }

        // Iterate over the variants
        for (auto v: input.variants) {
            size_t pos = std::get<0> (v);
            int8_t nuc = std::get<1> (v);
            std::string nid = (*input.variant_ids)[pos];
            // If variant is for the tree sample
            if (input.bfs_idx->find(nid) != input.bfs_idx->end()) {
                size_t idx= (*input.bfs_idx)[nid];
                // Initialize scores for bases corresponding to alternate
                // alleles to 0 and remaining to a high value (num_nodes)
                for (int j=0; j<4; j++) {
                    scores[idx][j] = (int) num_nodes;
                    if (((1 << j) & nuc) != 0) {
                        scores[idx][j] = 0;
                    }
                }
            }
            // If variant is for the missing sample to be placed, simply add the
            // variant to the sample mutation list
            else {
                //auto mutations_iter = input.missing_sample_mutations->begin() + (iter - input.missing_samples->begin());
                MAT::Mutation m;
                m.chrom = input.chrom;
                m.position = input.variant_pos;
                m.ref_nuc = input.ref_nuc;
                if (nuc == 15) {
                    m.is_missing = true;
                } else {
                    m.is_missing = false;
                    assert ((nuc > 0) && (nuc < 15));
                    m.mut_nuc = nuc;
                }
                auto iter = std::find(input.missing_samples->begin(), input.missing_samples->end(), nid);
                data_lock.lock();
                (*iter).mutations.emplace_back(m);
                data_lock.unlock();
            }
        }

        // Sankoff: forward pass
        for (auto it=(*input.bfs).rbegin(); it!=(*input.bfs).rend(); it++) {
            auto node = (*it);
            auto node_idx = (*input.bfs_idx)[node->identifier];

            if (!node->is_leaf()) {
                for (auto child: (*node).children) {
                    auto c_id = child->identifier;
                    auto c_idx = (*input.bfs_idx)[c_id];
                    for (int j=0; j<4; j++) {
                        int min_s = (int) num_nodes+1;
                        for (int k=0; k<4; k++) {
                            int c_s;
                            if (k==j) {
                                c_s = scores[c_idx][k];
                            } else {
                                c_s = scores[c_idx][k]+1;
                            }
                            if (c_s < min_s) {
                                min_s = c_s;
                            }
                        }
                        scores[node_idx][j] += min_s;
                    }
                }
            }
        }

        // Sankoff: backward pass
        for (auto it=(*input.bfs).begin(); it!=(*input.bfs).end(); it++) {
            auto node = (*it);
            auto node_idx = (*input.bfs_idx)[node->identifier];

            int8_t par_state = 0;
            if (node->parent != NULL) {
                auto par = node->parent;
                auto par_id = par->identifier;
                auto par_idx = (*input.bfs_idx)[par_id];
                par_state = states[par_idx];
            } else {
                par_state = MAT::get_nt(input.ref_nuc);
                assert (par_state >= 0);
            }

            int8_t state = par_state;
            int min_s = scores[node_idx][par_state];
            for (int j=0; j<4; j++) {
                if (scores[node_idx][j] < min_s) {
                    min_s = scores[node_idx][j];
                    state = j;
                }
            }
            if (state != par_state) {
                if (scores[node_idx][par_state] == min_s) {
                    state = par_state;
                }
            }
            states[node_idx] = state;

            if (state != par_state) {
                MAT::Mutation m;
                m.chrom = input.chrom;
                m.position = input.variant_pos;
                m.ref_nuc = input.ref_nuc;
                m.par_nuc = (1 << par_state);
                m.mut_nuc  = (1 << state);

                data_lock.lock();
                node->add_mutation(m);
                data_lock.unlock();
            }
        }

    }

    return 1;
}

// Used to do a parallel search for the parsimony-optimal placement node. If
// compute_parsimony_scores is not set, the function can return early if the
// parsimony score at the current input node exceeds the smallest parsimony
// score encountered during the parallel search
void mapper2_body(mapper2_input& input, bool compute_parsimony_scores, bool compute_vecs) {
    //    TIMEIT();

    // Variable to store the number of parsimony-increasing mutations to
    // place sample at the current node
    int set_difference = 0;

    // Current smallest value of the number of parsimony-increasing mutations
    // during the parallel search to place the same at some node in the tree
    int best_set_difference = *input.best_set_difference;

    std::vector<int> anc_positions;
    std::vector<MAT::Mutation> ancestral_mutations;

    // if node has some unique mutations not in new sample, placement should be
    // done as a sibling
    bool has_unique = false;
    int node_num_mut = 0;
    int num_common_mut = 0;

    // For non-root nodes, add mutations common to current node (branch) to
    // excess mutations. Set has_unique to true if a mutation unique to current
    // node not in new sample is found.
    if (!input.node->is_root()) {
        size_t start_index = 0;
        for (auto m1: input.node->mutations) {
            node_num_mut++;
            auto anc_nuc = m1.mut_nuc;
            // if mutation is masked, treat it as a unique mutation (add as
            // sibling)
            if (m1.is_masked()) {
                has_unique = true;
                break;
            }
            assert (((anc_nuc-1) & anc_nuc) == 0);
            bool found = false;
            bool found_pos = false;
            for (size_t k = start_index; k < input.missing_sample_mutations->size(); k++) {
                auto m2 = (*input.missing_sample_mutations)[k];
                start_index = k;
                if (m1.position == m2.position) {
                    found_pos = true;
                    if (m2.is_missing) {
                        found = true;
                        num_common_mut++;
                    } else {
                        auto nuc = m2.mut_nuc;
                        if ((nuc & anc_nuc) != 0) {
                            MAT::Mutation m;
                            m.chrom = m1.chrom;
                            m.position = m1.position;
                            m.ref_nuc = m1.ref_nuc;
                            m.par_nuc = m1.par_nuc;
                            m.mut_nuc = anc_nuc;

                            ancestral_mutations.emplace_back(m);
                            anc_positions.emplace_back(m.position);
                            assert((m.mut_nuc & (m.mut_nuc-1)) == 0);
                            if (compute_vecs) {
                                (*input.excess_mutations).emplace_back(m);
                            }

                            // Ambiguous base
                            //if ((nuc & (nuc-1)) != 0) {
                            //    (*input.imputed_mutations).emplace_back(m);
                            //}
                            found = true;
                            num_common_mut++;
                            break;
                        }
                    }
                }
                if (m1.position < m2.position) {
                    break;
                }
            }
            if (!found) {
                if (!found_pos && (anc_nuc == m1.ref_nuc)) {
                    MAT::Mutation m;
                    m.position = m1.position;
                    m.chrom = m1.chrom;
                    m.ref_nuc = m1.ref_nuc;
                    m.par_nuc = m1.par_nuc;
                    m.mut_nuc = anc_nuc;

                    ancestral_mutations.emplace_back(m);
                    anc_positions.emplace_back(m.position);
                    assert((m.mut_nuc & (m.mut_nuc-1)) == 0);
                    if (compute_vecs) {
                        (*input.excess_mutations).emplace_back(m);
                    }

                    num_common_mut++;
                } else {
                    has_unique = true;
                }
            }
        }
    } else {
        for (auto m: input.node->mutations) {
            ancestral_mutations.emplace_back(m);
            anc_positions.emplace_back(m.position);
        }
    }

    // Add ancestral mutations to ancestral mutations. When multiple mutations
    // at same position are found in the path leading from the root to the
    // current node, add only the most recent mutation to the vector
    {
        auto n = input.node;
        while (n->parent != NULL) {
            n = n->parent;
            for (auto m: n->mutations) {
                if (!m.is_masked() && (std::find(anc_positions.begin(), anc_positions.end(), m.position) == anc_positions.end())) {
                    ancestral_mutations.emplace_back(m);
                    anc_positions.emplace_back(m.position);
                }
            }
        }
    }

    // sort by position. This helps speed up the search
    std::sort(ancestral_mutations.begin(), ancestral_mutations.end());

    // Iterate over missing sample mutations
    for (auto m1: (*input.missing_sample_mutations)) {
        // Missing bases (Ns) are ignored
        if (m1.is_missing) {
            continue;
        }
        size_t start_index = 0;
        bool found_pos = false;
        bool found = false;
        bool has_ref = false;
        auto anc_nuc = m1.ref_nuc;
        if ((m1.mut_nuc & m1.ref_nuc) != 0) {
            has_ref = true;
        }
        // Check if mutation is found in ancestral_mutations
        for (size_t k = start_index; k < ancestral_mutations.size(); k++) {
            auto m2 = ancestral_mutations[k];
            // Masked mutations don't match anything
            if (m2.is_masked()) {
                continue;
            }
            start_index = k;
            if (m1.position == m2.position) {
                found_pos = true;
                anc_nuc = m2.mut_nuc;
                if ((m1.mut_nuc & anc_nuc) != 0) {
                    found = true;
                }
                break;
            }
        }
        if (found) {
            // If mutation is found in ancestral_mutations
            // and if the missing sample base was ambiguous,
            // add it to imputed_mutations
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0)) {
                MAT::Mutation m;
                m.chrom = m1.chrom;
                m.position = m1.position;
                m.ref_nuc = m1.ref_nuc;
                m.par_nuc = anc_nuc;
                m.mut_nuc = anc_nuc;
                input.imputed_mutations->emplace_back(m);
            }
        }
        // If neither the same mutation nor another mutation at the same
        // position is found in ancestor but if the missing sample can carry
        // the reference allele, add a mutation with reference allele to
        // imputed_mutations for the sample (it's not a parsimony-increasing
        // mutation)
        else if (!found_pos && has_ref) {
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0)) {
                MAT::Mutation m;
                m.chrom = m1.chrom;
                m.position = m1.position;
                m.ref_nuc = m1.ref_nuc;
                m.par_nuc = anc_nuc;
                m.mut_nuc = m1.ref_nuc;
                input.imputed_mutations->emplace_back(m);
            }
        }
        // In all other cases, it is a parsimony-increasing mutation. Return
        // early if number of parsimony-increasing mutations exceeds the current
        // best. Otherwise add the mutation to excess_mutations and to
        // imputed_mutations, if base was originally ambiguous
        else {
            MAT::Mutation m;
            m.chrom = m1.chrom;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            m.par_nuc = anc_nuc;
            if (has_ref) {
                m.mut_nuc = m1.ref_nuc;
            } else {
                for (int j=0; j < 4; j++) {
                    if (((1 << j) & m1.mut_nuc) != 0) {
                        m.mut_nuc = (1 << j);;
                        break;
                    }
                }
            }
            assert((m.mut_nuc & (m.mut_nuc-1)) == 0);
            // If the missing sample base is ambiguous, add it to
            // imputed_mutations
            if (compute_vecs && ((m1.mut_nuc & (m1.mut_nuc - 1)) != 0)) {
                input.imputed_mutations->emplace_back(m);
            }
            if (m.mut_nuc != m.par_nuc) {
                if (compute_vecs) {
                    input.excess_mutations->emplace_back(m);
                }
                set_difference += 1;
                if (!compute_parsimony_scores && (set_difference > best_set_difference)) {
                    return;
                }
            }
        }
    }

    // For loop to add back-mutations for cases in which a mutation from the
    // root to the current node consists of a non-reference allele but no such
    // variant is found in the missing sample
    for (auto m1: ancestral_mutations) {
        size_t start_index = 0;
        bool found = false;
        bool found_pos = false;
        auto anc_nuc = m1.mut_nuc;
        for (size_t k = start_index; k < input.missing_sample_mutations->size(); k++) {
            // If ancestral mutation is masked, terminate the search for
            // identical mutation
            if (m1.is_masked()) {
                break;
            }
            auto m2 = (*input.missing_sample_mutations)[k];
            start_index = k;
            if (m1.position == m2.position) {
                found_pos = true;
                // Missing bases (Ns) are ignored
                if (m2.is_missing) {
                    found = true;
                    break;
                }
                if ((m2.mut_nuc & anc_nuc) != 0) {
                    found = true;
                }
            }
        }
        // If ancestral mutation is found, do nothing. Else, if last ancestor
        // with mutation in the same position is found having reference allele, do
        // nothing. If same position is found but not with the same allele, do
        // nothing since the mutation must be added to excess mutations already
        // in the previous loop. In all remaining cases, add the mutation to
        // excess_mutations.
        if (found) {
        } else if (!found_pos && !m1.is_masked() && (anc_nuc == m1.ref_nuc)) {
        } else if (found_pos && !found) {
        } else {
            MAT::Mutation m;
            m.chrom = m1.chrom;
            m.position = m1.position;
            m.ref_nuc = m1.ref_nuc;
            m.par_nuc = anc_nuc;
            m.mut_nuc = m1.ref_nuc;
            assert(m.is_masked() || ((m.mut_nuc & (m.mut_nuc-1)) == 0));
            if (m.mut_nuc != m.par_nuc) {
                set_difference += 1;
                if (!compute_parsimony_scores && (set_difference > best_set_difference)) {
                    return;
                }
                if (compute_vecs) {
                    (*input.excess_mutations).emplace_back(m);
                }
            }
        }
    }

    // Set the number of parsimony-increasing mutations
    if (compute_parsimony_scores) {
        *input.set_difference = set_difference;
    }

    // if sibling of internal node or leaf, ensure it is not equivalent to placing under parent
    // if child of internal node, ensure all internal node mutations are present in the sample
    if (input.node->is_root() || ((has_unique && !input.node->is_leaf() && (num_common_mut > 0) && (node_num_mut != num_common_mut)) || \
                                  (input.node->is_leaf() && (num_common_mut > 0)) || (!has_unique && !input.node->is_leaf() && (node_num_mut == num_common_mut)))) {
        rd_wr_lock.lock_read();
        if (set_difference > *input.best_set_difference) {
            rd_wr_lock.unlock();
            return;
        }
        rd_wr_lock.unlock();

        rd_wr_lock.lock();
        size_t num_leaves = input.T->get_num_leaves(input.node);
        if (set_difference < *input.best_set_difference) {
            *input.best_set_difference = set_difference;
            *input.best_node = input.node;
            *input.best_node_num_leaves = num_leaves;
            *input.best_j = input.j;
            *input.num_best = 1;
            *input.has_unique = has_unique;
            *input.best_distance = input.distance;
            (*input.node_has_unique)[input.j] = has_unique;
            input.best_j_vec->clear();
            input.best_j_vec->emplace_back(input.j);
        } else if (set_difference == *input.best_set_difference) {
            // Tie breaking strategy when multiple parsimony-optimal placements
            // are found. it picks the node with a greater number of descendant
            // leaves for placement. However, if the choice is between a parent
            // and its child node, it picks the parent node if the number of
            // descendant leaves of the parent that are not shared with the child
            // node exceed the number of descendant leaves of the child.
            if (((input.distance == *input.best_distance) &&
                    ((num_leaves > *input.best_node_num_leaves) ||
                     ((num_leaves == *input.best_node_num_leaves) && (*input.best_j < input.j))))
                    || (input.distance < *input.best_distance)) {
                *input.best_set_difference = set_difference;
                *input.best_node = input.node;
                *input.best_node_num_leaves = num_leaves;
                *input.best_j = input.j;
                *input.has_unique = has_unique;
                *input.best_distance = input.distance;
            }
            *input.num_best += 1;
            (*input.node_has_unique)[input.j] = has_unique;
            input.best_j_vec->emplace_back(input.j);
        }
        rd_wr_lock.unlock();
    } else if (compute_parsimony_scores) {
        // Add 1 to the parsimony score for this node
        // as its current best placement is equivalent
        // to placing at parent/child
        *input.set_difference = set_difference + 1;
    }
}

