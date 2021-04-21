#include "process_each_node.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>
    int get_parsimmony_score_only(MAT::Node* src, MAT::Node* dst);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
static void check_LCA(MAT::Node* src, MAT::Node* dst,MAT::Node* LCA_to_check){
    std::unordered_set<int> path;
    while (src!=LCA_to_check) {
        path.insert(src->bfs_index);
        src=src->parent;
    }
    while (dst!=LCA_to_check) {
        assert(!path.count(dst->bfs_index));
        dst=dst->parent;
    }
    assert(dst==LCA_to_check);
}
#endif
static void update_debug(MAT::Node* cur, MAT::Node* LCA,std::vector<MAT::Node*>& dst_stack
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
 ,std::vector<state_change_hist_dbg>& debug
#endif
){
    while (cur!=LCA) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        auto cur_mut_iter=cur->mutations.begin();
        auto cur_mut_end=cur->mutations.end();
        for(state_change_hist_dbg& ele:debug){
            while (cur_mut_iter!=cur_mut_end&&cur_mut_iter->get_position()<ele.position) {
                cur_mut_iter++;
            }
            if (cur_mut_iter!=cur_mut_end&&cur_mut_iter->get_position()==ele.position) {
                ele.count_change.push_back(Mutation_Count_Change());
                ele.mutation_score_change.push_back(0);
                ele.par_nuc.push_back(cur_mut_iter->get_par_one_hot());
                ele.major_allele.push_back(cur_mut_iter->get_mut_one_hot()|cur_mut_iter->get_tie_one_hot());
            }else {
                ele.count_change.push_back(Mutation_Count_Change());
                ele.mutation_score_change.push_back(0);
                ele.par_nuc.push_back(ele.par_nuc.back());
                ele.major_allele.push_back(ele.par_nuc.back());
            }
        }
#endif
        dst_stack.push_back(cur);
        cur=cur->parent;
    }
}

static int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    // const New_Tie_Collection_t &unresolved_ties_from_LCA,
    int parsimony_score_change,output_t& output,
    const std::vector<MAT::Node *> node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> debug_from_src
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    check_LCA(src, dst, LCA);
    assert(src->parent!=dst);
    assert(src!=dst);
    std::vector<state_change_hist_dbg> debug_from_dst;
#endif
    std::vector<MAT::Node *> node_stack_from_dst({});

    // No removing mutation even if adding to leaf, as all mutations in the
    // original leaf will be present, and it will at least be tie, so contribute
    // to count by 1. However, removing mutations shared with original leaf is
    // still needed.
    // New_Tie_Collection_t ties;
    MAT::Node *this_node = dst;
    assert(dst);
    Mutation_Count_Change_Collection parent_added;
    Mutation_Count_Change_Collection parent_of_parent_added;
    // New_Tie_Collection_t new_ties;
    if (LCA != dst) {

        node_stack_from_dst.push_back(this_node);
        if (this_node != LCA) {
            get_parsimony_score_change_from_add(
                this_node, mutations, parent_added, true, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,
                debug_from_dst, node_stack_from_dst
#endif
            );
            this_node = this_node->parent;
            while (this_node != LCA) {

                get_intermediate_nodes_mutations(
                    this_node, parent_added, parent_of_parent_added,
                    parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                    ,
                    debug_from_dst, node_stack_from_dst
#endif
                );
                node_stack_from_dst.push_back(this_node);
                parent_added = std::move(parent_of_parent_added);
                if (parent_added.empty()) {
                    update_debug(this_node->parent, LCA, node_stack_from_dst
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
,debug_from_dst
#endif
                    );
                    break;
                }

                this_node = this_node->parent;
                // ties = std::move(new_ties);
            }
        }
    }

#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(debug_from_dst.empty() ||
           debug_from_dst[0].mutation_score_change.size() ==
               node_stack_from_dst.size());
    std::vector<LCA_merged_states> merged_states;
    std::vector<state_change_hist_dbg> debug_above_LCA;
    prep_LCA_checker(debug_from_src, debug_from_dst, merged_states);
#endif
    std::vector<MAT::Node *> node_stack_above_LCA;

    MAT::Node *ancestor = LCA->parent;
    bool is_src_terminal = src->parent == LCA;
    bool is_dst_terminal = dst == LCA;
    if ((!(root_mutations_altered.empty() && parent_added.empty())) ||
        is_src_terminal || is_dst_terminal) {
        get_LCA_mutation(LCA, is_src_terminal?src:0, is_dst_terminal?dst:0,
                         root_mutations_altered, is_dst_terminal?mutations:parent_added,
                         parent_of_parent_added, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                         ,debug_above_LCA, merged_states, node_stack_above_LCA,
                         is_src_terminal ? src : node_stack_from_src.back(),
                         is_dst_terminal?dst:node_stack_from_dst.back()
#endif
        );
        node_stack_above_LCA.push_back(LCA);
        ancestor = LCA->parent;
        parent_added = std::move(parent_of_parent_added);
        while (ancestor && (!parent_added.empty())) {
            get_intermediate_nodes_mutations(
                ancestor, parent_added, parent_of_parent_added,
                parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                ,
                debug_above_LCA, node_stack_above_LCA
#endif
            );
            node_stack_above_LCA.push_back(ancestor);
            parent_added = std::move(parent_of_parent_added);
            ancestor=ancestor->parent;
        }
        if (!ancestor) {
            for(auto & a:parent_added){
                parsimony_score_change+=a.get_default_change_internal();
            }
        }
    }
    else {
        node_stack_above_LCA.push_back(LCA);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        auto iter=LCA->mutations.begin();
        auto end=LCA->mutations.end();
        for(auto& in:merged_states){
            while (iter!=end&&iter->get_position()<in.position) {
                iter++;
            }
            if(iter!=end&&iter->get_position()==in.position){
                debug_above_LCA.emplace_back(iter->get_position(),iter->get_par_one_hot(), iter->get_mut_one_hot()|iter->get_tie_one_hot(), 0,Mutation_Count_Change());
            }else {
                debug_above_LCA.emplace_back(in.position, in.par_allele,in.par_allele, 0,Mutation_Count_Change());
            }
        }
#endif
    }
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    //fprintf(stderr, "LCA idx: %zu",LCA_a->bfs_index);
    int ref_parsimony_score = get_parsimmony_score_dumb(
        ancestor ? ancestor : node_stack_above_LCA.back(), src, dst,
        debug_from_src, node_stack_from_src, debug_from_dst,
        node_stack_from_dst, debug_above_LCA, node_stack_above_LCA);
    assert(ref_parsimony_score == parsimony_score_change);
#endif
    assert((dst==LCA&&node_stack_above_LCA[0]==LCA)||node_stack_from_dst[0]==dst);
    if (parsimony_score_change<=output.score_change) {
        if (parsimony_score_change<output.score_change) {
            output.score_change=parsimony_score_change;
            for(auto t:output.moves){
                delete t;
            }
            output.moves.clear();
        }
        Profitable_Moves_ptr_t new_move=new Profitable_Moves;
        new_move->score_change=parsimony_score_change;
        new_move->src_to_LCA=std::move(node_stack_from_src);
        new_move->dst_to_LCA=std::move(node_stack_from_dst);
        new_move->src=src;
        new_move->LCA=LCA;
        assert(dst==LCA||new_move->dst_to_LCA.back()->parent==LCA);
        assert(src->parent==LCA||new_move->src_to_LCA.back()->parent==LCA);
        for(const auto node:node_stack_above_LCA){
            new_move->dst_to_LCA.push_back(node);
        }
        output.moves.push_back(new_move);
    }
    return parsimony_score_change;
}

/*
static int check_move_profitable(
    MAT::Node *src, MAT::Node *dst, MAT::Node *LCA,
    const Mutation_Count_Change_Collection &mutations,
    const Mutation_Count_Change_Collection &root_mutations_altered,
    // const New_Tie_Collection_t &unresolved_ties_from_LCA,
    int parsimony_score_change,output_t& output,
    const std::vector<MAT::Node *> node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,
    const std::vector<state_change_hist_dbg> debug_from_src
#endif
) {
    int parsimony_score_change__=check_move_profitable_internal(src, dst, LCA,mutations,
    root_mutations_altered,parsimony_score_change,output,node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,debug_from_src
#endif
    );
    int ref=individual_move(src,dst,LCA);
    assert(parsimony_score_change__==ref);
    return parsimony_score_change__;
}*/
static void
merge_mutation_LCA_to_dst(MAT::Node *child,
                          const Mutation_Count_Change_Collection &mutations,
                          Mutation_Count_Change_Collection &child_mutations) {
    auto child_mutation_iter = child->mutations.begin();
    auto child_mutation_end = child->mutations.end();
    for (const Mutation_Count_Change &m : mutations) {
        while (child_mutation_iter != child_mutation_end &&
               !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }
        while (child_mutation_iter != child_mutation_end &&
               *child_mutation_iter < m) {
            child_mutations.push_back(*child_mutation_iter);
            child_mutations.back().set_change(
                0, child_mutation_iter->get_par_one_hot());
            child_mutations.back().set_par_nuc(
                child_mutation_iter->get_mut_one_hot());
            child_mutation_iter++;
        }
        if (child_mutation_iter != child_mutation_end &&
            child_mutation_iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_allele = child_mutation_iter->get_mut_one_hot();
                child_mutations.emplace_back(m);
                child_mutations.back().set_par_nuc(new_par_allele);
            child_mutation_iter++;
        } else {
            child_mutations.emplace_back(m);
        }
    }
    while (child_mutation_iter != child_mutation_end) {
        child_mutations.push_back(*child_mutation_iter);
        child_mutations.back().set_change(
            0, child_mutation_iter->get_par_one_hot());
        child_mutations.back().set_par_nuc(
            child_mutation_iter->get_mut_one_hot());
        child_mutation_iter++;
    }
}
void search_subtree(
    MAT::Node *root, MAT::Node *child_exclude, MAT::Node *LCA, MAT::Node *src,
    Mutation_Count_Change_Collection mutations,
    output_t &profitable_moves,
    const Mutation_Count_Change_Collection &root_mutations_altered,int radius,
    int parsimony_score_change_from_removal,const std::vector<MAT::Node *> &node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    ,const std::vector<state_change_hist_dbg> &debug_from_src
#endif
) {
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(debug_from_src.empty() ||
           debug_from_src[0].mutation_score_change.size() ==
               node_stack_from_src.size());
               #endif
    if(!(root->children.size()==1&&root->children[0]->is_leaf())&&radius>0){
    for (MAT::Node *child : root->children) {
        if (child != child_exclude) {
            Mutation_Count_Change_Collection child_mutations;
            merge_mutation_LCA_to_dst(child, mutations, child_mutations);
            search_subtree(child, nullptr, LCA, src, child_mutations,
                           profitable_moves, root_mutations_altered,radius-1,
                           parsimony_score_change_from_removal, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                           ,debug_from_src
#endif
            );
        }
    }}
    if (src->parent != root) {
        auto score_change = check_move_profitable(
            src, root, LCA, mutations, root_mutations_altered,
            parsimony_score_change_from_removal,profitable_moves, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,debug_from_src
#endif
        );
    }
}
static void
merge_mutation_src_to_LCA(const MAT::Node *root,
                          Mutation_Count_Change_Collection &mutations) {
    Mutation_Count_Change_Collection merged_mutations;
    merged_mutations.reserve(mutations.size() + root->mutations.size());
    auto iter = mutations.begin();
    auto end = mutations.end();
    for (const auto &m : root->mutations) {
        if (!m.is_valid()) {
            continue;
        }
        while (iter != end && iter->get_position() < m.get_position()) {
            merged_mutations.push_back(*iter);
            iter++;
        }
        if (iter != end && iter->get_position() == m.get_position()) {
            nuc_one_hot new_par_nuc = m.get_par_one_hot();
            if (new_par_nuc != iter->get_incremented()) {
                merged_mutations.push_back(*iter);
                merged_mutations.back().set_par_nuc(new_par_nuc);
            }
            iter++;
        } else {
            merged_mutations.emplace_back(m);
            merged_mutations.back().set_change(0, m.get_mut_one_hot());
        }
    }
    while (iter!=end) {
        merged_mutations.push_back(*iter);
        iter++;
    }
    mutations = std::move(merged_mutations);
}
static void init_mutation_change(MAT::Node* src, Mutation_Count_Change_Collection& mutations){
    mutations.reserve(src->mutations.size());
    for (const auto &m : src->mutations) {
        if (m.is_valid()||m.get_tie_one_hot()) {
            mutations.emplace_back(m);
            mutations.back().set_change(0, m.get_mut_one_hot() |
                                               m.get_tie_one_hot());
        }
    }
}
void find_profitable_moves(MAT::Node *src, output_t &out,int radius) {
    if (src->is_root()) {
        return;
    }
    if(src->parent->children.size()<=1){
        return;
    }
    MAT::Node *root = src->parent;
    MAT::Node *exclude = src;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection alter_mutations;
    Mutation_Count_Change_Collection new_alter_mutations;
    // New_Tie_Collection_t unresolved_ties_from_LCA;
    // New_Tie_Collection_t new_tie;
    std::vector<MAT::Node *> node_stack;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<state_change_hist_dbg> debug;
#endif

    // Do searching on the current parent of src directly
    {
        search_subtree(root, exclude, root, src, mutations, out,
                       Mutation_Count_Change_Collection(),radius, 0,
                       std::vector<MAT::Node *>({})
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,std::vector<state_change_hist_dbg>()
                       #endif
    );
    }

    get_parent_altered_remove(alter_mutations, src, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                              ,
                              debug, node_stack
#endif
    );

    node_stack.push_back(src->parent);
    merge_mutation_src_to_LCA(root, mutations);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(std::max(src->mutations.size(), src->parent->mutations.size()) <=
           debug.size());
    size_t old_debug_size = debug.size();
#endif
    exclude=root;
    root=root->parent;
    while (root&&radius) {
        search_subtree(root, exclude, root, src, mutations, out,
                       alter_mutations,radius, parsimony_score_change,
                       std::vector<MAT::Node *>(node_stack)
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                       ,
                       std::vector<state_change_hist_dbg>(debug)
#endif
        );

        get_intermediate_nodes_mutations(
            root, alter_mutations, new_alter_mutations, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,
            debug, node_stack
#endif
        );
        node_stack.push_back(root);
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        assert(debug.size() == old_debug_size);
#endif
        // mutations.merge(root->mutations,
        // MAT::Mutations_Collection::KEEP_SELF);
        merge_mutation_src_to_LCA(root, mutations);
        alter_mutations = std::move(new_alter_mutations);
        // unresolved_ties_from_LCA = std::move(new_tie);
        exclude = root;
        root = root->parent;
        radius--;
    }
}
int individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA){
MAT::Node *root = src->parent;
    Mutation_Count_Change_Collection mutations;
    init_mutation_change(src, mutations);
    Mutation_Count_Change_Collection root_mutations_altered;
    Mutation_Count_Change_Collection new_alter_mutations;
    // New_Tie_Collection_t unresolved_ties_from_LCA;
    // New_Tie_Collection_t new_tie;
    int parsimony_score_change = 0;
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    std::vector<state_change_hist_dbg> debug_from_src;
#endif
    std::vector<MAT::Node *> node_stack_from_src;
    size_t old_debug_size;
    if(src->parent->children.size()<=1){
        return 0;
    }
    if(root!=LCA){
    get_parent_altered_remove(root_mutations_altered, src, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
                              ,debug_from_src, node_stack_from_src
#endif
    );
    node_stack_from_src.push_back(src->parent);
    #ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    assert(std::max(src->mutations.size(), src->parent->mutations.size()) <=
           debug_from_src.size());
    old_debug_size = debug_from_src.size();
    merge_mutation_src_to_LCA(root, mutations);
    root=root->parent;
#endif
    }


    while (root&&root!=LCA) {
        get_intermediate_nodes_mutations(
            root, root_mutations_altered, new_alter_mutations, parsimony_score_change
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
            ,
            debug_from_src, node_stack_from_src
#endif
        );
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
        assert(debug_from_src.size() == old_debug_size);
#endif
        node_stack_from_src.push_back(root);
        // mutations.merge(root->mutations,
        // MAT::Mutations_Collection::KEEP_SELF);
        merge_mutation_src_to_LCA(root, mutations);
        root_mutations_altered = std::move(new_alter_mutations);
        // unresolved_ties_from_LCA = std::move(new_tie);
        root = root->parent;
    }
    std::vector<MAT::Node*> dst_path;
    root=dst;
    while (root!=LCA) {
        dst_path.push_back(root);
        root=root->parent;
    }
    for (int i=dst_path.size()-1; i>=0; i--) {
        Mutation_Count_Change_Collection child_mutations;
        merge_mutation_LCA_to_dst(dst_path[i], mutations, child_mutations);
        mutations=std::move(child_mutations);
    }
    output_t out;
    return check_move_profitable(src, dst, LCA, mutations, root_mutations_altered, parsimony_score_change,out, node_stack_from_src
#ifdef DEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
    , debug_from_src
#endif
    );
}