#include "src/mutation_annotated_tree.hpp"
#include "src/tree_rearrangements/check_samples.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

struct Merge_Discriptor{
    MAT::Node* to_merge_with;
    MAT::Node* this_child;
    MAT::Mutations_Collection shared_mutations;
    MAT::Mutations_Collection src_unique_mutations;
    MAT::Mutations_Collection to_merge_width_unique_mutations;
};

struct Merge_Comparator{
    bool operator()(const Merge_Discriptor& first, const Merge_Discriptor& second)const {
        return first.shared_mutations.size()>second.shared_mutations.size();
    }
};

static bool check_mergable(MAT::Node* this_node, MAT::Node* merge_with,std::vector<Merge_Discriptor>& possible_merges){
    Merge_Discriptor& merge=possible_merges.back();
    this_node->mutations.set_difference(merge_with->mutations, merge.src_unique_mutations, merge.to_merge_width_unique_mutations, merge.shared_mutations);
    if(merge.shared_mutations.size()==0){
        merge.src_unique_mutations.clear();
        merge.to_merge_width_unique_mutations.clear();
        return false;
    }
    possible_merges.emplace_back();
    return true;
}

static void remove_child(MAT::Node* child_to_remove,MAT::Node* parent){
        auto iter=std::find(parent->children.begin(),parent->children.end(),child_to_remove);
        if(iter!=parent->children.end()) parent->children.erase(iter);
}
static MAT::Node* execute_merge(Merge_Discriptor& merge, MAT::Tree* tree,MAT::Node* parent){
    MAT::Node* merged_node=tree->create_node(std::to_string(++tree->curr_internal_node),parent);
    merged_node->mutations.swap(merge.shared_mutations);
    remove_child(merge.to_merge_with,parent);

    merged_node->children.push_back(merge.this_child);
    merge.this_child->parent=merged_node;
    merge.this_child->mutations.swap(merge.src_unique_mutations);

    merged_node->children.push_back(merge.to_merge_with);
    merge.to_merge_with->parent=merged_node;
    merge.to_merge_with->mutations.swap(merge.to_merge_width_unique_mutations);

    return merged_node;
}

void finalize_children(MAT::Node* parent,ConfirmedMove& edits,MAT::Tree* tree,const Original_State_t& checker,std::vector<std::pair<MAT::Node*, MAT::Node*>>& deleted_map){
    if(parent->is_leaf()){
        assert(edits.removed.empty());
        assert(edits.added.front()->identifier==parent->identifier);
        tree->rename_node(parent->identifier, std::to_string(++tree->curr_internal_node));
        MAT::Node* sample_node=edits.added.front();
        sample_node->tree=tree;
        auto result=tree->all_nodes.emplace(sample_node->identifier,sample_node);
        assert(result.second);
        assert(tree->get_node(sample_node->identifier)==sample_node);
    }

    //remove children
    for(MAT::Node* child_to_remove:edits.removed){
        //assert(child_to_remove->parent==parent);
        remove_child(child_to_remove,parent);
    }

    #ifdef DETAIL_DEBUG_CHILD_MERGER_CHECK
    Original_State_t copy(checker);
    Mutation_Set parent_mutations;
    get_mutation_set(parent, parent_mutations);
    #endif 

    //plain add, forgoing merge for testing
    for(MAT::Node* to_add:edits.added){
        to_add->parent=parent;
        parent->children.push_back(to_add);
#ifdef DETAIL_DEBUG_CHILD_MERGER_CHECK
        check_samples_worker(to_add, parent_mutations, copy);
#endif
    }
    
    while(parent->children.empty()){
        tree->all_nodes.erase(parent->identifier);
        std::vector<MAT::Node*>& parent_children=parent->parent->children;
        parent_children.erase(std::find(parent_children.begin(),parent_children.end(),parent));
        deleted_map.emplace_back(parent,parent->parent);
        MAT::Node* next_parent=parent->parent;
        delete parent;
        parent=next_parent;
    }

    /*
    //add children
    Merge_Discriptor merge;
    std::vector<Merge_Discriptor> possible_merges;
    possible_merges.reserve(parent->children.size());
    possible_merges.emplace_back();
    std::vector<MAT::Node*> not_mergeable;
    for (auto iter=edits.added.rbegin(); iter<edits.added.rend(); iter++) {
        bool this_mergable=false;
        for (auto other_added_iter=iter+1;other_added_iter<edits.added.rend(); other_added_iter++) {
            this_mergable|=check_mergable(*iter, *other_added_iter, possible_merges);
        }
        for(MAT::Node* c:parent->children){
            this_mergable|=check_mergable(*iter,c, possible_merges);
        }
        if (!this_mergable) {
            parent->add_child(*iter);
            #ifndef NDEBUG
            check_samples_worker(*iter, parent_mutations, copy);
            #endif
        }
    }
    possible_merges.pop_back();

    if(possible_merges.size()>=1){
        std::unordered_map<MAT::Node*,MAT::Node*> merged_nodes;
        std::unordered_set<MAT::Node*> finished_nodes;
        std::vector<MAT::Node*> merge_conflicted_children;
        std::make_heap(possible_merges.begin(),possible_merges.end(),Merge_Comparator());
        while(!possible_merges.empty()){
            std::pop_heap(possible_merges.begin(),possible_merges.end(),Merge_Comparator());
            Merge_Discriptor& m=possible_merges.back();
            //either nodes are already merged
            if (finished_nodes.count(m.this_child)||finished_nodes.count(m.to_merge_with)) {
                possible_merges.pop_back();
                continue;
            }
            auto iter=merged_nodes.find(m.to_merge_with);
            if(iter==merged_nodes.end()){
                MAT::Node* new_node=execute_merge(m, tree,parent);
#ifndef NDEBUG
                check_samples_worker(new_node, parent_mutations, copy);
#endif
                finished_nodes.insert(m.this_child);
                finished_nodes.insert(m.to_merge_with);
                merged_nodes.emplace(m.to_merge_with,new_node);
                possible_merges.pop_back();
                continue;
            }
            MAT::Node* merged_node=iter->second;
            m.this_child->mutations.set_difference(merged_node->mutations, possible_merges.back().src_unique_mutations, possible_merges.back().to_merge_width_unique_mutations, possible_merges.back().shared_mutations);
            if(possible_merges.back().shared_mutations.size()){
                std::push_heap(possible_merges.begin(),possible_merges.end(),Merge_Comparator());
            }else {
                merge_conflicted_children.push_back(m.this_child);
                possible_merges.pop_back();
            }
        }
        for(MAT::Node* c:merge_conflicted_children){
            if (!finished_nodes.count(c)) {
                parent->add_child(c);
#ifndef NDEBUG
                check_samples_worker(c, parent_mutations, copy);
#endif

            }
        }
    }
    {
        Mutation_Set parent_mutations;
        Sample_Mut_Type copy(checker);
        get_mutation_set(parent->parent, parent_mutations);
        check_samples_worker(parent, parent_mutations, copy);
    }
    */
}