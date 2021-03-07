#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <tbb/blocked_range.h>
#include "check_samples.hpp"
#include <tbb/flow_graph.h>
#include <utility>
#include <vector>
#include "Twice_Bloom_Filter.hpp"
#include "src/tree_rearrangement.hpp"
tbb::concurrent_vector<MAT::Node*> postponed;
extern uint32_t num_cores;
static void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output){
    Twice_Bloom_Filter filter;
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            filter.insert(m.position);
        }
    }
    for(MAT::Node* n:all_nodes){
        for(const MAT::Mutation& m:n->mutations){
            if(filter.query(m.position)){
                output.push_back(n);
                break;
            }
        }
    }
}
static void fix_condensed_nodes(MAT::Tree* tree){
    std::vector<MAT::Node*> nodes_to_fix;
    for(auto iter:tree->all_nodes){
        if (tree->condensed_nodes.count(iter.first)&&(!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for(auto node:nodes_to_fix){
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier, std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier,node);
    }
}
static void push_all_nodes(MAT::Tree* tree,std::vector<MAT::Node*>& nodes){
    nodes.clear();
    for(auto a:tree->all_nodes){
        nodes.push_back(a.second);
    }
}
static void feed_nodes(std::vector<MAT::Node *> &to_feed,
                         std::deque<MAT::Node*>& out_queue,
                         std::mutex& out_mutex,
                         std::condition_variable& out_pushed,
                         const std::vector<MAT::Node *> &dfs_ordered_nodes) {
    std::vector<size_t> this_round_idx;
    std::vector<size_t> next_round_idx;

#define output_node(this_node) \
{std::lock_guard<std::mutex> lk(out_mutex);\
out_queue.push_back(this_node);\
out_pushed.notify_one();}\
this_node=this_node->parent;\
if(this_node->parent){\
    next_round_idx.push_back(this_node->index);\
}

    {
        this_round_idx.reserve(to_feed.size());
        for (auto node : to_feed) {
            auto parent = node->parent;
            if (parent)
                this_round_idx.push_back(node->index);
        }
    }
    while (!this_round_idx.empty()) {
    std::sort(this_round_idx.begin(), this_round_idx.end());
    for (auto iter = this_round_idx.begin(); iter < this_round_idx.end()-1; iter++) {
        auto next_ele=*(iter+1);
        if(next_ele==*iter){
            continue;
        }
        if (!check_grand_parent(dfs_ordered_nodes[next_ele],dfs_ordered_nodes[*iter])) {
            MAT::Node* this_node=dfs_ordered_nodes[*iter];
            output_node(this_node);
        } else {
            next_round_idx.push_back(*iter);
        }
    }
    MAT::Node* last_node=dfs_ordered_nodes[this_round_idx.back()];
    output_node(last_node);
    this_round_idx.swap(next_round_idx);
    next_round_idx.clear();
    }

    {
        std::lock_guard<std::mutex> lk(out_mutex);
        out_queue.push_back(nullptr);
        out_pushed.notify_one();
    }
}

typedef tbb::flow::multifunction_node<char,std::tuple<MAT::Node*>> Seach_Source_Node_t;
struct Search_Source{
    std::deque<MAT::Node*>& buffer;
    std::mutex& queue_mutex;
    std::condition_variable& ready;
    unsigned int conflict_pct_limit;
    unsigned int min_batch_size;
    mutable unsigned int conflicting_count;
    mutable unsigned int inflight_left;
    unsigned int inflight_limit;
    mutable size_t all_count;
    std::atomic_bool& all_nodes_done;
    bool& graph_reseted;
    Search_Source(std::deque<MAT::Node*>& buff, std::mutex& queue_mutex, std::condition_variable& ready,unsigned int inflight_limit,unsigned int conflict_pct_limit,unsigned int min_batch_size,std::atomic_bool& all_nodes_done,bool& graph_reseted):buffer(buff),queue_mutex(queue_mutex),ready(ready),conflict_pct_limit(conflict_pct_limit),min_batch_size(min_batch_size),conflicting_count(0),inflight_left(inflight_limit),inflight_limit(inflight_limit),all_count(0),all_nodes_done(all_nodes_done),graph_reseted(graph_reseted){}
    void reset()const{
        all_count=0;
        inflight_left=inflight_limit;
        conflicting_count=0;
        graph_reseted=false;
    }
    void operator()(char in,Seach_Source_Node_t::output_ports_type& out) const{
        inflight_left++;
        if(graph_reseted){
            reset();
        }
        //fprintf(stderr,"%d conflicting moves over %zu moves so far\n",conflicting_count,all_count);
        
        if(in&MOVE_FOUND_MASK){
        all_count++;
        if(!(in&NONE_CONFLICT_MOVE_MASK)){
            conflicting_count++;
        }        
        if ((conflicting_count>min_batch_size)&&(all_count*conflict_pct_limit<conflicting_count*100)) {
           return;
        }
        }
        if(all_nodes_done){
           //fprintf(stderr,"all_nodes_done");
            return;
        }
        while(inflight_left>0){
            std::unique_lock<std::mutex> lk(queue_mutex);
            while (buffer.empty()) {
                ready.wait(lk);
            }
            if (buffer.front()==nullptr) {
                buffer.pop_front();
                all_nodes_done.store(true);
           	fprintf(stderr,"all_nodes_done");
                return;
            }else {
                std::get<0>(out).try_put(buffer.front());
                buffer.pop_front();
                inflight_left--;
            }
        }
    }
};

void Tree_Rearrangement::refine_trees(std::vector<MAT::Tree> &optimal_trees,int radius) {

    for (MAT::Tree& this_tree : optimal_trees) {
        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        Original_State_t ori;
        std::vector<MAT::Node*>& to_optimize=this_tree.new_nodes;
        check_samples(this_tree.root, ori,&this_tree);
        Pending_Moves_t pending_moves;
        tbb::concurrent_vector<Profitable_Move*> profitable_moves;
        //building pipeline
        tbb::queuing_rw_mutex mutex;
        tbb::flow::graph search_graph;
        std::deque<MAT::Node*> search_queue;
        std::mutex queue_mutex;
        std::condition_variable queue_ready;
        std::atomic_bool all_nodes_done;
        bool graph_reseted;
        Search_Source input_functor(search_queue,queue_mutex,queue_ready,num_cores,50,800,all_nodes_done,graph_reseted);
        Seach_Source_Node_t input(search_graph,1,input_functor);

        tbb::flow::function_node<MAT::Node*,Possible_Moves*> neighbors_finder(search_graph,tbb::flow::unlimited,Neighbors_Finder(radius));
        tbb::flow::make_edge(std::get<0>(input.output_ports()),neighbors_finder);

        tbb::flow::function_node<Possible_Moves*,Candidate_Moves*> parsimony_score_calculator(search_graph,tbb::flow::unlimited,Parsimony_Score_Calculator{ori,dfs_ordered_nodes});
        tbb::flow::make_edge(neighbors_finder,parsimony_score_calculator);
        
        tbb::flow::function_node<Candidate_Moves*,Profitable_Moves_From_One_Source*> profitable_move_enumerator(search_graph,tbb::flow::unlimited,Profitable_Moves_Enumerator{dfs_ordered_nodes,ori});
        tbb::flow::make_edge(parsimony_score_calculator,profitable_move_enumerator);

        std::vector<Profitable_Move*> non_conflicting_moves;
        Cross_t potential_crosses;
        Mut_t repeatedly_mutation_loci;
        tbb::flow::function_node<Profitable_Moves_From_One_Source*,char> conflict_resolver(search_graph,1,Conflict_Resolver{non_conflicting_moves,potential_crosses,repeatedly_mutation_loci,to_optimize});
        tbb::flow::make_edge(profitable_move_enumerator,conflict_resolver);
        tbb::flow::make_edge(conflict_resolver,input);

        bool have_improvement=true;
        while (have_improvement) {
            have_improvement=false;
            find_nodes_with_recurrent_mutations(dfs_ordered_nodes, to_optimize);
	    fprintf(stderr,"next_round\n");
            //push_all_nodes(&this_tree, to_optimize);
        while (!to_optimize.empty()) {
	    feed_nodes(to_optimize,search_queue,queue_mutex,queue_ready,dfs_ordered_nodes);
            all_nodes_done.store(false);
            while (!all_nodes_done.load()) {
            graph_reseted=true;
            potential_crosses.clear();
            repeatedly_mutation_loci.clear();
	    non_conflicting_moves.clear();

            input.try_put(0);
            search_graph.wait_for_all();
            to_optimize.clear();
	    fprintf(stderr,"%zu moves profitable\n",non_conflicting_moves.size());
            if(!non_conflicting_moves.empty()){
                Pending_Moves_t tree_edits;
                // tbb::parallel_for(tbb::blocked_range<size_t>(0,
                // non_conflicting_moves.size()),
                // Move_Executor{dfs_ordered_nodes,this_tree,non_conflicting_moves,tree_edits});
                Move_Executor temp{dfs_ordered_nodes, this_tree,
                                   non_conflicting_moves, tree_edits, ori};
                tbb::blocked_range<size_t> range_temp(
                    0, non_conflicting_moves.size());
                temp(range_temp);
                non_conflicting_moves.clear();
                // this_tree.finalize();
                /*
                tbb::parallel_for_each(
                    tree_edits.begin(), tree_edits.end(),
                    [&this_tree,
                     &ori](const std::pair<MAT::Node *, ConfirmedMove> &in) {
                   });
                */
                std::vector<std::pair<MAT::Node*, MAT::Node*>> deleted_map;
                for (Pending_Moves_t::const_iterator iter=tree_edits.begin(); iter!=tree_edits.end(); iter++) {
                    finalize_children(reinterpret_cast<MAT::Node *>(iter->first),
                                      const_cast<ConfirmedMove &>(iter->second),
                                      &this_tree, ori,deleted_map);
                }
                for(MAT::Node*& next_node:to_optimize){
                    for(const auto& deleted:deleted_map){
                        if(next_node==deleted.first){
                            next_node=deleted.second;
                        }
                    }
                }
                for(MAT::Node*& next_node:search_queue){
                    for(const auto& deleted:deleted_map){
                        if(next_node==deleted.first){
                            next_node=deleted.second;
                        }
                    }
                }
 
                dfs_ordered_nodes = this_tree.depth_first_expansion();
//#ifndef NDEBUG
                //Original_State_t copy(ori);
                //check_samples(this_tree.root, copy,&this_tree);
                fprintf(stderr, "Between refinement: %zu \n",this_tree.get_parsimony_score());
//#endif
                have_improvement=true;
            }
            }
        }
        }

        check_samples(this_tree.root, ori,&this_tree);
        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        fix_condensed_nodes(&this_tree);
        this_tree.reassign_level();
    }
}
