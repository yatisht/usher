#include "../src/tree_rearrangements/tree_rearrangement_internal.hpp"
#include "src/mutation_annotated_tree.hpp"
#include <bits/types/FILE.h>
#include <tbb/concurrent_vector.h>
#include <tbb/queuing_rw_mutex.h>
#include <vector>
#include <random>

int main(int argc, char** argv){
    tbb::concurrent_vector<Profitable_Move*> test_in;
    std::vector<Profitable_Move> to_match(100);
    std::default_random_engine rng;
    std::uniform_int_distribution<int> size_dist(2,22);
    std::uniform_int_distribution<int> nuc_dist(0,3);
    tbb::queuing_rw_mutex mutex;
    Profitable_Moves_Cacher cacher(test_in,mutex);
    cacher.run();
    assert(to_match.size()==100);
    //generate test cases
    for (int score_change=0; score_change>-100; score_change--) {
        Profitable_Move* temp=&to_match[99+score_change];
        size_t path_size=size_dist(rng);
        for (size_t i=0; i<path_size; i++) {
            temp->path.push_back((MAT::Node*)rng());
        }
        temp->src=(MAT::Node*)rng();
        temp->dst=(MAT::Node*)rng();
        temp->LCA=(MAT::Node*)rng();
        temp->score_change=score_change;
        size_t states_size=size_dist(rng);
        size_t subtree_size=size_dist(rng);
        size_t node_start=rng();
        temp->range.first=node_start;
        temp->range.second=node_start+subtree_size;
        for (size_t state_idx=0; state_idx<states_size; state_idx++) {
            Fitch_Sankoff_Result_Final state;
            state.LCA_parent_state=1<<(nuc_dist(rng));
            state.mutation.is_missing=1&rng();
            state.mutation.ref_nuc=1<<(nuc_dist(rng));
            state.mutation.par_nuc=1<<(nuc_dist(rng));
            state.mutation.mut_nuc=1<<(nuc_dist(rng));
            state.mutation.position=rng();
            state.mutation.chrom="1";
            for (size_t node_idx=0; node_idx<subtree_size; node_idx++) {
                state.scores.emplace_back((MAT::Node*)rng());
                for (int i=0; i<4; i++) {
                    state.scores.back()[i]=1<<(nuc_dist(rng));
                }
            }
            temp->states.push_back(state);
        }
        {
        tbb::queuing_rw_mutex::scoped_lock lock(mutex,false);
        test_in.push_back(new Profitable_Move(*temp));
        }
    }
    cacher.finish();
    //Test output
    for(const Profitable_Move& gold:to_match){
        assert(!cacher.eof());
        MAT::Node** path;
        size_t n_nodes=cacher.get_path(&path);
        assert(gold.path.size()==n_nodes);
        for (size_t i=0; i<gold.path.size(); i++) {
            assert(gold.path[i]==path[i]);
        }
        Profitable_Move_Deserialized* move=*cacher;
        assert(gold.src==move->src);
        assert(gold.dst==move->dst);
        assert(gold.LCA==move->LCA);
        assert(gold.range==move->range);
        assert(gold.states.size()==move->states.size());
        size_t dist=move->range.second-move->range.first;
        for (size_t state_idx=0; state_idx<gold.states.size(); state_idx++) {
            assert(gold.states[state_idx].mutation==move->states[state_idx].mutation);
            assert(gold.states[state_idx].LCA_parent_state==move->states[state_idx].LCA_parent_state);
            for (size_t node_idx=0; node_idx<dist; node_idx++) {
                assert(gold.states[state_idx].scores[node_idx].score==move->states[state_idx].scores[node_idx].score);
            }
        }
        ++cacher;
        delete move;
    }
    assert(cacher.eof());
}