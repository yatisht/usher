#include "place_sample.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <climits>
#include <csignal>
#include <cstdio>
#include <tbb/pipeline.h>
#include <thread>
#include <mpi.h>
static void check_parent(MAT::Node* root,MAT::Tree& tree){
    if (root!=tree.get_node(root->node_id)) {
        fprintf(stderr, "dict mismatch\n");
        raise(SIGTRAP);
    }
    for (auto child : root->children) {
        if (child->parent!=root) {
            fprintf(stderr, "parent mismatch\n");
            raise(SIGTRAP);
        }
    }
}

static std::string
serialize_move(move_type *in) {
    Mutation_Detailed::search_result result;
    result.set_sample_id(std::get<1>(*in));
    result.mutable_place_targets()->Reserve(std::get<0>(*in).size());
    for (const auto &place_target : std::get<0>(*in)) {
        auto new_target = result.add_place_targets();
        new_target->set_target_node_id(place_target.target_node->node_id);
        new_target->set_parent_node_id(place_target.parent_node->node_id);
        fill_mutation_vect(new_target->mutable_sample_mutation_positions(),
                           new_target->mutable_sample_mutation_other_fields(),
                           place_target.sample_mutations);
        fill_mutation_vect(new_target->mutable_split_mutation_positions(),
                           new_target->mutable_split_mutation_other_fields(),
                           place_target.splited_mutations);
        fill_mutation_vect(new_target->mutable_shared_mutation_positions(),
                           new_target->mutable_shared_mutation_other_fields(),
                           place_target.shared_mutations);
        if (place_target.sample_mutations.size()==0) {
            fprintf(stderr, "follower send size mismatch\n");
            raise(SIGTRAP);
        }
        if (new_target->sample_mutation_positions_size()!=place_target.sample_mutations.size()) {
            fprintf(stderr, "follower send size mismatch\n");
            raise(SIGTRAP);
        }
    }
    return result.SerializeAsString();
}
struct Send_Main_Tree_Target {
    void
    operator()(move_type *in)const {
        auto buffer = serialize_move(in);
        mpi_trace_print( "follower sent placement \n");
        MPI_Send(buffer.c_str(), buffer.size(), MPI_BYTE, 0, PROPOSED_PLACE,
                 MPI_COMM_WORLD);
        delete in;
    }
};
void check_order_node(MAT::Node* node){
    int prev_pos=-1;
    for (const auto& mut : node->mutations) {
            if (mut.get_position()>=MAT::Mutation::refs.size()) {
                fprintf(stderr, "%zu: placement_check strange size\n",node->node_id);
                raise(SIGTRAP);
            }
        if (mut.get_position()<=prev_pos) {
            fprintf(stderr, "%zu:placement out of order\n",node->node_id);
            raise(SIGTRAP);
        }
        prev_pos=mut.get_position();
    }
}
void check_order(MAT::Mutations_Collection& in){
    int prev_pos=-1;
    for (const auto& mut : in) {
        if (mut.get_position()<=prev_pos) {
            if (mut.get_position()>=MAT::Mutation::refs.size()) {
                fprintf(stderr, "placement_check strange size\n");
                raise(SIGTRAP);
            }
            fprintf(stderr, "placement out of order\n");
            raise(SIGTRAP);
        }
        prev_pos=mut.get_position();
    }
}
static void recv_and_place_follower(MAT::Tree &tree,
                           std::vector<MAT::Node *> &deleted_nodes) {
    fprintf(stderr, "Place Recievier started \n");
    while (true) {
        size_t bcast_size;
        MPI_Bcast(&bcast_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        mpi_trace_print("Recieving move of size %zu \n",bcast_size);
        if (bcast_size == 0) {
            return;
        }
        auto buffer = new uint8_t[bcast_size];
        MPI_Bcast(buffer, bcast_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        Mutation_Detailed::placed_target parsed_target;
        parsed_target.ParseFromArray(buffer, bcast_size);
        //fprintf(stderr,"Recieved move target :%zu, split %zu, sample: %zu\n",parsed_target.target_node_id(),parsed_target.split_node_id(),parsed_target.sample_id());
        MAT::Mutations_Collection sample_mutations;
        MAT::Mutations_Collection shared_mutations;
        MAT::Mutations_Collection splitted_mutations;
        load_mutations(parsed_target.sample_mutation_positions(),
                       parsed_target.sample_mutation_other_fields(),
                       sample_mutations.mutations);
        load_mutations(parsed_target.split_mutation_positions(),
                       parsed_target.split_mutation_other_fields(),
                       splitted_mutations.mutations);
        load_mutations(parsed_target.shared_mutation_positions(),
                       parsed_target.shared_mutation_other_fields(),
                       shared_mutations.mutations);
        /*check_order(sample_mutations);
        check_order(splitted_mutations);
        check_order(shared_mutations);*/
        auto target_node=tree.get_node(parsed_target.target_node_id());
        if (target_node==nullptr) {
            fprintf(stderr, "Node not found %zu \n",parsed_target.target_node_id());
            std::raise(SIGTRAP);
        }
        auto out = update_main_tree(
            sample_mutations, splitted_mutations, shared_mutations,
            target_node,
            parsed_target.sample_id(), tree, parsed_target.split_node_id());
        target_node=tree.get_node(parsed_target.target_node_id());
        check_parent(tree.root, tree);
        /*check_order(target_node->mutations);
        check_order(tree.get_node(parsed_target.sample_id())->mutations);
        if (out.splitted_node) {
            check_order(out.splitted_node->mutations);
        }
        auto dfs=tree.depth_first_expansion();
        for (auto node : dfs) {
            check_order_node(node);
        }*/
        if (out.deleted_nodes) {
            deleted_nodes.push_back(out.deleted_nodes);
        }
        delete[] buffer;
    }
}


struct Fetcher {
    Sample_Muts* operator()(tbb::flow_control& fc) const {
        int res_size;
        mpi_trace_print("follower send work req \n");
        MPI_Send(&res_size, 0, MPI_BYTE, 0, WORK_REQ_TAG, MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Probe(0, WORK_RES_TAG, MPI_COMM_WORLD, &status);
        mpi_trace_print("follower recieve work res \n");
        MPI_Get_count(&status, MPI_BYTE, &res_size);
        if (res_size==0) {
            fc.stop();
            return nullptr;
        }
        auto buffer=new char[res_size];
        MPI_Recv(buffer, res_size, MPI_BYTE, 0, WORK_RES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        mpi_trace_print("follower recieved work res \n");
        Mutation_Detailed::sample_to_place parsed;
        parsed.ParseFromArray(buffer, res_size);
        auto out=new Sample_Muts;
        out->sample_idx=parsed.sample_id();
        load_mutations(parsed.sample_mutation_positions(),parsed.sample_mutation_other_fields(),out->muts);
        mpi_trace_print( "follower finished parsing \n");
        delete[] buffer;
        return out;
    }
};

void follower_place_sample(MAT::Tree &main_tree,int batch_size){
    std::vector<MAT::Node *> deleted_nodes;
    std::thread tree_update_thread(recv_and_place_follower,std::ref(main_tree),std::ref(deleted_nodes));
    MPI_Barrier(MPI_COMM_WORLD);
    tbb::parallel_pipeline(batch_size,
        tbb::make_filter<void,Sample_Muts*>(tbb::filter::serial,Fetcher())&
        tbb::make_filter<Sample_Muts*,move_type *>(
            tbb::filter::parallel,Finder{main_tree,true})&
        tbb::make_filter<move_type *,void>(
            tbb::filter::serial,Send_Main_Tree_Target())
        );
    for (auto node : deleted_nodes) {
        delete node;
    }
    int ignored;
    MPI_Send(&ignored, 0, MPI_BYTE, 0, PROPOSED_PLACE, MPI_COMM_WORLD);
    tree_update_thread.join();
}
