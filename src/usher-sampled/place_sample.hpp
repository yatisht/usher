#include "mutation_detailed.pb.h"
#include "mapper.hpp"
#define DELETED_NODE_THRESHOLD 1000
#define PROPOSED_PLACE 0
#define CONFIRMED_PLACE 1
#define WORK_REQ_TAG 2
#define WORK_RES_TAG 3
struct update_main_tree_output{
    MAT::Node* splitted_node;
    MAT::Node* deleted_nodes;
};
update_main_tree_output update_main_tree(const MAT::Mutations_Collection& sample_mutations,
                                    const MAT::Mutations_Collection& splitted_mutations,
                                    const MAT::Mutations_Collection& shared_mutations,
                                    MAT::Node* target_node,
                                    size_t node_idx, MAT::Tree& tree,size_t split_node_idx) ;
struct Finder{
    MAT::Tree& tree;
    bool do_free;
    std::tuple<std::vector<Main_Tree_Target>,size_t>* operator()(Sample_Muts* in)const {
        auto output=new std::tuple<std::vector<Main_Tree_Target>,size_t>;
        std::get<1>(*output)= in->sample_idx;
        const auto& condensed_muts =in->muts;
        auto main_tree_out=place_main_tree(condensed_muts, tree);
        std::get<0>(*output)=std::move(std::get<0>(main_tree_out));
        if (do_free) {
            delete in;
        }
        return output;
    }    
};
template <typename pos_field_type, typename other_field_type, typename mut_type>
static void fill_mutation_vect(pos_field_type *pos_field,
                               other_field_type *other_field,
                               const mut_type &to_ser) {
    pos_field->Reserve(to_ser.size());
    other_field->Reserve(to_ser.size());
    for (const auto &mut : to_ser) {
        pos_field->Add(*((uint32_t *)(&mut)));
        other_field->Add(*((uint32_t *)(&mut) + 1));
    }
}
template <typename pos_field_type, typename other_field_type, typename MUT_TYPE>
static void load_mutations(const pos_field_type &pos_fields,
                           const other_field_type &other_fields,
                           std::vector<MUT_TYPE> &out) {
    auto mut_count = pos_fields.size();
    out.reserve(mut_count);
    for (int idx = 0; idx < mut_count; idx++) {
        MUT_TYPE mut;
        *((uint32_t *)(&mut)) = pos_fields.Get(idx);
        *((uint32_t *)(&mut) + 1) = other_fields.Get(idx);
    }
}