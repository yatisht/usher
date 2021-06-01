#include "mutation_annotated_tree.hpp"
#include "check_samples.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include "tbb/parallel_for_each.h"
#include "Fitch_Sankoff.hpp"
#include <cstdio>
#include "import_vcf.hpp"
#include "apply_move/apply_move.hpp"
namespace MAT = Mutation_Annotated_Tree;

int main(int argc, char** argv){
    MAT::Tree t;
    t.load_detatiled_mutations("with_boundary_mut.pb");
    fprintf(stderr, "%zu\n",t.get_parsimony_score());
    t.condense_leaves();
    Original_State_t origin_state;
    
    printf("%d condensed nodes\n",t.condensed_nodes.size());
}