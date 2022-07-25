#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/usher-sampled/usher.hpp"
#include "src/usher-sampled/mapper.hpp"
#include "src/usher-sampled/place_sample.hpp"
#include <vector>
#include <array>
#define INDEX_END_POSITION INT_MAX
struct index_ele {
    int dfs_idx;
    int dfs_idx_end;
    int children_idx;
};
typedef std::vector<std::array<std::vector<index_ele>, 4>> index_tree;
struct traversal_track_elem {
    int mutation_size;
    int dfs_end_idx;
    const Mutation_Annotated_Tree::Mutation* mutation_start;
};
struct Traversal_Info {
    index_tree indexes;
    std::vector<traversal_track_elem> traversal_track;
    int tree_height;
};
Traversal_Info build_idx(MAT::Tree& tree);
move_type* place_sample_fixed_idx(const Traversal_Info &in,
                                  Sample_Muts* to_search,
                                  const std::vector<MAT::Node*>& dfs_ordered_nodes);