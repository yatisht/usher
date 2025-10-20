#include "ripples.hpp"
size_t check_parallelizable(const MAT::Node *root,
                            std::vector<bool> &do_parallel,
                            size_t parallel_threshold, size_t check_threshold,
                            unsigned short &tree_height,
                            std::vector<Mapper_Info> &traversal_track,
                            unsigned short level) {
    size_t child_counted_size = 0;
    if ((root->dfs_end_idx - root->dfs_idx) >= check_threshold) {
        child_counted_size++;
        auto cur_idx = traversal_track.size();
        traversal_track.push_back(
            Mapper_Info{root->mutations.data(),
                        root->mutations.data() + root->mutations.size(),
                        UINT32_MAX, level, root->children.empty()});
        for (const auto child : root->children) {
            child_counted_size += check_parallelizable(
                child, do_parallel, parallel_threshold, check_threshold,
                tree_height, traversal_track, level + 1);
        }
        if (child_counted_size > parallel_threshold) {
            do_parallel[root->dfs_idx] = true;
        }
        traversal_track[cur_idx].sibling_start_idx = traversal_track.size();
        tree_height = std::max(tree_height, level);
    }
    return child_counted_size;
}
