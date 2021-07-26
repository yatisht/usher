#define LOAD
#include "../check_samples.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "transpose_vcf.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <unordered_map>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;
struct Condesed_Muts {
    std::vector<MAT::Mutation> not_Ns;
    std::vector<MAT::Mutation> Ns;
};
void check_consistent(Mutation_Set *old_muts, Mutation_Set &new_muts,
                      const std::string *sample) {
    for (const auto &mut : *old_muts) {
        if (!new_muts.count(mut)) {
            fprintf(
                stderr,
                "Inconsistent allele at %d of %s, %c in protobuf, %c in VCF\n",
                mut.position, sample->c_str(),
                MAT::get_nuc(mut.get_all_major_allele() & 0xf),
                MAT::get_nuc(MAT::Mutation::refs[mut.position]));
        }
    }
    for (auto &mut : new_muts) {
        auto iter = old_muts->find(mut);
        if (iter == old_muts->end()) {
            if (!(mut.boundary1_all_major_allele &
                  MAT::Mutation::refs[mut.position])) {
                fprintf(stderr,
                        "Inconsistent allele at %d of %s, %c in protobuf, %c "
                        "in VCF\n",
                        mut.position, sample->c_str(),
                        MAT::get_nuc(MAT::Mutation::refs[mut.position]),
                        MAT::get_nuc(mut.get_all_major_allele() & 0xf));
            }
        } else {
            if (!(iter->boundary1_all_major_allele &
                  mut.boundary1_all_major_allele)) {
                if (!(mut.boundary1_all_major_allele &
                      MAT::Mutation::refs[mut.position])) {
                    fprintf(stderr,
                            "Inconsistent allele at %d of %s, %c in protobuf, "
                            "%c in VCF\n",
                            mut.position, sample->c_str(),
                            MAT::get_nuc(iter->get_all_major_allele()),
                            MAT::get_nuc(mut.get_all_major_allele() & 0xf));
                }
            }
        }
    }
}
struct Adder {
    Condesed_Muts *condensed_add;
    Mutation_Set *uncondensed_add;
    Mutation_Set staging;
    std::unordered_set<int> set_positions;
    const std::string *sample;
    Adder(Condesed_Muts *condensed_add)
        : condensed_add(condensed_add), uncondensed_add(nullptr) {}
    Adder(Mutation_Set *uncondensed_add,const std::string* sample)
        : condensed_add(nullptr), uncondensed_add(uncondensed_add),sample(sample) {}
    Adder() : condensed_add(nullptr), uncondensed_add(nullptr) {}
    void add_mut(MAT::Mutation &mutation,
                 std::vector<MAT::Mutation> &condensed_out) {
        if (mutation.position >= (int) MAT::Mutation::refs.size() ||
            (!MAT::Mutation::refs[mutation.position])) {
            return;
        }
        if (condensed_add) {
            condensed_out.push_back(mutation);
        } else if (uncondensed_add) {
            // set_positions.insert(mutation.position);
            staging.insert(mutation);
            // insert_checking_consistent(uncondensed_add, mutation);
        }
    }
    void add_Not_N(int position, char mut) {
        MAT::Mutation mutation;
        mutation.position = position;
        mutation.boundary1_all_major_allele = mut;
        add_mut(mutation, condensed_add->not_Ns);
    }
    void add_N(int first, int second) {
        for (int pos = first; pos <= second; pos++) {
            MAT::Mutation mutation;
            mutation.position = pos;
            mutation.boundary1_all_major_allele = 0xf;
            add_mut(mutation, condensed_add->Ns);
        }
    }
    ~Adder() {
        if (uncondensed_add) {
            check_consistent(uncondensed_add, staging, sample);
            uncondensed_add->swap(staging);
        }
    }
};

struct Sample_Adder {
    const Original_State_t &not_condensed;
    const std::unordered_map<std::string, Condesed_Muts *> &condensed;
    Adder set_name(std::string &&name) {
        auto not_condensed_iter = not_condensed.find(name);
        if (not_condensed_iter != not_condensed.end()) {
            return Adder{&(not_condensed_iter->second),&(not_condensed_iter->first)};
        }
        auto condensed_iter = condensed.find(name);
        if (condensed_iter == condensed.end()) {
            //fprintf(stderr, "Sample %s not found in tree\n", name.c_str());
            return Adder();
        }
        return Adder(condensed_iter->second);
    }
};
struct condesed_container {
    Mutation_Set *mut_set;
    std::vector<Condesed_Muts *> children;
    const std::string* sample;
    condesed_container(
        std::unordered_map<std::string, Condesed_Muts *> &out,
        Original_State_t *ori_state,
        const std::pair<std::string, std::vector<std::string>> &samp_name) {
        sample=&samp_name.first;
        mut_set = &(*ori_state)[samp_name.first];
        for (const auto &condensed_children : samp_name.second) {
            children.push_back(new Condesed_Muts);
            out.emplace(condensed_children, children.back());
        }
    }
};
struct iter_heap {
    struct Mut_Iter {
        std::vector<MAT::Mutation>::const_iterator iter;
        std::vector<MAT::Mutation>::const_iterator end;
        Mut_Iter(const std::vector<MAT::Mutation> &muts) {
            iter = muts.begin();
            end = muts.end();
        }
        bool operator!() const { return iter == end; }
        operator const MAT::Mutation &() const { return *iter; }
        void operator++() { iter++; }
        bool operator==(const MAT::Mutation &mut) const {
            return iter->position == mut.position;
        }
        bool operator<(const Mut_Iter &other) const {
            return iter->position > other.iter->position;
        }
        uint8_t get_allele() const { return iter->get_all_major_allele(); }
    };
    std::vector<Mut_Iter> iters;
    size_t size;
    int lst;
    iter_heap(std::vector<std::vector<MAT::Mutation>> &to_merge) {
        lst = 0;
        size = to_merge.size();
        iters.reserve(to_merge.size());
        for (const auto &one : to_merge) {
            iters.emplace_back(one);
        }
        std::make_heap(iters.begin(), iters.end());
        std::pop_heap(iters.begin(), iters.end());
    }
    void get_one(Mutation_Set *mut_set) {
        MAT::Mutation mut = iters.back();
        assert(mut.position > lst);
        lst = mut.position;
        size_t count = 0;
        while (iters.back() == mut) {
            auto &inc = iters.back();
            count++;
            mut.boundary1_all_major_allele &= inc.get_allele();
            ++inc;
            if (!inc) {
                iters.pop_back();
                if (iters.empty()) {
                    break;
                }
            } else {
                std::push_heap(iters.begin(), iters.end());
            }
            std::pop_heap(iters.begin(), iters.end());
        }
        assert(count <= size);
        if (count == size && mut.boundary1_all_major_allele !=
                                 MAT::Mutation::refs[mut.position]) {
            mut_set->insert(mut);
        }
    }
    operator bool() const { return !iters.empty(); }
};
void add_ambuiguous_mutations(const char *path, Original_State_t &to_patch,
                              Mutation_Annotated_Tree::Tree &tree) {
    std::vector<condesed_container> condensed_children;
    condensed_children.reserve(tree.condensed_nodes.size());
    std::unordered_map<std::string, Condesed_Muts *> condensed_map;
    for (const auto &condensed_one : tree.condensed_nodes) {
        condensed_children.emplace_back(condensed_map, &to_patch,
                                        condensed_one);
    }
    Sample_Adder adder{to_patch, condensed_map};
    load_mutations(path, num_threads, adder);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, condensed_children.size()),
        [&condensed_children](const tbb::blocked_range<size_t> &range) {
            for (size_t idx = range.begin(); idx < range.end(); idx++) {
                std::vector<std::vector<MAT::Mutation>> to_merge;
                for (auto &c : condensed_children[idx].children) {
                    if (!c) {
                        return;
                    }
                    if (c->not_Ns.empty() && c->Ns.empty()) {
                        return;
                    }
                    std::vector<MAT::Mutation> this_merged;
                    std::merge(c->not_Ns.begin(), c->not_Ns.end(),
                               c->Ns.begin(), c->Ns.end(),
                               std::back_inserter(this_merged));
                    to_merge.push_back(std::move(this_merged));
                }
                iter_heap heap(to_merge);
                Mutation_Set mut_set_saging;
                while (heap) {
                    heap.get_one(&mut_set_saging);
                }
                check_consistent(condensed_children[idx].mut_set,
                                 mut_set_saging,condensed_children[idx].sample);
                condensed_children[idx].mut_set->swap(mut_set_saging);
            }
        });
}
