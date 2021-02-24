#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/mman.h>
#include <sys/types.h>
#include <tbb/concurrent_vector.h>
#include <tbb/queuing_rw_mutex.h>
#include <unistd.h>
#define POS_OFFSET 0
#define MUT_NUC_OFFSET 4
#define REF_NUC_OFFSET 5
#define PAR_NUC_OFFSET 6
#define LCA_PARENT_STATE_OFFSET 7
#define META_SIZE 8

#define SRC_OFFSET 0
#define DST_OFFSET 1
#define LCA_OFFSET 2
#define RANGE_FIRST_OFFSET 3
#define RANGE_SECOND_OFFSET 4
#define MOVE_META_SIZE 5
static void serialize_FS_Results(const std::vector<Fitch_Sankoff_Result_Final>& to_serialize,FILE* fd,size_t& offset){
    char temp[META_SIZE];
    size_t n_states=to_serialize.size();
    offset+=fwrite(&n_states, 1,sizeof(n_states),fd);
    size_t dist=to_serialize[0].scores.size();
    size_t ori_offset=offset;
    for(size_t state_idx=0;state_idx<n_states;state_idx++){
        const Fitch_Sankoff_Result_Final& result=to_serialize[state_idx];
        //0-3 bytes position
        *((int*)temp+POS_OFFSET)=result.mutation.position;
        //4-6 byte mut_nuc
        temp[MUT_NUC_OFFSET]=result.mutation.mut_nuc;
        temp[REF_NUC_OFFSET]=result.mutation.is_missing?(0x80|result.mutation.ref_nuc):result.mutation.ref_nuc;
        temp[PAR_NUC_OFFSET]=result.mutation.par_nuc;
        //7 byte LCA
        temp[LCA_PARENT_STATE_OFFSET]=result.LCA_parent_state;
        offset+=fwrite(temp, 1,META_SIZE,fd);
        assert(dist==result.scores.size());
    }
    size_t meta_end=ori_offset+META_SIZE*n_states;
    assert(meta_end==offset);
    for(size_t state_idx=0;state_idx<n_states;state_idx++){
        const Fitch_Sankoff_Result_Final& result=to_serialize[state_idx];
        assert(result.scores.size()==dist);
        assert(offset==meta_end+(dist*sizeof(Fitch_Sankoff::Score_Type)*state_idx));
        offset+=fwrite(result.scores.data(),1, sizeof(Fitch_Sankoff::Score_Type)*dist,fd);
    }
}

static void deserialize_FS_Results(std::vector<Fitch_Sankoff_Result_Deserialized>& out,size_t dist,const char* in,std::string& chrom){
    size_t n_states=*((size_t*) in);
    in+=sizeof(size_t);
    const char* meta_end=in+META_SIZE*n_states;
    out.reserve(n_states);
    for(size_t state_idx=0;state_idx<n_states;state_idx++){
        out.emplace_back();
        //0-3 bytes position
        MAT::Mutation& mut=out.back().mutation;
        mut.position=*((int*)in+POS_OFFSET);
        //4-6 byte mut_nuc
        mut.mut_nuc=in[MUT_NUC_OFFSET];
        mut.is_missing=in[REF_NUC_OFFSET]&0x80;
        mut.ref_nuc=in[REF_NUC_OFFSET]&0x7f;
        mut.par_nuc=in[PAR_NUC_OFFSET];
        mut.chrom=chrom;
        out.back().LCA_parent_state=in[LCA_PARENT_STATE_OFFSET];
        in+=META_SIZE;
        out.back().scores=(Fitch_Sankoff::Score_Type*)(meta_end+(dist*sizeof(Fitch_Sankoff::Score_Type)*state_idx));
    }
}


struct Move_Sorter{
    bool operator()(const std::pair<int, size_t>& first, const std::pair<int, size_t>& second) const{
        return first.first<second.first;
    }
};

static void serialize_profitable_moves(tbb::concurrent_vector<Profitable_Move*>& to_serialize,size_t& offset,std::vector<std::pair<int, size_t>>& file_offsets,FILE* fd,std::string& chrom){
    size_t temp[MOVE_META_SIZE];
    for(const Profitable_Move* ind:to_serialize){
        //save offset
        file_offsets.emplace_back(ind->score_change,offset);
        assert(offset==ftell(fd));
        //save path first
        size_t n_hope=ind->path.size();
        offset+=fwrite( &n_hope, 1,sizeof(size_t),fd);
        offset+=fwrite(ind->path.data(),1, (sizeof(MAT::Node*))*n_hope,fd);
        //save src, dst and LCA
        temp[SRC_OFFSET]=(size_t)ind->src;
        temp[DST_OFFSET]=(size_t)ind->dst;
        temp[LCA_OFFSET]=(size_t)ind->LCA;
        temp[RANGE_FIRST_OFFSET]=ind->range.first;
        temp[RANGE_SECOND_OFFSET]=ind->range.second;
        chrom=ind->states.begin()->mutation.chrom;
        offset+=fwrite(temp, 1,(sizeof(size_t))*MOVE_META_SIZE,fd);
        //save Fitch Sankoff results
        serialize_FS_Results(ind->states,fd,offset);
        delete ind;
    }
    to_serialize.clear();
}

static void deserialize_profitable_moves(char* start,Profitable_Move_Deserialized* out,std::string& chrom){
    //assuming seeked pass path
    MAT::Node** recasted=(MAT::Node**)start;
    out->src=recasted[SRC_OFFSET];
    out->dst=recasted[DST_OFFSET];
    out->LCA=recasted[LCA_OFFSET];
    out->range.first=(size_t)recasted[RANGE_FIRST_OFFSET];
    out->range.second=(size_t)recasted[RANGE_SECOND_OFFSET];
    deserialize_FS_Results(out->states,out->range.second-out->range.first, start+MOVE_META_SIZE*sizeof(size_t),chrom);
}

size_t Profitable_Moves_Cacher::get_path(MAT::Node*** out){
    *out=(MAT::Node**)(mapped_address+file_offsets_iter->second+sizeof(size_t));
    return *((size_t*)(mapped_address+file_offsets_iter->second));
}

Profitable_Move_Deserialized* Profitable_Moves_Cacher::operator*(){
    MAT::Node** start;
    size_t path_len=get_path(&start);
    start+=path_len;
    Profitable_Move_Deserialized* out=new Profitable_Move_Deserialized;
    deserialize_profitable_moves((char*)start, out,chrom);
    return out;
}

void Profitable_Moves_Cacher::operator()(){
    size_t offset=fseek(fd, 0, SEEK_SET);
    tbb::concurrent_vector<Profitable_Move*> to_swap;
    while (true) {
        std::unique_lock<std::mutex> finished_lock(finish_mutex);
        finish_cv.wait_for(finished_lock,std::chrono::seconds(1));
        if(finished){
            finished_lock.unlock();
            break;
        }
        finished_lock.unlock();
        if (to_monitor.size()>20) {
            {
                tbb::queuing_rw_mutex::scoped_lock lock(swap_lock,true);
                to_swap.swap(to_monitor);
                
                lock.release();
            }
            serialize_profitable_moves(to_swap, offset, file_offsets, fd,chrom);
        }
    }
    serialize_profitable_moves(to_monitor, offset, file_offsets, fd,chrom);
    length=offset;
    std::sort(file_offsets.begin(),file_offsets.end(),Move_Sorter());
    file_offsets_iter=file_offsets.begin();
    fflush(fd);
    mapped_address=(char*)mmap(nullptr, offset, PROT_READ, MAP_SHARED, raw_fd, 0);
    fclose(fd);
}

void Profitable_Moves_Cacher::run(){
        this_thread=std::thread(std::ref(*this));
}
Profitable_Moves_Cacher::Profitable_Moves_Cacher(tbb::concurrent_vector<Profitable_Move*>& to_monitor,tbb::queuing_rw_mutex& rw_mutex):to_monitor(to_monitor),swap_lock(rw_mutex),finished(false){
        filename=(char*)malloc(15);
        strcpy(filename, "aXXXXXXXXX");
        raw_fd=mkstemp(filename);
        fd=fdopen(raw_fd,"w");
    }
void Profitable_Moves_Cacher::finish(){
    {
        std::lock_guard<std::mutex> finish_lock(finish_mutex);
        finished=true;
    }
    finish_cv.notify_all();
    this_thread.join();
}
Profitable_Moves_Cacher::~Profitable_Moves_Cacher(){
        munmap(mapped_address, length);
        unlink(filename);
        free(filename);
}