#include <cstdint>
#include <tbb/parallel_pipeline.h>
#include <unordered_map>
#include <zlib.h>
#include <cassert>
#include <string>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>
//#include <vector>
#define MAX_SIZ 0x30000

struct partitioner {
    const uint8_t*& last_out;
    const uint8_t* const end;
    const uint8_t* operator()(tbb::flow_control& fc)const {
        if (last_out>=end) {
            fc.stop();
        }
        //assert(last_out<end);
        auto out=last_out;
        //fprintf(stderr, "%d\n",*(int*)out);
        unsigned int item_len=*(int*) out;
        last_out+=(item_len+4);
        return out;
    }
};
static unsigned int loadVariant(const uint8_t*& in) {
    unsigned int out=(*in&0x7f);
    int shamt=7;
    while (*in&0x80) {
        in++;
        out|=((*in&0x7f)<<shamt);
        shamt+=7;
    }
    in++;
    return out;
}
template<typename output_t>
const uint8_t* parse_buffer(const uint8_t* in,output_t& out_all) {
    std::string sample_name;
    while(*in!=0) {
        sample_name.push_back(*in);
        in++;
    }
    auto out=out_all.set_name(std::move(sample_name));
    in++;
    while(*in) {
        int pos1=loadVariant(in);
        assert(*in);
        if (*(in+1)) {
            int pos2=loadVariant(in);
            out.add_Not_N(pos1,(*in)&0xf);
            out.add_Not_N(pos2,0xf&((*in)>>4));
        } else {
            out.add_Not_N(pos1,(*in)&0xf);
        }
        in++;
    }
    in++;

    while (*in) {
        int first=loadVariant(in);
        auto after_first=in;
        if (!(*in)) {
            out.add_N(first,first);
            break;
        }
        int second=loadVariant(in);
        if (first>second) {
            out.add_N(second,first);
        } else {
            out.add_N(first,first);
            in=after_first;
        }
    }
    in++;
    return in;
}
template<typename output_t>
struct printer {
    output_t& out;
    void operator()(const uint8_t* in) const {
        uint8_t buffer[MAX_SIZ];
        unsigned int item_len=*(int*) in;
        size_t out_len=MAX_SIZ;
        auto uncompress_out=uncompress(buffer, &out_len, in+4, item_len);
        if (uncompress_out!=Z_OK) exit(1);
        const uint8_t* start=buffer;
        auto end=buffer+out_len;
        while(start!=end) {
            start=parse_buffer(start,out);
        }
    }
};
class mapped_file {
    size_t size;
    const uint8_t *map_start;

  public:
    mapped_file(const char *path) {
        auto fh = open(path, O_RDONLY);
        struct stat stat_buf;
        fstat(fh, &stat_buf);
        size = stat_buf.st_size;
        map_start =
            (uint8_t *)mmap(nullptr, size, PROT_READ, MAP_SHARED, fh, 0);
        if (fh==-1) {
            map_start=nullptr;
            return;
        }
    }
    void get_mapped_range(const uint8_t *&start, const uint8_t *&end) const {
        start = map_start;
        end = map_start + size;
    }
    operator bool() const {
        return map_start!=0;
    }
    ~mapped_file() {
        if(map_start) munmap((void *)map_start, size);
    }
};
template <typename output_t>
static bool load_mutations(const char *path, int nthread, output_t &out) {
    mapped_file f(path);
    if (!f) {
        return false;
    }
    const uint8_t *last_out;
    const uint8_t *end;
    f.get_mapped_range(last_out, end);
    tbb::parallel_pipeline(
        nthread, tbb::make_filter<void, const uint8_t *>(
            tbb::filter::serial_in_order, partitioner{last_out, end}) &
        tbb::make_filter<const uint8_t *, void>(
            tbb::filter::parallel, printer<output_t> {out}));
    return true;
}
static void parse_rename_file(const std::string&  in_file_name, std::unordered_map<std::string,std::string>& mapping) {
    FILE* fd=fopen(in_file_name.c_str(),"r");
    char sample_name[BUFSIZ];
    char rename[BUFSIZ];
    int ret_val;
    while((ret_val=fscanf(fd,"%s	%s\n",sample_name,rename))!=EOF) {
        if(ret_val==2) {
            mapping[sample_name]=rename;
        }
    }
    fclose(fd);
}
