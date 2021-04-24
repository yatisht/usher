#include <fstream>
#include <numeric>
#include <algorithm>
#include <boost/program_options.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/tbb.h>
#include "../mutation_annotated_tree.hpp"

namespace po = boost::program_options;
namespace MAT = Mutation_Annotated_Tree;

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

int main(int argc, char** argv) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;
    std::string inp_vcf_filename;
    //std::string out_vcf_filename;
    
    po::options_description desc{"Options"};
    
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    desc.add_options()
        ("in-vcf,v", po::value<std::string>(&inp_vcf_filename)->required(), "Input VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
//        ("out-vcf,o", po::value<std::string>(&out_vcf_filename)->required(), "Output VCF file (in uncompressed or gzip-compressed .gz format) [REQUIRED]")
        ("threads,T", po::value<uint32_t>(&num_threads)->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
    
    po::options_description all_options;
    all_options.add(desc);
    
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    // Boost library used to stream the contents of the input VCF file in
    // uncompressed or compressed .gz format
    std::ifstream infile(inp_vcf_filename, std::ios_base::in | std::ios_base::binary);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", inp_vcf_filename.c_str());
        exit(1);
    }
    boost::iostreams::filtering_istream instream;
    try {
        if (inp_vcf_filename.find(".gz\0") != std::string::npos) {
            instream.push(boost::iostreams::gzip_decompressor());
        }
        instream.push(infile);
    }
    catch(const boost::iostreams::gzip_error& e) {
        std::cout << e.what() << '\n';
    }

    bool header_found = false;
    std::vector<size_t> indices;
    while (true) {
        //check if reached end-of-file
        int curr_char = instream.peek();
        if(curr_char == EOF)
            return false;

        std::string s;
        std::getline(instream, s);
        std::vector<std::string> words;
        MAT::string_split(s, words);
        
        // Header is found when "POS" is the second word in the line
        if (not header_found) {
            if ((words.size()>1) && (words[1] == "POS")) {
                std::vector<std::string> sample_names = {words.begin()+9, words.end()};
                indices = sort_indexes(sample_names);
                //for (auto i: indices) {
                //fprintf(stderr, "%zu\n", i);
                //}
                header_found = true;
            }
            else {
                fprintf(stdout, "%s\n", s.c_str());
            }
        }
        if (header_found) {
            for (size_t i=0; i<9; i++) {
                fprintf(stdout, "%s\t", words[i].c_str()); 
            }
            for (auto& i: indices) {
                fprintf(stdout, "%s", words[i+9].c_str()); 
                if (&i!=&indices.back()) {
                    fprintf(stdout, "\t");
                }
            }
            fprintf(stdout, "\n");
        }
    }

    return 0;
}

