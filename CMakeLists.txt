
cmake_minimum_required (VERSION 3.8) 
cmake_policy(SET CMP0069 NEW)
project(usher)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()
set(usher_VERSION "?.?.?")
find_library(ISAL_LIB isal)
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if ("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(usher_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
message(STATUS "${usher_VERSION}")

configure_file(src/version.hpp.in version.hpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations") 
#-DDEBUG_PARSIMONY_SCORE_CHANGE_CORRECT -fsanitize=address -DCHECK_PAR_MAIN  -DEASY_DEBUG -DSTOP_ON_ERROR  -DDETAILED_MERGER_CHECK -O0
set(CMAKE_CXX_FLAGS_DEBUG "  -fno-eliminate-unused-debug-symbols -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 ")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

if(DEFINED Protobuf_PATH) 
    find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
    find_package(Protobuf REQUIRED) 
endif()

# Print version messages
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${Protobuf_INCLUDE_DIRS})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

#find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
AUX_SOURCE_DIRECTORY(src/matOptimize/Profitable_Moves_Enumerators New_Profitable_Moves_Enumerators) 
AUX_SOURCE_DIRECTORY(src/matOptimize/apply_move patch_tree) 
AUX_SOURCE_DIRECTORY(src/check_samples_place check_samples_place) 
file(GLOB MATUTIL_SRCS "src/matUtils/*.cpp" "src/matUtils/*.hpp")
file(GLOB RIPPLES_SRCS "src/ripples/*.cpp" "src/ripples/*.hpp")
file(GLOB USHER_SAMPLED_SRCS "src/usher-sampled/static_tree_mapper/*" "src/usher-sampled/*.cpp" "src/usher-sampled/*.hpp")
file(GLOB RIPPLES_FAST_SRCS "src/ripples/ripples_fast/*.cpp" "src/ripples/ripples_fast/*.hpp")
file(GLOB RIPPLES_UTILS_SRCS "src/ripples/util/*.cpp" "src/ripples/util/*.hpp")

set_source_files_properties(src/mutation_annotated_tree.cpp PROPERTIES COMPILE_FLAGS -O3)
set_source_files_properties(src/usher-sampled/wait_debug.cpp PROPERTIES COMPILE_FLAGS -O0)
#set_source_files_properties(src/matOptimize/Profitable_Moves_Enumerators/range_tree.cpp PROPERTIES COMPILE_FLAGS -O0)

#set_source_files_properties(src/matOptimize/detailed_mutations_store.cpp PROPERTIES COMPILE_FLAGS -O0)
#set_source_files_properties(src/matOptimize/detailed_mutations_load.cpp PROPERTIES COMPILE_FLAGS -O0)
#set_source_files_properties(src/matOptimize/serializer_test.cpp PROPERTIES COMPILE_FLAGS -O0)
find_package(ZLIB)
if(DEFINED Protobuf_PATH)
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher_common.cpp
        src/usher.cpp
        )

    add_executable(matUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATUTIL_SRCS}
        )

    add_executable(ripples
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_SRCS}
        )
    add_executable(ripples-fast
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_FAST_SRCS}
        )
        add_executable(usher-sampled
        ${USHER_SAMPLED_SRCS}
        src/matOptimize/main_helper.cpp
        src/matOptimize/detailed_mutations_store.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/check_samples.cpp
        src/matOptimize/optimize_tree.cpp
        src/matOptimize/condense.cpp
        src/matOptimize/Fitch_Sankoff.cpp
        src/matOptimize/main_load_tree.cpp
        src/matOptimize/priority_conflict_resolver.cpp
        src/usher-sampled/driver/main.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
    )

    add_executable(ripplesUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/matUtils/describe.cpp
				${RIPPLES_UTILS_SRCS}
        )

    add_executable(ripplesInit
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/ripples/init/init_pipeline.cpp
        src/ripples/init/main.cpp
        )

    add_executable(matOptimize
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/detailed_mutations_store.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/main.cpp
        src/matOptimize/optimize_tree.cpp
        src/matOptimize/import_vcf_fast.cpp
        src/matOptimize/condense.cpp
        src/matOptimize/VCF_load_tree.cpp
        src/matOptimize/main_load_tree.cpp
        src/matOptimize/main_helper.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/Fitch_Sankoff.cpp
        src/matOptimize/check_samples.cpp
        #src/matOptimize/priority_conflict_resolver_cross_only.cpp
        src/matOptimize/priority_conflict_resolver.cpp
        src/matOptimize/transpose_vcf/transposed_vcf_patch.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
        )
    add_executable(compareVCF
        src/mutation_annotated_tree.cpp
        src/compareVCF.cpp
        src/usher_mapper.cpp
        )

    add_executable(output_final_protobuf
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/main_helper.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/output_final_protobuf.cpp
    )
    add_executable(transpose_vcf
        src/matOptimize/transpose_vcf/transpose_vcf_encode.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
    add_executable(transposed_vcf_to_vcf
        src/matOptimize/transpose_vcf/transposed_vcf_to_vcf.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
    add_executable(transposed_vcf_to_fa
    src/matOptimize/mutation_annotated_tree_nuc_util.cpp
    src/matOptimize/transpose_vcf/transposed_vcf_to_fa.cpp
    )
    add_executable(transposed_vcf_print_name
    src/matOptimize/transpose_vcf/transposed_vcf_print_name.cpp
    src/matOptimize/mutation_annotated_tree_nuc_util.cpp
    )
    if(USHER_SERVER)
        add_executable(usher_server
            src/mutation_annotated_tree.cpp
            src/usher_mapper.cpp
            src/usher_common.cpp
            src/usher_server.cpp
            )
    endif()

    protobuf_generate(
        LANGUAGE cpp
        TARGET ripples 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET ripples-fast 
        PROTOS parsimony.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET compareVCF 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET usher 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET matUtils
        PROTOS parsimony.proto)
    
    protobuf_generate(
        LANGUAGE cpp
        TARGET usher-sampled
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET ripplesUtils
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET ripplesInit
        PROTOS parsimony.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET matOptimize
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET matOptimize
        PROTOS mutation_detailed.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET usher-sampled
        PROTOS mutation_detailed.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET output_final_protobuf
        PROTOS parsimony.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET output_final_protobuf
        PROTOS mutation_detailed.proto)
    protobuf_generate(
        LANGUAGE cpp
        TARGET matUtils
        PROTOS taxodium.proto)
    
    if(USHER_SERVER)
        protobuf_generate(
            LANGUAGE cpp
            TARGET usher_server
            PROTOS parsimony.proto)
    endif()

else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    protobuf_generate_cpp(
        DETAILED_MUTATIONS_PROTO_SRCS
        DETAILED_MUTATIONS_PROTO_HDRS
        mutation_detailed.proto)

    protobuf_generate_cpp(
        TAXO_SRCS TAXO_HDRS
        taxodium.proto)
    
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher_common.cpp
        src/usher.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(matUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATUTIL_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${TAXO_SRCS}
        ${TAXO_HDRS}
        )

    add_executable(ripples
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(ripples-fast
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_FAST_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )
    add_executable(usher-sampled
        ${USHER_SAMPLED_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
        src/matOptimize/main_helper.cpp
        src/matOptimize/detailed_mutations_store.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/check_samples.cpp
        src/matOptimize/optimize_tree.cpp
        src/matOptimize/condense.cpp
        src/matOptimize/Fitch_Sankoff.cpp
        src/matOptimize/main_load_tree.cpp
        src/matOptimize/priority_conflict_resolver.cpp
        src/usher-sampled/driver/main.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
    )

    add_executable(ripplesUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/matUtils/describe.cpp
				${RIPPLES_UTILS_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(ripplesInit
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/ripples/init/init_pipeline.cpp
        src/ripples/init/main.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(compareVCF
        src/mutation_annotated_tree.cpp
        src/compareVCF.cpp
        src/usher_mapper.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(matOptimize
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/detailed_mutations_store.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/main.cpp
        src/matOptimize/optimize_tree.cpp
        src/matOptimize/import_vcf_fast.cpp
        src/matOptimize/condense.cpp
        src/matOptimize/VCF_load_tree.cpp
        src/matOptimize/main_load_tree.cpp
        src/matOptimize/main_helper.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/Fitch_Sankoff.cpp
        src/matOptimize/check_samples.cpp
        #src/matOptimize/priority_conflict_resolver_cross_only.cpp
        src/matOptimize/priority_conflict_resolver.cpp
        src/matOptimize/transpose_vcf/transposed_vcf_patch.cpp
        src/usher-sampled/utils.cpp
        src/usher-sampled/Min_back_FS.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
        )
        
        add_executable(check_samples_place
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/Mutation_Collection.cpp
        ${check_samples_place}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
        )
    
    #[[add_executable(output_final_protobuf
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/detailed_mutations_load.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/main_helper.cpp
        src/matOptimize/Mutation_Collection.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/output_final_protobuf.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
    )]]
    add_executable(transpose_vcf
        src/matOptimize/transpose_vcf/transpose_vcf_encode.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
    add_executable(transposed_vcf_to_vcf
        src/matOptimize/transpose_vcf/transposed_vcf_to_vcf.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
        add_executable(transposed_vcf_to_fa
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        src/matOptimize/transpose_vcf/transposed_vcf_to_fa.cpp
        )
        add_executable(transposed_vcf_print_name
        src/matOptimize/transpose_vcf/transposed_vcf_print_name.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
    if(USHER_SERVER)
        add_executable(usher_server
            src/mutation_annotated_tree.cpp
            src/usher_mapper.cpp
            src/usher_common.cpp
            src/usher_server.cpp
            ${PROTO_SRCS}
            ${PROTO_HDRS}
            )
            add_executable(usher-sampled-server
            ${USHER_SAMPLED_SRCS}
            ${PROTO_SRCS}
            ${PROTO_HDRS}
            ${DETAILED_MUTATIONS_PROTO_SRCS}
            ${DETAILED_MUTATIONS_PROTO_HDRS}
            src/matOptimize/main_helper.cpp
            src/matOptimize/detailed_mutations_store.cpp
            src/matOptimize/detailed_mutations_load.cpp
            src/matOptimize/mutation_annotated_tree.cpp
            src/matOptimize/mutation_annotated_tree_node.cpp
            src/matOptimize/mutation_annotated_tree_load_store.cpp
            src/matOptimize/mutation_annotated_tree_nuc_util.cpp
            src/matOptimize/Mutation_Collection.cpp
            src/matOptimize/check_samples.cpp
            src/matOptimize/optimize_tree.cpp
            src/matOptimize/condense.cpp
            src/matOptimize/Fitch_Sankoff.cpp
            src/matOptimize/main_load_tree.cpp
            src/matOptimize/priority_conflict_resolver.cpp
            src/usher-sampled/driver/socket.cpp
            ${patch_tree}
            ${New_Profitable_Moves_Enumerators}
        )
    endif()
endif()

if(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1 -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1)
    endif(DEBUG)
else(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    endif(DEBUG)
endif(SAVE_PROFILE)
#set_property(TARGET matOptimize PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
TARGET_LINK_LIBRARIES(compareVCF PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)


TARGET_LINK_LIBRARIES(check_samples_place PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB} ) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(matOptimize PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB} ) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(usher-sampled PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB} ) # OpenMP::OpenMP_CXX)
#TARGET_LINK_LIBRARIES(output_final_protobuf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transpose_vcf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transposed_vcf_to_vcf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transposed_vcf_to_fa PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transposed_vcf_print_name PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(usher PUBLIC "${PROJECT_BINARY_DIR}")

TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(ripplesUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(ripplesInit PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(ripples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(ripples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(ripples-fast PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

if(USHER_SERVER)
    TARGET_COMPILE_OPTIONS(usher_server PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    TARGET_LINK_LIBRARIES(usher_server PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
    TARGET_LINK_LIBRARIES(usher-sampled-server PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB} ) # OpenMP::OpenMP_CXX)
    install(TARGETS usher matUtils matOptimize ripples usher_server DESTINATION bin)
else()
    install(TARGETS usher matUtils matOptimize ripples DESTINATION bin)
endif()
