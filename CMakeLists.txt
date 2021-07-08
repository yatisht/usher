
cmake_minimum_required (VERSION 3.8) 
cmake_policy(SET CMP0069 NEW)
project(usher)

execute_process(COMMAND git describe --tags --abbrev=0 OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
set(usher_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

configure_file(src/version.hpp.in version.hpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.1 -std=c++11 -Wall -g -Wno-unused-function -Wno-deprecated-declarations")
#-DDEBUG_PARSIMONY_SCORE_CHANGE_CORRECT
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -DCHECK_PLACEMENT")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
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

include_directories(${Protobuf_INCLUDE_DIRS})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

#find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
AUX_SOURCE_DIRECTORY(src/new_tree_rearrangements/Profitable_Moves_Enumerators New_Profitable_Moves_Enumerators) 
AUX_SOURCE_DIRECTORY(src/new_tree_rearrangements/apply_move patch_tree) 
file(GLOB MATUTIL_SRCS "src/matUtils/*.cpp" "src/matUtils/*.hpp")
file(GLOB MATOPTIMIZE_SRCS "src/matOptimize/*.cpp" "src/matOptimize/*.hpp")
find_package(ZLIB)
if(DEFINED Protobuf_PATH)
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher.cpp
        )

    add_executable(matUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATUTIL_SRCS}
        )

    add_executable(matOptimize
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATOPTIMIZE_SRCS}
        )

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
        TARGET matOptimize
        PROTOS parsimony.proto)
    
else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    protobuf_generate_cpp(
        DETAILED_MUTATIONS_PROTO_SRCS
        DETAILED_MUTATIONS_PROTO_HDRS
        mutation_detailed.proto)

    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
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
        )

    add_executable(matOptimize
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATOPTIMIZE_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(compareVCF
        src/mutation_annotated_tree.cpp
        src/compareVCF.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(rename_fasta
        src/mutation_annotated_tree.cpp
        src/new_tree_rearrangements/rename_fasta.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )
    add_executable(tree_rearrange_new
        src/new_tree_rearrangements/mutation_annotated_tree.cpp
        src/new_tree_rearrangements/mutation_annotated_tree_node.cpp
        src/new_tree_rearrangements/mutation_annotated_tree_load_store.cpp
        src/new_tree_rearrangements/detailed_mutations_load_store.cpp
        src/new_tree_rearrangements/mutation_annotated_tree_nuc_util.cpp
        src/new_tree_rearrangements/main.cpp
        src/new_tree_rearrangements/optimize_tree.cpp
        src/new_tree_rearrangements/import_vcf_fast.cpp
        src/new_tree_rearrangements/condense.cpp
        src/new_tree_rearrangements/VCF_load_tree.cpp
        src/new_tree_rearrangements/main_load_tree.cpp
        src/new_tree_rearrangements/main_helper.cpp
        src/new_tree_rearrangements/Mutation_Collection.cpp
        src/new_tree_rearrangements/Fitch_Sankoff.cpp
        src/new_tree_rearrangements/check_samples.cpp
        #src/new_tree_rearrangements/priority_conflict_resolver_cross_only.cpp
        src/new_tree_rearrangements/priority_conflict_resolver.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
        )
        
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

#set_property(TARGET tree_rearrange_new PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
TARGET_LINK_LIBRARIES(compareVCF PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(tree_rearrange_new PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(rename_fasta PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(usher PUBLIC "${PROJECT_BINARY_DIR}")

TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(matOptimize PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matOptimize PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

install(TARGETS usher matUtils matOptimize DESTINATION bin)
