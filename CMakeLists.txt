
cmake_minimum_required (VERSION 3.8) 
cmake_policy(SET CMP0069 NEW)
project(usher)

set(usher_VERSION "?.?.?")
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if ("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(usher_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
message(STATUS "${usher_VERSION}")

configure_file(src/version.hpp.in version.hpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.1 -std=c++11 -Wall -g -Wno-unused-function -Wno-deprecated-declarations")
#-DDEBUG_PARSIMONY_SCORE_CHANGE_CORRECT 
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -fsanitize=address")
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
AUX_SOURCE_DIRECTORY(src/matOptimize/Profitable_Moves_Enumerators New_Profitable_Moves_Enumerators) 
AUX_SOURCE_DIRECTORY(src/matOptimize/apply_move patch_tree) 
file(GLOB MATUTIL_SRCS "src/matUtils/*.cpp" "src/matUtils/*.hpp")
file(GLOB RIPPLES_SRCS "src/ripples/*.cpp" "src/ripples/*.hpp")

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

    add_executable(ripples
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_SRCS}
        )

    add_executable(matOptimize
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/detailed_mutations_load_store.cpp
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
        src/matOptimize/transpose_vcf/transposed_vcf_patch.cpp
        #src/matOptimize/priority_conflict_resolver_cross_only.cpp
        src/matOptimize/priority_conflict_resolver.cpp
        ${patch_tree}
        ${New_Profitable_Moves_Enumerators}
        )
    add_executable(compareVCF
        src/mutation_annotated_tree.cpp
        src/compareVCF.cpp)


        add_executable(transpose_vcf
        src/matOptimize/transpose_vcf/transpose_vcf_encode.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
        add_executable(transposed_vcf_to_vcf
        src/matOptimize/transpose_vcf/transposed_vcf_to_vcf.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )


        protobuf_generate(
        LANGUAGE cpp
        TARGET ripples 
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
        TARGET matOptimize
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET matOptimize
        PROTOS mutation_detailed.proto)
    
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
    add_executable(ripples
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${RIPPLES_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )


    add_executable(compareVCF
        src/mutation_annotated_tree.cpp
        src/compareVCF.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(matOptimize
        src/matOptimize/mutation_annotated_tree.cpp
        src/matOptimize/mutation_annotated_tree_node.cpp
        src/matOptimize/mutation_annotated_tree_load_store.cpp
        src/matOptimize/detailed_mutations_load_store.cpp
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
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        ${DETAILED_MUTATIONS_PROTO_SRCS}
        ${DETAILED_MUTATIONS_PROTO_HDRS}
        )
        add_executable(transpose_vcf
        src/matOptimize/transpose_vcf/transpose_vcf_encode.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
        )
        add_executable(transposed_vcf_to_vcf
        src/matOptimize/transpose_vcf/transposed_vcf_to_vcf.cpp
        src/matOptimize/mutation_annotated_tree_nuc_util.cpp
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

#set_property(TARGET matOptimize PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
TARGET_LINK_LIBRARIES(compareVCF PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(matOptimize PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transpose_vcf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(transposed_vcf_to_vcf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(usher PUBLIC "${PROJECT_BINARY_DIR}")

TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(ripples PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(ripples PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

install(TARGETS usher matUtils matOptimize ripples DESTINATION bin)
