
cmake_minimum_required (VERSION 3.8) 

project(strain_phylogenetics)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -g -Wno-unused-function")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -fsanitize=address  -fsanitize=undefined -fsanitize=leak")

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

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem REQUIRED)

#find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

AUX_SOURCE_DIRECTORY(src/tree_rearrangements tree_rearrange_src) 
file(GLOB MATUTIL_SRCS "src/matUtils/*.cpp" "src/matUtils/*.hpp")

if(DEFINED Protobuf_PATH)
    add_executable(usher
        src/Mutation_Collection.cpp
        ${tree_rearrange_src}
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher.cpp
        )

    add_executable(matUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${MATUTIL_SRCS}
        )
    
    #add_executable(serializer_test
    #    src/serializer_check.cpp
    #    src/tree_rearrangements/Profitable_Move_Cacher.cpp
    #)
    
    add_executable(possible_muts_test
        src/tree_rearrangements/get_possible_muts_at_tip.cpp
    )
    
    protobuf_generate(
        LANGUAGE cpp
        TARGET usher 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET matUtils
        PROTOS parsimony.proto)
else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/Mutation_Collection.cpp
        ${tree_rearrange_src}
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(matUtils
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/Mutation_Collection.cpp
        ${MATUTIL_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )
        
    #add_executable(serializer_test
    #    src/serializer_check.cpp
    #    src/tree_rearrangements/Profitable_Move_Cacher.cpp
    #)
    
    add_executable(possible_muts_test
    src/tree_rearrangements/get_possible_muts_at_tip.cpp
    src/mutation_annotated_tree.cpp
    src/Mutation_Collection.cpp
    src/tree_rearrangements/Fitch_Sankoff.cpp
    src/tree_rearrangements/check_samples.cpp
    ${PROTO_SRCS}
    ${PROTO_HDRS}
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

TARGET_COMPILE_OPTIONS(possible_muts_test PRIVATE -DTEST -DNDEBUG)

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(possible_muts_test PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
