
cmake_minimum_required (VERSION 3.8) 

project(strain_phylogenetics)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -g -Wno-unused-function")
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

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(matOptimize PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matOptimize PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
