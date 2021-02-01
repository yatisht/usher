
cmake_minimum_required (VERSION 3.8) 

project(strain_phylogenetics)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -Wno-unused-function -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG " -O0 -fsanitize=address -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE " -DNDEBUG")
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

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem REQUIRED)

#find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

if(DEFINED Protobuf_PATH)
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/tree_rearrangement.cpp
        src/fitch-sankoff.cpp
        src/NNI.cpp
        src/usher.cpp
        )

    add_executable(matToVcf
        src/mutation_annotated_tree.cpp
        src/matToVcf.cpp
        )

    if(COMPILE_UTILS)
        add_executable(matUtils
            src/mutation_annotated_tree.cpp
            src/matUtils.cpp
            )
    endif(COMPILE_UTILS)

    protobuf_generate(
        LANGUAGE cpp
        TARGET usher 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET matToVcf
        PROTOS parsimony.proto)

    if(COMPILE_UTILS)
        protobuf_generate(
            LANGUAGE cpp
            TARGET matUtils
            PROTOS parsimony.proto)
    endif(COMPILE_UTILS)
else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/Fitch_Sankoff.cpp
        src/check_samples.cpp
        src/tree_rearrangement.cpp
        src/NNI.cpp
        src/usher.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(matToVcf
        src/mutation_annotated_tree.cpp
        src/matToVcf.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    if(COMPILE_UTILS)
        add_executable(matUtils
            src/mutation_annotated_tree.cpp
            src/matUtils.cpp
            ${PROTO_SRCS}
            ${PROTO_HDRS}
            )
    endif(COMPILE_UTILS)
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

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

TARGET_COMPILE_OPTIONS(matToVcf PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(matToVcf PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)

if(COMPILE_UTILS)
    TARGET_COMPILE_OPTIONS(matUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    TARGET_LINK_LIBRARIES(matUtils PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
endif(COMPILE_UTILS)
