cmake_minimum_required(VERSION 3.18)
project(ImplementingSearch VERSION 1.0.0 LANGUAGES CXX)

## BUILD TYPE
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

## OUTPUT DIRS
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

## C++ STANDARD (SeqAn3 requires C++20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

## Dependency: SeqAn3
find_package(seqan3 3.2 REQUIRED
    HINTS "${CMAKE_CURRENT_LIST_DIR}/lib/seqan3/build_system"
)

## cereal shipped with SeqAn3 submodules (header-only)
set(CEREAL_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/lib/seqan3/submodules/cereal/include")

## Dependency: libdivsufsort
option(BUILD_EXAMPLES "" OFF)
add_subdirectory(lib/libdivsufsort)

## Your targets
add_subdirectory(src)

