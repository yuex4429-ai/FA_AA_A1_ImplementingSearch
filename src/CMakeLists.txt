# src/CMakeLists.txt

add_library("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)

target_include_directories("${PROJECT_NAME}_interface" INTERFACE
    "${CMAKE_CURRENT_LIST_DIR}/../include"
    "${CEREAL_INCLUDE_DIR}"
)

target_compile_options("${PROJECT_NAME}_interface" INTERFACE
    "$<$<COMPILE_LANGUAGE:CXX>:-pedantic;-Wall;-Wextra>"
)

# -----------------------
# executables
# -----------------------
add_executable(naive_search naive_search.cpp)
target_link_libraries(naive_search PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_construct fmindex_construct.cpp)
target_link_libraries(fmindex_construct PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_search fmindex_search.cpp)
target_link_libraries(fmindex_search PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_pigeon_search fmindex_pigeon_search.cpp)
target_link_libraries(fmindex_pigeon_search PRIVATE "${PROJECT_NAME}_interface")

# -----------------------
# suffix array: divsufsort
# -----------------------
# In your setup, divsufsort headers are generated/placed in the BUILD tree:
#   ${CMAKE_BINARY_DIR}/lib/libdivsufsort/include/divsufsort.h
# Some environments also place headers under:
#   ${CMAKE_BINARY_DIR}/_deps/libdivsufsort-src/include
# So we add both as fallbacks.
set(DIVSUFSORT_BUILD_INCLUDE "${CMAKE_BINARY_DIR}/lib/libdivsufsort/include")
set(DIVSUFSORT_DEPS_INCLUDE  "${CMAKE_BINARY_DIR}/_deps/libdivsufsort-src/include")

add_executable(suffixarray_construct suffixarray_construct.cpp)
target_link_libraries(suffixarray_construct PRIVATE "${PROJECT_NAME}_interface" divsufsort)
target_include_directories(suffixarray_construct PRIVATE
    "${DIVSUFSORT_BUILD_INCLUDE}"
    "${DIVSUFSORT_DEPS_INCLUDE}"
)

add_executable(suffixarray_search suffixarray_search.cpp)
target_link_libraries(suffixarray_search PRIVATE "${PROJECT_NAME}_interface" divsufsort)
target_include_directories(suffixarray_search PRIVATE
    "${DIVSUFSORT_BUILD_INCLUDE}"
    "${DIVSUFSORT_DEPS_INCLUDE}"
)

