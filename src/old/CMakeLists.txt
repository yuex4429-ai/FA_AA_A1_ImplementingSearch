# src/CMakeLists.txt

add_library("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)

target_include_directories("${PROJECT_NAME}_interface" INTERFACE
    "${CMAKE_CURRENT_LIST_DIR}/../include"
)

target_compile_options("${PROJECT_NAME}_interface" INTERFACE
    "$<$<COMPILE_LANGUAGE:CXX>:-pedantic;-Wall;-Wextra>"
)

add_executable(naive_search naive_search.cpp)
target_link_libraries(naive_search PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_construct fmindex_construct.cpp)
target_link_libraries(fmindex_construct PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_search fmindex_search.cpp)
target_link_libraries(fmindex_search PRIVATE "${PROJECT_NAME}_interface")

add_executable(fmindex_pigeon_search fmindex_pigeon_search.cpp)
target_link_libraries(fmindex_pigeon_search PRIVATE "${PROJECT_NAME}_interface")

add_executable(suffixarray_search suffixarray_search.cpp)
target_link_libraries(suffixarray_search PRIVATE "${PROJECT_NAME}_interface" divsufsort)
target_include_directories(suffixarray_search PRIVATE /home/mi/xiaoy95/opt/include)
