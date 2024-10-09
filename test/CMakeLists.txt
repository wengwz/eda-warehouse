
add_executable(test_FM test_FM.cpp)
target_include_directories(test_FM PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_FM PUBLIC ${PROJECT_SOURCE_DIR}/src/partition)
target_link_libraries(test_FM partition)


add_executable(test_BK test_BucketArray.cpp)
target_include_directories(test_BK PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_BK PUBLIC ${PROJECT_SOURCE_DIR}/src/partition)
target_link_libraries(test_BK partition)