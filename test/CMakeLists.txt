
add_executable(test_FM test_FM.cpp)
target_include_directories(test_FM PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_FM PUBLIC ${PROJECT_SOURCE_DIR}/src/partition)
target_link_libraries(test_FM partition)


add_executable(test_BK test_BucketArray.cpp)
target_include_directories(test_BK PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_BK PUBLIC ${PROJECT_SOURCE_DIR}/src/partition)
target_link_libraries(test_BK partition)

find_package(nlohmann_json 3.2.0 REQUIRED)
add_executable(test_PositionAwarePartition test_PositionAwarePartition.cpp)
target_link_libraries(test_PositionAwarePartition partition nlohmann_json::nlohmann_json)
target_include_directories(test_PositionAwarePartition PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_PositionAwarePartition PUBLIC ${PROJECT_SOURCE_DIR}/src/partition)