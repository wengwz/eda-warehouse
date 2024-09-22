aux_source_directory(. SRC_FILES)
add_library(partition ${SRC_FILES})
target_include_directories(partition PUBLIC ${PROJECT_SOURCE_DIR}/include)
