file(GLOB_RECURSE SRC_CPP_FILES "*.cpp")

add_library(osiris_lib ${SRC_CPP_FILES})
set_property(TARGET osiris_lib PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(osiris_lib PUBLIC ${CMAKE_SOURCE_DIR}/osiris_lib)
target_include_directories(osiris_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(osiris_lib PUBLIC nlohmann_json::nlohmann_json ${CMAKE_THREAD_LIBS_INIT} xtensor)
set_target_properties(osiris_lib PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
