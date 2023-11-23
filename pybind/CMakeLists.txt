pybind11_add_module(osiris main.cpp)
target_link_libraries(osiris PUBLIC pybind11::module xtensor-python Python::NumPy osiris_lib)
target_include_directories(osiris PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_definitions(osiris PRIVATE VERSION_INFO=${OSIRIS_VERSION_INFO})
