find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
Message("fetching pybind11")
FetchContent_MakeAvailable(pybind11)

FetchContent_Declare(
  xtensor-python
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-python.git
  GIT_TAG 0.26.1
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
Message("fetching xtensor-stack/xtensor-python")
FetchContent_MakeAvailable(xtensor-python)

python_add_library(osiris_core MODULE osiris_core_module.cpp WITH_SOABI)
target_link_libraries(osiris_core PUBLIC pybind11::headers xtensor-python Python::NumPy osiris_lib)
target_include_directories(osiris_core PUBLIC ${CMAKE_SOURCE_DIR}/src/osiris_lib)
target_compile_definitions(osiris_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
set_target_properties(osiris_core PROPERTIES CXX_VISIBILITY_PRESET hidden)
install(TARGETS osiris_core DESTINATION osiris)
