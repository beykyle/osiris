# cliutils.gitlab.io/modern-cmake/chapters/basics.html
# https://docs.travis-ci.com/user/languages/cpp/#cmake
cmake_minimum_required(VERSION 3.18)

# TODO this is a bit hacky
if(DEFINED PROJECT_NAME)
 set(SUBPROJECT ON)
else()
  set(SUBPROJECT OFF)
  project("osiris" VERSION 0.1 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
# valloric.github.io/YouCompleteMe/#c-family-semantic-completion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# stackoverflow.com/a/6595001/5101335
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif ()

find_package(Git REQUIRED)
find_package(Threads)

# stackoverflow.com/a/24414345
if (MSVC)
  add_compile_options(/permissive- /std:c++17 /W4 /WX)
  add_compile_definitions(_USE_MATH_DEFINES)
else()
  #add_compile_options(-Wall -Wextra -Wpedantic -Werror)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage")
endif()

### External projects
include(FetchContent)

FetchContent_Declare(
  json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
Message("fetching nlohmann/json")
FetchContent_MakeAvailable(json)

FetchContent_Declare(
  xtl
  GIT_REPOSITORY https://github.com/xtensor-stack/xtl
  GIT_TAG 0.7.5
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
Message("fetching xtensor-stack/xtl")
FetchContent_MakeAvailable(xtl)

FetchContent_Declare(
  xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
  GIT_TAG 0.24.7
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
Message("fetching xtensor-stack/xtensor")
FetchContent_MakeAvailable(xtensor)

# directories
add_subdirectory(src)
add_subdirectory(tests/unit_tests)
if (NOT SUBPROJECT)
  add_subdirectory(exec)
endif()

option(BUILD_PY_MODULE "BUILD_PY_MODULE" OFF)
if(BUILD_PY_MODULE)

  find_package(Python REQUIRED  COMPONENTS Interpreter Development NumPy)
  find_package(pybind11 REQUIRED CONFIG)
  find_package(xtensor-python REQUIRED)

  if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
  endif()

  add_subdirectory(pybind)
endif()
