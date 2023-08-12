cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

project(psol)

add_library(psol INTERFACE)
target_include_directories(psol INTERFACE ${PROJECT_SOURCE_DIR}/psol/)

if (PROJECT_IS_TOP_LEVEL)
  find_package(doctest CONFIG REQUIRED)
  file(GLOB test_src "${CMAKE_SOURCE_DIR}/test/**.cpp")
  add_executable(psol_test ${test_src})
  target_link_libraries(psol_test PUBLIC doctest::doctest psol)
endif()



