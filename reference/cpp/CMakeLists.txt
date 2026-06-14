cmake_minimum_required(VERSION 3.16)
project(mapa_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(mapa mapa.cpp)
target_include_directories(mapa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()

add_executable(test_mapa test_mapa.cpp)
target_link_libraries(test_mapa PRIVATE mapa)
target_compile_definitions(test_mapa PRIVATE FIXTURES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../fixtures")

add_test(NAME test_mapa COMMAND test_mapa)
