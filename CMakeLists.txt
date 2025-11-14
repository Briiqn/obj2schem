cmake_minimum_required(VERSION 3.18)
project(ObjToSchem LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

find_package(ZLIB REQUIRED)

add_executable(obj2schem main.cu)

target_include_directories(obj2schem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(obj2schem PRIVATE ZLIB::ZLIB)

set_target_properties(obj2schem PROPERTIES CUDA_ARCHITECTURES "native")
