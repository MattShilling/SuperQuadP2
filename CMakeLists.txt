cmake_minimum_required(VERSION 2.8)

# set the project name
project(MonteCarlo CXX)

add_compile_options(-std=c++11)

find_package(OpenMP REQUIRED)

# add the executable
add_executable(MonteCarlo main.cc test_rig.cc)

target_link_libraries(MonteCarlo PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(MonteCarlo PRIVATE "${OpenMP_CXX_FLAGS}")