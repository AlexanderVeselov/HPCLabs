cmake_minimum_required(VERSION 2.8)

project(HPCLabs)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(lab6)
add_subdirectory(lab7_1)
add_subdirectory(lab7_2)
add_subdirectory(lab9)
add_subdirectory(lab10)
