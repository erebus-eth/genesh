cmake_minimum_required(VERSION 3.25)
project(Gens)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

enable_language(CUDA)

add_executable(Gens
        defs.h
        gen_compare.cpp
        gen_compare_gpu.cu
        gen_compare_gpu.h)
