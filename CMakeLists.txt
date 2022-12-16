cmake_minimum_required(VERSION 3.23)
project(cpp8v)

set(CMAKE_CXX_STANDARD 14)
add_executable(cpp8v main.cpp SparseMatrix.h SquareSparseMatrix.h)