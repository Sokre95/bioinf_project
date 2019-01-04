cmake_minimum_required(VERSION 3.13)
project(bioinf)

set(CMAKE_CXX_STANDARD 14)

add_executable(bioinf main.cpp FastaParser.h FastaParser.cpp Sequence.cpp Sequence.h)