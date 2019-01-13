cmake_minimum_required(VERSION 3.10.3)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_INSTALL_PREFIX ./cmake-build-debug)

project(bioinf)

add_executable(bioinf main.cpp FastaParser.h FastaParser.cpp Sequence.cpp Sequence.h Viterbi.cpp Viterbi.h IViterbi.h ViterbiLogOdds.cpp ViterbiLogOdds.h MleEstimator.cpp MleEstimator.h)

install(TARGETS bioinf DESTINATION ./)