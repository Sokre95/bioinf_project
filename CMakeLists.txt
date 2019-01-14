cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_INSTALL_PREFIX ./build)

project(bioinf_project)

include_directories(include)
file(GLOB SOURCES "src/*")

add_executable(bioinf ${SOURCES})

install(TARGETS bioinf DESTINATION /)