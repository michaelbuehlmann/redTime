cmake_minimum_required(VERSION 3.11)
project(TimeRG LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

find_package(OpenMP)
find_package(GSL REQUIRED)

add_subdirectory(src)