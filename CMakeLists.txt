cmake_minimum_required(VERSION 3.13.0)
project(pycold LANGUAGES C CXX)


# Backend C library
add_subdirectory("src")

# Python bindings
add_subdirectory("tool/python/pycold")
