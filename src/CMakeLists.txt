# LITMUS: An Automated Molecular Integrals Generator
# Copyright 2022 Z. Rinkevicius, KTH, Sweden.

# All Litmus headers are referenced by their bare file name and the source
# modules cross-include each other freely (e.g. recursions uses the *_utils
# headers from generators). Expose every module directory through a single
# INTERFACE target so each module and the executable see the full header set.
add_library(litmus_headers INTERFACE)
target_include_directories(litmus_headers INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/general
    ${CMAKE_CURRENT_SOURCE_DIR}/algebra
    ${CMAKE_CURRENT_SOURCE_DIR}/recursions
    ${CMAKE_CURRENT_SOURCE_DIR}/generators)

# Each module is an OBJECT library. Because recursions and generators include
# one another's headers, OBJECT libraries are used so all objects are linked
# directly into the executable, avoiding static-archive link-order issues.
add_subdirectory(general)
add_subdirectory(algebra)
add_subdirectory(recursions)
add_subdirectory(generators)

# Litmus driver executable.
add_executable(litmus.x litmus.cpp)
target_link_libraries(litmus.x PRIVATE
    litmus_headers
    ltm_general
    ltm_algebra
    ltm_recursions
    ltm_generators)
