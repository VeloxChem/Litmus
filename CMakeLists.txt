# LITMUS: An Automated Molecular Integrals Generator
# Copyright 2022 Z. Rinkevicius, KTH, Sweden.

cmake_minimum_required(VERSION 3.16)

project(Litmus
        VERSION 1.0.0
        DESCRIPTION "Litmus: an automated molecular integrals generator"
        LANGUAGES CXX)

# Default to a Release build when the user has not chosen a build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

# Language requirements (Litmus uses C++17 features such as std::optional).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Emit compile_commands.json for editors and language servers.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Place the litmus.x executable in the top of the build tree.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)

# Unit tests (GoogleTest, fetched at configure time). Disable with
# -DLITMUS_BUILD_TESTS=OFF (e.g. for offline or packaging builds).
option(LITMUS_BUILD_TESTS "Build the Litmus unit tests" ON)
if(LITMUS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
