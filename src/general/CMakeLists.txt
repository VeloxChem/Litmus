# General-purpose utilities (file streams, string formatting, provenance).
file(GLOB GENERAL_SOURCES CONFIGURE_DEPENDS *.cpp)

# Regenerate litmus_version.hpp (the git revision stamped into every generated
# file) on each build; the script leaves the header untouched if unchanged.
set(LITMUS_VERSION_HEADER ${CMAKE_CURRENT_BINARY_DIR}/litmus_version.hpp)

add_custom_target(litmus_version
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT=${LITMUS_VERSION_HEADER}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/litmus_version.cmake
    BYPRODUCTS ${LITMUS_VERSION_HEADER}
    COMMENT "Recording Litmus git revision"
    VERBATIM)

add_library(ltm_general OBJECT ${GENERAL_SOURCES})
target_link_libraries(ltm_general PUBLIC litmus_headers)
target_include_directories(ltm_general PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(ltm_general litmus_version)
