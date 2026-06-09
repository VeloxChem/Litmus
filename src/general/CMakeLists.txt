# General-purpose utilities (file streams, string formatting).
file(GLOB GENERAL_SOURCES CONFIGURE_DEPENDS *.cpp)

add_library(ltm_general OBJECT ${GENERAL_SOURCES})
target_link_libraries(ltm_general PUBLIC litmus_headers)
