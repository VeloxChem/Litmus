# Recursion drivers that derive Obara-Saika-type relations per integral type.
file(GLOB RECURSIONS_SOURCES CONFIGURE_DEPENDS *.cpp)

add_library(ltm_recursions OBJECT ${RECURSIONS_SOURCES})
target_link_libraries(ltm_recursions PUBLIC litmus_headers)
