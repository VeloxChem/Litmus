# Code generators that emit C++ integral routines from the derived recursions.
file(GLOB GENERATORS_SOURCES CONFIGURE_DEPENDS *.cpp)

add_library(ltm_generators OBJECT ${GENERATORS_SOURCES})
target_link_libraries(ltm_generators PUBLIC litmus_headers)
