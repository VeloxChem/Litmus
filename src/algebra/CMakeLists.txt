# Symbolic algebra core (integrals, operators, tensors, recursion terms).
file(GLOB ALGEBRA_SOURCES CONFIGURE_DEPENDS *.cpp)

add_library(ltm_algebra OBJECT ${ALGEBRA_SOURCES})
target_link_libraries(ltm_algebra PUBLIC litmus_headers)
