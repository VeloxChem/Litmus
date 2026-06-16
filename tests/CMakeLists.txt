# Litmus unit tests.
#
# GoogleTest is fetched at configure time. This requires network access on the
# first configure; it is then cached under the build tree.
include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)

# Keep GoogleTest from clobbering the parent project's settings (Windows CRT).
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Tests for the symbolic algebra core (concrete value types and templates).
add_executable(algebra_tests
    algebra/test_axes.cpp
    algebra/test_components.cpp
    algebra/test_fraction.cpp
    algebra/test_spherical_harmonics.cpp
    algebra/test_tensor_component.cpp
    algebra/test_tensor.cpp
    algebra/test_factor.cpp
    algebra/test_operator.cpp
    algebra/test_operator_component.cpp
    algebra/test_one_center.cpp
    algebra/test_one_center_component.cpp
    algebra/test_two_center_pair.cpp
    algebra/test_two_center_pair_component.cpp
    algebra/test_integral.cpp
    algebra/test_integral_component.cpp
    algebra/test_recursion_term.cpp
    algebra/test_recursion_expansion.cpp
    algebra/test_recursion_group.cpp)

target_link_libraries(algebra_tests PRIVATE
    GTest::gtest_main
    litmus_headers
    ltm_algebra)

gtest_discover_tests(algebra_tests)

# Tests for the general-purpose utilities (file streams, string formatting).
add_executable(general_tests
    general/test_string_formater.cpp
    general/test_file_stream.cpp
    general/test_config.cpp
    general/test_run_configuration.cpp)

target_link_libraries(general_tests PRIVATE
    GTest::gtest_main
    litmus_headers
    ltm_general)

gtest_discover_tests(general_tests)

# Tests for the code generators (file-emitting kernels). Globbed so new
# generator tests are picked up automatically on reconfigure.
file(GLOB GENERATOR_TEST_SOURCES CONFIGURE_DEPENDS generators/*.cpp)

add_executable(generator_tests ${GENERATOR_TEST_SOURCES})

target_link_libraries(generator_tests PRIVATE
    GTest::gtest_main
    litmus_headers
    ltm_generators
    ltm_recursions
    ltm_algebra
    ltm_general)

gtest_discover_tests(generator_tests)

# Tests for the Obara-Saika recursion drivers. One test file per driver; globbed
# so new driver tests are picked up automatically on reconfigure.
file(GLOB RECURSION_TEST_SOURCES CONFIGURE_DEPENDS recursions/*.cpp)

add_executable(recursion_tests ${RECURSION_TEST_SOURCES})

target_link_libraries(recursion_tests PRIVATE
    GTest::gtest_main
    litmus_headers
    ltm_recursions
    ltm_algebra
    ltm_general)

gtest_discover_tests(recursion_tests)
