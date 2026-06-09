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

# Tests for the symbolic algebra core (concrete value types).
add_executable(algebra_tests
    algebra/test_fraction.cpp
    algebra/test_tensor_component.cpp
    algebra/test_tensor.cpp
    algebra/test_factor.cpp
    algebra/test_operator.cpp
    algebra/test_integral_component.cpp)

target_link_libraries(algebra_tests PRIVATE
    GTest::gtest_main
    litmus_headers
    ltm_algebra)

gtest_discover_tests(algebra_tests)
