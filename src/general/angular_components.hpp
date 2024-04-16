#ifndef AngularComponents_hpp
#define AngularComponents_hpp

#include <array>
#include <cstddef>

namespace ten {  // ten

/// Gets number of Cartesian components in canonical tensor.
/// - Parameter order: the order of canonical tensor.
inline auto
number_of_cartesian_components(const int order) -> int
{
    return (order + 1) * (order + 2) / 2;
}

/// Gets number of spherical components in canonical tensor.
/// - Parameter order: the order of canonical tensor.
inline auto
number_of_spherical_components(const int order) -> int
{
    return 2 * order + 1;
}

/// Gets compound number of Cartesian components of canonical tensors array.
/// - Parameter order: the array of orders of canonical tensors.
template <std::size_t N>
auto
number_of_cartesian_components(const std::array<int, N>& orders) -> int
{
    int ncomps = 1;

    for (std::size_t i = 0; i < N; i++)
    {
        ncomps *= ten::number_of_cartesian_components(orders[i]);
    }

    return ncomps;
}

/// Gets compound number of spherical components of canonical tensors array.
/// - Parameter order: the array of orders of canonical tensors.
template <std::size_t N>
auto
number_of_spherical_components(const std::array<int, N>& orders) -> int
{
    int ncomps = 1;

    for (std::size_t i = 0; i < N; i++)
    {
        ncomps *= ten::number_of_spherical_components(orders[i]);
    }

    return ncomps;
}

}  // namespace ten

#endif /* AngularComponents_hpp */
