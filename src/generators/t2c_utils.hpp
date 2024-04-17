// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef t2c_utils_hpp
#define t2c_utils_hpp

#include <map>
#include <string>
#include <array>
#include <set>
#include <vector>
#include <utility>

#include "t2c_defs.hpp"

namespace t2c { // t2c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base two center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I2CIntegral& integral);

/// Gets recursion namespace label of standart integral.
/// @param integral The base two center integral.
/// @return The recursion namespace label of standart integral.
std::string namespace_label(const I2CIntegral& integral);

/// Gets label of standart integrand.
/// @param integrand the integrand operator.
/// @return The label of standart integrand.
std::string integrand_label(const Operator& integrand);

/// Gets standart label of integral prefixes.
/// @param integral The base two center integral.
/// @return The standart capitalized label of integral.
std::pair<std::string, std::string> prefixes_label(const I2CIntegral& integral);

/// Gets all labels of integrand with specific prefix.
/// @param integral The base two center integral.
/// @param prefix The prefix to label of integrand.
/// @return The vector of integrand labels.
std::vector<std::string> integrand_labels(const I2CIntegral& integral,
                                          const std::string& prefix);

/// Generates compute function  name.
/// @param integral The base two center integral.
/// @param rec_form The recursion form for two center integrals (summation, convolution flags).
/// @return The compute function name.
std::string compute_func_name(const I2CIntegral&           integral,
                              const std::pair<bool, bool>& rec_form);

/// Generates primitive file name.
/// @param integral The base two center integral.
/// @return The primitive file name.
std::string prim_file_name(const I2CIntegral& integral);

/// Gets number of Cartesian components in canonical tensor.
/// @param order the order of canonical tensor.
/// @return The number of Cartesian components.
inline auto
number_of_cartesian_components(const int order) -> int
{
    return (order + 1) * (order + 2) / 2;
}

/// Gets number of spherical components in canonical tensor.
/// @param order the order of canonical tensor.
/// @return The number of spherical components.
inline auto
number_of_spherical_components(const int order) -> int
{
    return 2 * order + 1;
}

/// Gets compound number of Cartesian components of canonical tensors array.
/// @param orders the array of orders of canonical tensor.
/// @return The number of Cartesian components.
template <std::size_t N>
auto
number_of_cartesian_components(const std::array<int, N>& orders) -> int
{
    int ncomps = 1;

    for (std::size_t i = 0; i < N; i++)
    {
        ncomps *= t2c::number_of_cartesian_components(orders[i]);
    }

    return ncomps;
}

/// Gets compound number of spherical components of canonical tensors array.
/// @param orders the array of orders of canonical tensor.
/// @return The number of spherical components.
template <std::size_t N>
auto
number_of_spherical_components(const std::array<int, N>& orders) -> int
{
    int ncomps = 1;

    for (std::size_t i = 0; i < N; i++)
    {
        ncomps *= t2c::number_of_spherical_components(orders[i]);
    }

    return ncomps;
}

/// Generates integral buffer label.
/// @param integral The base two center integral.
/// @return The string with integral label.
std::string get_buffer_label(const I2CIntegral& integral,
                             const std::string& prefix);

/// Generates primitive compute function name.
/// @param integral The base two center integral.
/// @return The primitive compute function name.
std::string prim_compute_func_name(const I2CIntegral& integral);

/// Gets arguments list for primitive function call.
/// @param integral The base two center integral.
SI2CIntegrals get_integrals(const I2CIntegral& integral);

} // t2c namespace

#endif /* t2c_utils_hpp */
