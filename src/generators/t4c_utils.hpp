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

#ifndef t4c_utils_hpp
#define t4c_utils_hpp

#include <string>
#include <utility>

#include "t4c_defs.hpp"

namespace t4c { // t4c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base four center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I4CIntegral& integral);

/// Generates compute function  name.
/// @param integral The base diagonal four center integral.
/// @return The compute function name.
std::pair<size_t, std::string> diag_compute_func_name(const I4CIntegral& integral);

/// Generates compute function  name.
/// @param integral The base diagonal four center integral.
/// @return The compute function name.
std::pair<size_t, std::string> full_compute_func_name(const I4CIntegral& integral);

/// Gets label of standart integrand.
/// @param integrand the integrand operator.
/// @return The label of standart integrand.
std::string integrand_label(const Operator& integrand);

/// Gets recursion namespace label of standart integral.
/// @param integral The base four center integral.
/// @return The recursion namespace label of standart integral.
std::string namespace_label(const I4CIntegral& integral);

/// Generates primitive diagonal compute function name.
/// @param component the integral component.
/// @param integral The base two center integral.
/// @return The primitive diagonal compute function name.
std::pair<size_t, std::string> prim_diag_compute_func_name(const T4CIntegral& component,
                                                           const I4CIntegral& integral);

/// Generates primitive diagonal compute function name.
/// @param component the integral component.
/// @param integral The base two center integral.
/// @return The primitive diagonal compute function name.
std::pair<size_t, std::string> prim_full_compute_func_name(const T4CIntegral& component,
                                                           const I4CIntegral& integral);

/// Generates primitive file name.
/// @param component the integral component.
/// @param integral The base two center integral.
/// @return The primitive file name.
std::string diag_prim_file_name(const T4CIntegral& component,
                                const I4CIntegral& integral);

/// Generates primitive file name.
/// @param component the integral component.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string full_prim_file_name(const T4CIntegral& component,
                                const I4CIntegral& integral);

/// Gets Boys function order for given integral.
/// @param integral The base two center integral.
/// @return The Boys function order of integral.
int boys_order(const I4CIntegral& integral);

/// Gets recursion factors label.
/// @param rterm The recursion term.
/// @param integral The base four center integral.
/// @param first The position of recursion term in code line.
/// @param diagonal The form of integral: diagonal or full.
std::string get_factor_label(const R4CTerm&     rterm,
                             const I4CIntegral& integral,
                             const bool         first,
                             const bool         diagonal);

/// Checks if factors is needed by recursion group.
/// @param rdist  The recursion distribution.
/// @param label The label of factor to find.
/// @return True if factor is found, False otherwise.
bool find_factor(const R4CDist&     rdist,
                 const std::string& label);

/// Prints debug infor for given recursion expansion.
/// @param rdist The recursion expansion.
void debug_info(const R4CDist& rdist);

} // t4c namespace

#endif /* t4c_utils_hpp */
