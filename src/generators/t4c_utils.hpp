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

#include "t4c_defs.hpp"

namespace t4c { // t4c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base two center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I4CIntegral& integral);

/// Gets standart split label of integral.
/// @param integral The base two center integral.
/// @return The standart split label of integral.
std::string integral_split_label(const I4CIntegral& integral);

/// Gets recursion namespace label of standart integral.
/// @param integral The base two center integral.
/// @return The recursion namespace label of standart integral.
std::string namespace_label(const I4CIntegral& integral);

/// Gets label of standart integrand.
/// @param integrand the integrand operator.
/// @return The label of standart integrand.
std::string integrand_label(const Operator& integrand);

/// Generates compute function  name.
/// @param integral The base four center integral.
/// @return The compute function name.
std::string compute_func_name(const I4CIntegral& integral);

/// Generates integral buffer label.
/// @param integral The base two center integral.
/// @return The string with integral label.
std::string get_buffer_label(const I4CIntegral& integral,
                             const std::string& prefix);

/// Generates primitive compute function name.
/// @param integral The base two center integral.
/// @return The primitive compute function name.
std::string prim_compute_func_name(const I4CIntegral& integral);

/// Gets arguments list for primitive vertical recursion function call.
/// @param integral The base two center integral.
SI4CIntegrals get_vrr_integrals(const I4CIntegral& integral);

/// Generates primitive file name.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string prim_file_name(const I4CIntegral& integral);

/// Generates ket horizontal recursion file name.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string ket_hrr_file_name(const I4CIntegral& integral);

/// Generates bra horizontal recursion file name.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string bra_hrr_file_name(const I4CIntegral& integral);

} // t4c namespace

#endif /* t4c_utils_hpp */