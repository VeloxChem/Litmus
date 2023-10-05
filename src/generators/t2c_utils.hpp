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

#include <string>

#include "t2c_defs.hpp"

using TMapOfStrings = std::map<Operator, std::string>; 

namespace t2c { // t2c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base two center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I2CIntegral& integral);

/// Gets label of standart integrand.
/// @param integrand the integrand operator.
/// @return The label of standart integrand.
std::string integrand_label(const Operator& integrand);

/// Gets vector of integrand component labels.
/// @param integrand the integrand operator.
/// @param label the base of label.
/// @return The vecotr of integrand component labels.
std::vector<std::string> integrand_components(const Operator&    integrand,
                                              const std::string& label);

/// Gets vector of integrand component labels.
/// @param bra_tensor the Cartessian tensor.
/// @param integrand the integrand operator.
/// @param label the base of label.
/// @return The vecotr of integrand component labels.
std::vector<std::string> integrand_components(const Tensor&      bra_tensor,
                                              const Operator&    integrand,
                                              const std::string& label);

/// Gets vector of integrand component labels.
/// @param bra_tensor the Cartessian tensor.
/// @param ket_tensor the Cartessian tensor.
/// @param integrand the integrand operator.
/// @param label the base of label.
/// @return The vecotr of integrand component labels.
std::vector<std::string> integrand_components(const Tensor&      bra_tensor,
                                              const Tensor&      ket_tensor,
                                              const Operator&    integrand,
                                              const std::string& label);

/// Gets vector of tensor component labels.
/// @param tensor the Cartessian tensor.
/// @param label the base of label.
/// @return The vecotr of tensor component labels.
std::vector<std::string> tensor_components(const Tensor&      tensor,
                                           const std::string& label);

/// Generates compute function  name.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @return The compute function name.
std::pair<size_t, std::string> compute_func_name(const I2CIntegral& integral,
                                                 const bool         sum_form);

/// Generates primitive compute function name.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @return The primitive compute function name.
std::pair<size_t, std::string> prim_compute_func_name(const I2CIntegral& integral,
                                                      const bool         sum_form);

/// Generates primitive compute function name.
/// @param component the integral component.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @param bra_first The flag to set bra as expansion point.
/// @return The primitive compute function name.
std::pair<size_t, std::string> prim_compute_func_name(const TensorComponent& component,
                                                      const I2CIntegral&     integral,
                                                      const bool             sum_form,
                                                      const bool             bra_first);

/// Generates primitive compute function name.
/// @param bra_component the integral component on bra side.
/// @param ket_component the integral component on ket side.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @return The primitive compute function name.
std::pair<size_t, std::string> prim_compute_func_name(const TensorComponent& bra_component,
                                                      const TensorComponent& ket_component,
                                                      const I2CIntegral&     integral,
                                                      const bool             sum_form);

/// Generates primitive file name.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @return The primitive file name.
std::string prim_file_name(const I2CIntegral& integral,
                           const bool         sum_form);

/// Generates primitive file name.
/// @param component the integral component.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @param bra_first The flag to set bra as expansion point.
/// @return The primitive file name.
std::string prim_file_name(const TensorComponent& component,
                           const I2CIntegral&     integral,
                           const bool             sum_form,
                           const bool             bra_first);

/// Generates primitive file name.
/// @param bra_component the integral component on bra side.
/// @param ket_component the integral component on ket side.
/// @param integral The base two center integral.
/// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
/// @return The primitive file name.
std::string prim_file_name(const TensorComponent& bra_component,
                           const TensorComponent& ket_component,
                           const I2CIntegral&     integral,
                           const bool             sum_form);

/// Gets recursion namespace label of standart integral.
/// @param integral The base two center integral.
/// @return The recursion namespace label of standart integral.
std::string namespace_label(const I2CIntegral& integral);


/// Gets index of tensor componnent in it's component wwise expansion.
/// @param component The tensor component to find index.
/// @return The index of tensor component.
int tensor_component_index(const TensorComponent& component);

/// Combines two symbolic factors into  one.
/// @param bra_factor The  symbolic factor on bra side.
/// @param ket_factor The  symbolic factor on ket side.
/// @return The combined  factor.
std::string combine_factors(const std::string& bra_factor,
                            const std::string& ket_factor);

/// Checks if factors is needed by recursion group.
/// @param rgroup The recursion group.
/// @param label The label of factor to find.
/// @return True if factor is found, False otherwise.
bool find_factor(const R2Group&     rgroup,
                 const std::string& label);

/// Gets recursion factors label.
/// @param rterm The recursion term.
/// @param first The position of recursion term in code line.
std::string get_factor_label(const R2CTerm& rterm,
                             const bool     first);

/// Gets Boys function order for given integral.
/// @param integral The base two center integral.
/// @return The Boys function order of integral.
int boys_order(const I2CIntegral& integral);

/// Checks if Boys function is needed for given integral.
/// @param integral The base two center integral.
/// @return True if Boys function is needed by given integral, False otherwise.
bool need_boys(const I2CIntegral& integral);

/// Prints debug infor for given recursion expansion.
/// @param rdist The recursion expansion.
void debug_info(R2CDist& rdist);

} // t2c namespace

#endif /* t2c_utils_hpp */
