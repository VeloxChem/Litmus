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

namespace t2c { // t2c namespace

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

/// Gets vector of tensor component labels.
/// @param tensor the Cartessian tensor.
/// @param label the base of label.
/// @return The vecotr of tensor component labels.
std::vector<std::string> tensor_components(const Tensor&      tensor,
                                           const std::string& label);

} // t2c namespace

#endif /* t2c_utils_hpp */
