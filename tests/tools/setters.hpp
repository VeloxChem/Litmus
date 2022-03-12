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

#ifndef setters_hpp
#define setters_hpp

#include <string>
#include <vector>

#include "tensor.hpp"
#include "tensor_component.hpp"
#include "two_center_pair_component.hpp"

namespace gset
{

/// Creates vector of tensor components for tensor of given order.
/// @param order The order of tensor
/// @return: A vector of tensor components.
VTensorComponents
tensor_components(const int order);

/// Creates vector of two center pair components for two center pair.
/// @param f_name The name of first expansion center.
/// @param f_angmom The angular momentum of first expansion centers.
/// @param s_name The name of second expansion center.
/// @param s_angmom The angular momentum of second expansion centers.
/// @return: A vector of two center pair components.
VTwoCenterPairComponents
two_center_pair_components(const std::string& f_name,
                           const int          f_angmom,
                           const std::string& s_name,
                           const int          s_angmom);

}

#endif /* setters_hpp */
