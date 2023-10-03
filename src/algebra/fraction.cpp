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

#include "fraction.hpp"

#include <numeric>

Fraction::Fraction()

    : _numerator(0)
    
    , _denominator(0)
{
    
}

Fraction::Fraction(const int numerator)

    : _numerator(numerator)

    , _denominator(1)
{
    
}

Fraction::Fraction(const int numerator,
                   const int denominator)

    : _numerator(numerator)

    , _denominator(denominator)
{
    _reduce();
}

bool
Fraction::operator==(const Fraction& other) const
{
    if (this == &other) return true;

    if (_numerator != other._numerator)
    {
        return false;
    }
    else
    {
        return _denominator == other._denominator;
    }
}

bool
Fraction::operator!=(const Fraction& other) const
{
    return !((*this) == other);
}

bool
Fraction::operator<(const Fraction& other) const
{
    auto cdenom = std::lcm(_denominator, other._denominator);
    
    return _numerator * (cdenom / _denominator) < other._numerator * (cdenom / other._denominator);
}

Fraction
Fraction::operator+(const Fraction& other) const
{
    auto cdenom = std::lcm(_denominator, other._denominator);
    
    auto cnumer = _numerator * (cdenom / _denominator) +
                
                   other._numerator * (cdenom / other._denominator);
    
    if (cnumer == 0)
    {
        return Fraction(0);
    }
    else
    {
        return Fraction(cnumer, cdenom);
    }
}

Fraction
Fraction::operator-(const Fraction& other) const
{
    auto cdenom = std::lcm(_denominator, other._denominator);
    
    auto cnumer = _numerator * (cdenom / _denominator) -
                
                  other._numerator * (cdenom / other._denominator);
    
    if (cnumer == 0)
    {
        return Fraction(0);
    }
    else
    {
        return Fraction(cnumer, cdenom);
    }
}

Fraction
Fraction::operator*(const Fraction& other) const
{
    return Fraction(_numerator * other._numerator,
                    _denominator * other._denominator);
}

Fraction
Fraction::operator/(const Fraction& other) const
{
    return Fraction(_numerator * other._denominator,
                    _denominator * other._numerator);
}

int
Fraction::numerator() const
{
    return _numerator;
}

int
Fraction::denominator() const
{
    return _denominator;
}

Fraction
Fraction::abs() const
{
    return Fraction(std::abs(_numerator), denominator());
}

bool
Fraction::is_negative() const
{
    return _numerator < 0;
}

bool
Fraction::isNan() const
{
    return _denominator == 0;
}

void
Fraction::_reduce()
{
    if (!isNan())
    {
        auto divisor = std::gcd(_numerator, _denominator);
            
        if (divisor < 0)  divisor *= -1;
        
        _numerator = _numerator / divisor;
    
        _denominator =  _denominator / divisor;
        
        if (_denominator < 0)
        {
            _numerator = -_numerator;
            
            _denominator = -_denominator;
        }
    }
}

std::string
Fraction::to_string() const
{
    return std::to_string(_numerator) + "/" + std::to_string(_denominator);
}

std::string
Fraction::label(std::string spacer) const
{
    if (_denominator == 1)
    {
        return std::to_string(_numerator) + ".0";
    }
    else
    {
        return std::to_string(_numerator) + ".0 / " + std::to_string(_denominator) + ".0";
    }
}
