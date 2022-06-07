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

#ifndef Fraction_hpp
#define Fraction_hpp

#include <string>

/// Fraction class.
class Fraction
{
    /// Numerator of fraction.
    int _numerator;
    
    /// Denominator of fraction.
    int _denominator;
    
    /// Reduces numerator and denominator into standard form.
    void _reduce();
        
public:
    /// Creates an empy fraction.
    Fraction();
    
    /// Creates a faction from only numerator.
    /// @param numerator The numerator of fraction.
    Fraction(const int numerator);
    
    /// Creates a faction from only nominator.
    /// @param numerator The numerator of fraction.
    /// @param denominator The nominator of fraction.
    Fraction(const int numerator,
             const int denominator);
        
    /// Compares this fraction with other fraction.
    /// @param other The other fraction to compare.
    /// @return true if fractions are equal, false otherwise.
    bool operator==(const Fraction& other) const;
    
    /// Compares this fraction with other fraction.
    /// @param other The other fraction to compare.
    /// @return true if fractions are not equal, false otherwise.
    bool operator!=(const Fraction& other) const;
    
    /// Compares this fraction with other fraction.
    /// @param other The other fraction to compare.
    /// @return true if this fraction is less than other fraction, false otherwise.
    bool operator<(const Fraction& other) const;
    
    /// Adds two fractions into single fraction.
    /// @param other The other fraction to be added.
    /// @return The fraction with two added fractions.
    Fraction operator+(const Fraction& other) const;
    
    /// Substracts two fractions into single fraction.
    /// @param other The other fraction to be added.
    /// @return The fraction with two substracted fractions.
    Fraction operator-(const Fraction& other) const;
    
    /// Multuplies two fractions into single fraction.
    /// @param other The other fraction to be added.
    /// @return The fraction with two multiplied fractions.
    Fraction operator*(const Fraction& other) const;
    
    /// Devides two fractions into single fraction.
    /// @param other The other fraction to be added.
    /// @return The fraction with two divided fractions.
    Fraction operator/(const Fraction& other) const;
    
    /// Returns numerator of fraction.
    /// @return The numerator of fraction.
    int numerator() const;
    
    /// Returns denominator of fraction.
    /// @return The numerator of fraction.
    int denominator() const;
    
    /// Determines if fraction is negative.
    /// @return True if fraction is negative, false otherwise.
    bool is_negative() const;
    
    /// Checks if fraction is not a number.
    /// @return true if this fraction is not a number.
    bool isNan() const;
    
    /// Creates frcaction with absolute value of this fraction.
    /// @return the absolute values of fraction. 
    Fraction abs() const;
    
    /// Creates primitive textual representation of this fraction.
    /// @return The string with primitive textual representation of fraction.
    std::string to_string() const;
    
    /// Creates label of this fraction.
    /// @return The string with label of fraction.
    std::string label() const;
};

#endif /* Fraction_hpp */
