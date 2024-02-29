#ifndef cold_npot_driver_hpp
#define cold_npot_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Two center nuclear potential integrals driver class.
class ColdNuclearPotentialDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a two center nuclear potential integrals driver.
    ColdNuclearPotentialDriver();
    
    /// Check if recursion term is for two-center nuclear potential integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to nuclear potential recursion, False otherwise.
    bool is_nuclear_potential(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> bra_vrr(const R2CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies vertical recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> ket_vrr(const R2CTerm& rterm,
                                   const char     axis) const;
        
    /// Applies vertical recursion to bra side recursion term containing nuclear potential
    /// integral.
    /// @param rterm The recursion term with nuclear potential integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_bra_vrr(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to ket side recursion term containing nuclear potential
    /// integral.
    /// @param rterm The recursion term with nuclear potential integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_ket_vrr(const R2CTerm& rterm) const;
        
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_recursion(R2CDist& rdist) const;
    
    /// Recursively applies vertical recursion to bra side of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_bra_vrr(R2CDist& rdist) const;
    
    /// Recursively applies vertical recursion to ket side of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_ket_vrr(R2CDist& rdist) const;
    
    /// Creates recursion group from vector of nuclear potential integral components.
    /// @param vints The  vector of nuclear potnetial integral components.
    /// @return The recursion group.
    R2Group create_recursion(const VT2CIntegrals& vints) const;
    
    /// Recursively applies Obara-Saika recursion to recursion group.
    /// @param rgroup The recursion group.
    void apply_recursion(R2Group& rgroup) const;
};

#endif /* cold_npot_driver_hpp */
