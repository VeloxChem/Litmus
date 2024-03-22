#ifndef cold_eri_driver_hpp
#define cold_eri_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t4c_defs.hpp"

/// Four center electron repulsion integrals driver class.
class ColdFullElectronRepulsionDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a four center electron repulsion integrals driver.
    ColdFullElectronRepulsionDriver();
    
    /// Check if recursion term is for four-center electron repulsion integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to electron repulsion recursion, False otherwise.
    bool is_electron_repulsion(const R4CTerm& rterm) const;
    
    /// Applies vertical recursion to bra side center A of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> bra_vrr_a(const R4CTerm& rterm,
                                     const char     axis) const;
    
    /// Applies vertical recursion to bra side center B of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> bra_vrr_b(const R4CTerm& rterm,
                                     const char     axis) const;
    
    /// Applies vertical recursion to ket side center C of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> ket_vrr_c(const R4CTerm& rterm,
                                     const char     axis) const;
    
    /// Applies vertical recursion to ket side center D of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> ket_vrr_d(const R4CTerm& rterm,
                                     const char     axis) const;
    
    /// Applies vertical recursion to bra side center A recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_vrr_a(const R4CTerm& rterm) const;
    
    /// Applies vertical recursion to bra side center B recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_vrr_b(const R4CTerm& rterm) const;
    
    /// Applies vertical recursion to ket side center C recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_ket_vrr_c(const R4CTerm& rterm) const;
    
    /// Applies vertical recursion to ket side center D recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_ket_vrr_d(const R4CTerm& rterm) const;
    
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_recursion(R4CDist& rdist) const;
    
    /// Recursively applies vertical recursion to bra side center A of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_bra_vrr_a(R4CDist& rdist) const;
    
    /// Recursively applies vertical recursion to bra side center B of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_bra_vrr_b(R4CDist& rdist) const;
    
    /// Recursively applies vertical recursion to ket side center C of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_ket_vrr_c(R4CDist& rdist) const;
    
    /// Recursively applies vertical recursion to ket side center D of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_ket_vrr_d(R4CDist& rdist) const;
    
    /// Creates recursion group from vector of electron repulsion integral components.
    /// @param vints The  vector of electron repulsion integral components.
    /// @return The recursion group.
    R4Group create_recursion(const VT4CIntegrals& vints) const;
    
    /// Recursively applies Obara-Saika recursion to recursion group.
    /// @param rgroup The recursion group.
    void apply_recursion(R4Group& rgroup) const;
};

#endif /* cold_eri_driver_hpp */
