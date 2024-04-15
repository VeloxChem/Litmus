#ifndef v2i_ovl_driver_hpp
#define v2i_ovl_driver_hpp

#include <optional>
#include <array>

#include "t2c_defs.hpp"

/// Two center overlap integrals driver class.
class V2IOverlapDriver
{
    /// Applies vertical recursion to bra side of overlap integral.
    /// @param integral The  overlap integral.
    /// @return The set of integrals.
    SI2CIntegrals _bra_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side of overlap integral.
    /// @param integral The  overlap integral.
    /// @return The set of integrals.
    SI2CIntegrals _ket_vrr(const I2CIntegral& integral) const;

public:
    /// Creates a two center overlap integrals driver.
    V2IOverlapDriver() = default;
    
    /// Check if integral is for two-center overlap integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to overlap recursion, False otherwise.
    bool is_overlap(const I2CIntegral& integral) const;
    
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param integrals The  integral to apply recursion.
    /// @return The set of integrals.
    SI2CIntegrals apply_recursion(const SI2CIntegrals& integrals) const;
    
    /// Applies vertical recursion to bra side of overlap integral.
    /// @param integral The  overlap integral.
    /// @return The recursion expansion of integral.
    SI2CIntegrals apply_bra_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side of overlap integral.
    /// @param integral The  overlap integral.
    /// @return The recursion expansion of integral.
    SI2CIntegrals apply_ket_vrr(const I2CIntegral& integral) const;
    
    /// Creates recursion expansion for set of integral.
    /// @param integrals The  integral to apply recursion.
    /// @return The set of integrals.
    SI2CIntegrals create_recursion(const SI2CIntegrals& integrals) const;
};

#endif /* v2i_ovl_driver_hpp */
