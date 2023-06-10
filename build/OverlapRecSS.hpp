#ifndef OverlapRecSS_hpp
#define OverlapRecSS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"

namespace ovlrec { // ovlrec namespace

/**
 Evaluates <S||S>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compOverlapSS(      CSubMatrix* matrix,
              const CGtoBlock&  gto_block,
              const int64_t     bra_first,
              const int64_t     bra_last) -> void;

/**
 Evaluates <S||S>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compOverlapSS(      CSubMatrix* matrix,
              const CGtoBlock&  bra_gto_block,
              const CGtoBlock&  ket_gto_block,
              const int64_t     bra_first,
              const int64_t     bra_last) -> void;

/**
 Evaluates block of primitive <S||S>  integrals.

 @param buffer the integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapSS(      TDoubleArray& buffer,
                       const double        bra_exp,
                       const double        bra_norm,
                       const TPoint3D&     bra_coord,
                       const TDoubleArray& ket_exps,
                       const TDoubleArray& ket_norms,
                       const TDoubleArray& ket_coords_x,
                       const TDoubleArray& ket_coords_y,
                       const TDoubleArray& ket_coords_z,
                       const int64_t       ket_dim) -> void;

} // ovlrec namespace

#endif /* OverlapRecSS_hpp */