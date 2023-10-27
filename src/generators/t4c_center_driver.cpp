//
//  t4c_center_driver.cpp
//  Litmus
//
//  Created by Zilvinas Rinkevicius on 2023-10-27.
//

#include "t4c_center_driver.hpp"


#include "axes.hpp"
#include "t2c_utils.hpp"

T4CCenterDriver::T4CCenterDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T4CCenterDriver::is_auxilary(const R4CTerm& rterm,
                             const int      index) const
{
    if (const auto nprefixes = rterm.prefixes().size(); index >= nprefixes)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<R4CDist>
T4CCenterDriver::bra_ket_vrr(const R4CTerm& rterm,
                             const char     axis,
                             const int      index) const
{
    if (is_auxilary(rterm, index)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, index, true))
    {
        R4CDist t4crt(rterm);
        
        if (const auto r1val = tval->shift(axis, 1, index))
        {
            auto x1val = *r1val;
            
            if (index == 0)
            {
                x1val.add(Factor("ba_e", "tba_e"), Fraction(2));
            }
            
            if (index == 1)
            {
                x1val.add(Factor("bb_e", "tbb_e"), Fraction(2));
            }
            
            if (index == 2)
            {
                x1val.add(Factor("kc_e", "tkc_e"), Fraction(2));
            }
            
            if (index == 3)
            {
                x1val.add(Factor("kd_e", "tkd_e"), Fraction(2));
            }
            
            t4crt.add(x1val);
        }
        
        if (const auto r2val = tval->shift(axis, -1, index))
        {
            auto x2val = *r2val;
            
            x2val.scale(Fraction(-(*tval)[index][axis]));
            
            t4crt.add(x2val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
T4CCenterDriver::apply_bra_ket_vrr(const R4CTerm& rterm,
                                   const int      index) const
{
    R4CDist t4crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_ket_vrr(rterm, axis, index))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t4crt;
}

void
T4CCenterDriver::apply_recursion(R4CDist& rdist) const
{
    // vertical recursions on D center
    
    apply_bra_ket_vrr(rdist, 3);
    
    // vertical recursions on C side
    
    apply_bra_ket_vrr(rdist, 2);
    
    // vertical recursions on B side
    
    apply_bra_ket_vrr(rdist, 1);
    
    // vertical recursions on A side
    
    apply_bra_ket_vrr(rdist, 0);
}

void
T4CCenterDriver::apply_bra_ket_vrr(      R4CDist& rdist,
                                   const int      index) const
{
    if (!is_auxilary(rdist.root(), index))
    {
        R4CDist new_dist(rdist.root());
        
        V4CTerms rec_terms;
        
        // set up initial terms for recursion expansion
        
        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                if (const auto rterm = rdist[i]; is_auxilary(rterm, index))
                {
                    new_dist.add(rterm);
                }
                else
                {
                    rec_terms.push_back(rterm);
                }
            }
        }
        else
        {
            if (const auto rterm = rdist.root(); !is_auxilary(rterm, index))
            {
                rec_terms.push_back(rterm);
            }
        }
        
        // apply recursion until only
            
        while (!rec_terms.empty())
        {
            V4CTerms new_terms;
            
            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_bra_ket_vrr(rec_terms[i], index);
                
                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; is_auxilary(rterm, index))
                        {
                            new_dist.add(rterm);
                        }
                        else
                        {
                            new_terms.push_back(rterm);
                        }
                    }
                }
            }
            
            rec_terms = new_terms;
        }
        
        // update recursion distribution
        
        rdist = new_dist;
    }
}

R4Group
T4CCenterDriver::create_recursion(const VT4CIntegrals& vints) const
{
    // create reference group
    
    R4Group r4group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R4CDist(R4CTerm(tcomp));
        
        apply_recursion(rdist);
                
        r4group.add(rdist);
    }
    
    r4group.simplify();
    
    return r4group;
}
