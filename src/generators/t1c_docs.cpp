#include "t1c_docs.hpp"

#include "file_stream.hpp"

void
T1CDocuDriver::write_doc_str(      std::ofstream& fstream,
                             const int            angmom,
                             const int            gdrv) const
{
    auto lines = VCodeLines();
        
    lines.push_back({0, 0, 1, "/**"});
        
    lines.push_back({0, 0, 2, _get_compute_str(angmom, gdrv)});
    
    for (const auto& label : _get_vars_str())
    {
        lines.push_back({0, 1, 1, label});
    }
        
    lines.push_back({0, 0, 1, "*/"});
        
    ost::write_code_lines(fstream, lines);
}

std::string
T1CDocuDriver::_get_compute_str(const int angmom,
                                const int gdrv) const
{
    return "Evaluates " + std::to_string(gdrv) + "-th order geometrical derivatives for "  + Tensor(angmom).label() + " type GTOs.";
}

std::vector<std::string>
T1CDocuDriver::_get_vars_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param gto_block the GTOs block.");
    
    vstr.push_back("@param grid_coords_x the vector of Cartesian X coordinates of grid.");
    
    vstr.push_back("@param grid_coords_y the vector of Cartesian Y coordinates of grid.");
    
    vstr.push_back("@param grid_coords_z the vector of Cartesian Z coordinates of grid.");
    
    vstr.push_back("@param gtos_mask the mask for GTOs (1 evaluate, 0 skip).");
    
    vstr.push_back("@return the matrix with GTO values on grid points.");
    
    return vstr;
}








