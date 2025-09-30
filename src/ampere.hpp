#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    Ampere(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // We calculate the current density using amperes law \mu_0 \mathbf{j} = \nabla \times \mathbf{B}, yee layout
            double mu0 = 1.0; // assuming normalised units. (keywords for grep): \mu \mu_0 mu permeability of vacuum 
            
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X)-1; ++ix)
            {
                J.x(ix) = 0; // in 1D (along x), Jx is zero
                // Jy (dual) -dBz/dx
                J.y(ix) = -(B.z(ix+1) - B.z(ix))/dx / mu0;
                // Jz (dual) dBy/dx
                J.z(ix) = (B.y(ix+1) - B.y(ix))/dx / mu0;
            }
        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
