#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension> const& B, VecField<dimension>& Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr(dimension==1)
        {
            // We calculate the change in B: \frac{\partial B}{\partial t} = -\nabla \times \E
            for (auto ix = m_grid->primal_dom_start(Direction::X); ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                Bnew.x(ix) = B.x(ix); // in 1D, Bx is constant and primal
            }
            for (auto ix = m_grid->dual_dom_start(Direction::X); ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                Bnew.y(ix) = B.y(ix) + (E.z(ix+1) - E.z(ix))/dx * dt; // dEz/dx (dual)
                Bnew.z(ix) = B.z(ix) - (E.y(ix+1) - E.y(ix))/dx * dt; // -dEy/dx (dual)
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");

    }
private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double dt;
};

#endif // HYBRIDIR_FARADAY_HPP
