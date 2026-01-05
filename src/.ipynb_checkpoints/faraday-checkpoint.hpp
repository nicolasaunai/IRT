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
    // TODO implement the Faraday class, hint - get inspiration from Ampere
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}
        , m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension>& B, VecField<dimension>& Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto& Bx = Bnew.x;
                Bx(ix) = B.x(ix);
            }
            // By, Bz are dual
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                auto const& Ex = E.x;
                auto const& Ey = E.y;
                auto const& Ez = E.z;

                
                auto& By = Bnew.y;
                auto& Bz = Bnew.z;

                
                By(ix) = B.y(ix) + (Ez(ix + 1) - Ez(ix))*m_dt/dx; 
                Bz(ix) = B.z(ix) - (Ey(ix + 1) - Ey(ix))*m_dt/dx;
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;
};

#endif // HYBRIDIR_FARADAY_HPP
