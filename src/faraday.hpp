#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

// Get updated B_n+1 from E 

template<std::size_t dimension>
class Faraday
{
    // TODO implement the Faraday class, hint - get inspiration from Ampere
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension> const& B, VecField<dimension>& Bnew)
    {
        // Ex is dual in x
        // Ey is primal in x, so is Ez
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            auto const& Ex = E.x;
            auto const& Ey = E.y;
            auto const& Ez = E.z;
            
            auto const& Bx = B.x;
            auto const& By = B.y;
            auto const& Bz = B.z;

            auto& Bnew_x = Bnew.x;
            auto& Bnew_y = Bnew.y;
            auto& Bnew_z = Bnew.z;
            
            // Bx is in primal
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                Bnew_x(ix) = Bx(ix); // Bx is constant in 1D
            }

            // By, Bz are in dual
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                Bnew_y(ix) = By(ix) + m_dt* (Ez(ix+1) - Ez(ix))/ (dx);
                Bnew_z(ix) = Bz(ix) - m_dt* (Ey(ix+1) - Ey(ix))/ (dx);
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
