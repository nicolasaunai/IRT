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
    Faraday(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }
    void operator()(VecField<dimension> const& B, VecField<dimension> const& E, VecField<dimension>& Bnew, double dt)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            //Bx is constant and primal
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto const& Bx = B.x;
                auto& Bnew_x = Bnew.x;
                Bnew_x(ix) = Bx(ix);
            }
            
            //By and Bz are both dual 
            //I should ask about ghost cells, because if not what happens when ix=last with ix+1??
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto const& Ey = E.y;
                auto const& Ez = E.z;
                
                auto& Bnew_y = Bnew.y;
                auto& Bnew_z = Bnew.z;
                
                Bnew_y(ix) = By(ix) + dt * (Ez(ix+1) - Ez(ix))/dx;
                Bnew_z(ix) = Bz(ix) - dt * (Ey(ix+1) - Ey(ix))/dx;
            }
        }    
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_FARADAY_HPP
