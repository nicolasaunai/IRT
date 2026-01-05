#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday // TODO implement the Faraday class, hint - get inspiration from Ampere
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
            
            auto N = B.y.data().size();
            for (int ix=1; ix< N -1 ; ++ix)
            {
                Bnew.x(ix) = B.x(ix);
                Bnew.y(ix) = B.y(ix) - dt * (-(E.z(ix+1) - E.z(ix-1)) / (2.0 * dx));
                Bnew.z(ix) = B.z(ix) - dt * ((E.y(ix+1) - E.y(ix-1)) / (2.0 * dx));
            }

            // Boundary Conditions
            Bnew.x(0) = B.x(0);
            Bnew.y(0) = B.y(0);
            Bnew.z(0) = B.z(0);
            Bnew.x(N-1) = B.x(N-1);
            Bnew.y(N-1) = B.y(N-1);
            Bnew.z(N-1) = B.z(N-1);   
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};  

#endif // HYBRIDIR_FARADAY_HPP
