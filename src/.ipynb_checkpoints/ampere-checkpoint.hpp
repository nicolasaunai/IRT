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
            //I think Jx must be primal, although always zero
            //I think Jy and Jz must both be primal
            //I should ask about ghost cells, because if not what happens when ix=first with ix-1 ??
            
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto& Jx = J.x;
                Jx(ix) = 0;
                
                auto& Jy = J.y;
                Jy(ix) = - (Bz(ix) - Bz(ix-1))/dx;

                auto& Jz = J.z;
                Jz(ix) = (By(ix) - By(ix-1))/dx;
            } 
        }    
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
