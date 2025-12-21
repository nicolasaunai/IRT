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
            // TODO your code here
            for (int ix=1; ix< B.y.data().size() -1 ; ++ix)
            {
                J.x(ix) = 0; 

                J.y(ix) = -(B.z(ix+1) - B.z(ix-1)) / (2.0*dx);
                J.z(ix) = (B.y(ix+1) - B.y(ix-1)) / (2.0*dx);
            }
            
            // Boundary Conditions
            J.x(0) = J.x(B.y.data().size()-1) = 0;
            J.y(0) = J.y(B.y.data().size()-1) = 0;
            J.z(0) = J.z(B.z.data().size()-1) = 0;
            
        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
