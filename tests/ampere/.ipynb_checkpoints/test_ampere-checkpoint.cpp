#include "ampere.hpp"
#include "gridlayout.hpp"
#include "vecfield.hpp"

#include <iostream>
#include <memory>
#include <cmath>

int main()
{
    constexpr std::size_t dimension = 1;

    std::array<std::size_t, dimension> grid_size = {10};
    std::array<double, dimension> cell_size = {1.0};
    constexpr std::size_t nbr_ghosts = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> B(layout, {Quantity::Bx, Quantity::By, Quantity::Bz});

    for (auto ix = layout->dual_dom_start(Direction::X); 
        ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        B.x(ix) = 0.0;
    }

    for (auto ix = layout->primal_dom_start(Direction::X); 
        ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        double x = ix * cell_size[0];
        B.y(ix) = x;
        B.z(ix) = 2.0 * x;
    }

    VecField<dimension> J(layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz});

    Ampere<dimension> ampere(layout);
    ampere(B, J);

    std::cout << "Jy and Jz from Ampere (expect Jy = -2, Jz = 1):\n";
    for (std::size_t ix = layout->primal_dom_start(Direction::X);
         ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        std::cout << "ix=" << ix
                  << " Jy=" << J.y(ix)
                  << " Jz=" << J.z(ix)
                  << "\n";
    }

    return 0;
}