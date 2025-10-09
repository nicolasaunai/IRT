#include "gridlayout.hpp"
#include "population.hpp"
#include "moments.hpp"

#include <iostream>
#include <cmath>
#include <vector>

void test_moments()
{
    std::cout << "Running population + moments test...\n";
    
    constexpr std::size_t dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};   
    std::array<double, dimension> cell_size = {0.01};       
    auto constexpr nbr_ghosts = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    //Population with a uniform density profile
    Population<dimension> electrons("electrons", layout);

    auto density_profile = [](double x) 
    {
        return 1.0; //Note that I chose it to be uniform, but I could have chosen a coordinate dependence
    };

    int nppc = 50; //Number of particles per cell
    electrons.load_particles(nppc, density_profile);
    electrons.deposit();

    //Now I compute total density and bulk velocity
    std::vector<Population<dimension>> populations{electrons};
    Field<dimension> N(layout->allocate(Quantity::N), {Quantity::N});
    VecField<dimension> V(layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz});

    total_density(populations, N);
    bulk_velocity(populations, N, V);

    
    std::cout << "# x   N(x)   Vx(x)   Vy(x)   Vz(x)\n";
    for (int ix = layout->dual_dom_start(Direction::X);
         ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        double x = layout->coordinate(Direction::X, Quantity::Bx, ix); //Taking everything as primal, including posisition
        std::cout << x << "  "
                  << N(ix) << "  "
                  << V.x(ix) << "  "
                  << V.y(ix) << "  "
                  << V.z(ix) << "\n";
    }

    std::cout << "test_moments completed successfully.\n";
}

int main()
{
        test_moments();
}

