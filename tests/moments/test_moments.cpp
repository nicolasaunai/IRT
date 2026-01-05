#include "highfive/highfive.hpp"
#include "moments.hpp"
#include "population.hpp"
#include "field.hpp"
#include "gridlayout.hpp"


#include <iostream>
#include <vector>



int main()
{
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 1;
    auto constexpr nppc                          = 100;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    // std::vector<Population<1>> populations;
    // populations.emplace_back("main", layout); //Appends a new element to the end of the container. 
    // for (auto& pop : populations)
    //     pop.load_particles(nppc, density);

    Population<1> pop1("ion1", layout);
    Population<1> pop2("ion2", layout);

    for (std::size_t i = 0; i < pop1.density().data().size(); ++i)
    {
        pop1.density()(i) = 2.0;
        pop2.density()(i) = 3.0;

        pop1.flux().x(i) = 4.0;
        pop2.flux().x(i) = 2.0;
    }
    
    Field<dimension> N{layout->allocate(Quantity::N), Quantity::N};
    std::vector<Population<1>> populations = {pop1, pop2};
    total_density<1>(populations, N);

    std::cout << "Total density:\n";
    for (auto n : N.data())
        std::cout << n << " ";
    std::cout << "\n";

    VecField<1> V(layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz});
    bulk_velocity<1>(populations, N, V);

    std::cout << "Bulk velocity (Vx):\n";
    for (auto v : V.x.data())
        std::cout << v << " ";
    std::cout << "\n";

    return 0;
}



