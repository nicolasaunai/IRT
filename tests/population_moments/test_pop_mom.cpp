#include "moments.hpp"
#include "population.hpp"
#include "gridlayout.hpp"
#include "field.hpp"
#include "vecfield.hpp"
#include "boundary_condition.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>
#include <cmath>

double density(double x)
{
    // Placeholder for a function that returns density based on x
    return (std::sin(x)+1.)*.5; // Example value
}

void moments()
{
    std::cout << "Running load particles and moments test...\n";

    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 1;
    auto constexpr nppc                          = 100;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    Field<dimension> N{layout->allocate(Quantity::N), Quantity::N};

    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);
    
    std::vector<Population<1>> populations;
    populations.emplace_back("main", layout);
    for (auto& pop : populations)
    {
        pop.load_particles(nppc, density);
        pop.deposit();
        boundary_condition->fill(pop.flux());
        boundary_condition->fill(pop.density());
    }

    total_density(populations, N);
    bulk_velocity<dimension>(populations, N, V);
    {
        
        std::string filename = "moments.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Vx", V.x.data());
        file.createDataSet("/Vy", V.y.data());
        file.createDataSet("/Vz", V.z.data());
        
        file.createDataSet("/N", N.data());
    }
}

int main()
{
  moments();
}
