#include "population.hpp"
#include "moments.hpp"
#include "highfive/highfive.hpp"
#include "field.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"
#include <iostream>
#include <vector>

double density(double x)
{
    // Placeholder for a function that returns density based on x
    return 1.0;
}

void test_population()
{
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 0;
    auto constexpr nppc                          = 10000;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    Population<1> pop("main", layout);
    pop.load_particles(nppc, density);
    pop.deposit();

    std::vector<Population<1>> populations = {pop};
    Field<1> N{layout->allocate(Quantity::N), Quantity::N};
    VecField<1> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};

    total_density(populations, N);
    bulk_velocity(populations, N, V);

    {
        std::string filename = "population.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/density", N.data());
        file.createDataSet("/vx", V.x.data());
        file.createDataSet("/vy", V.y.data());
        file.createDataSet("/vz", V.z.data());
    }
}

int main()
{
    test_population();
    return 0;
}
