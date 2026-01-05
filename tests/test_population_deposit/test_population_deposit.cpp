// test_population_deposit.cpp
#include "population.hpp"
#include "gridlayout.hpp"
#include "moments.hpp"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    constexpr std::size_t dim = 1;
    std::array<std::size_t, dim> grid_size = {10};
    std::array<double, dim> cell_size = {1.0};
    auto layout = std::make_shared<GridLayout<dim>>(grid_size, cell_size, 1);

    Population<dim> pop{"test_species", layout};

    Particle<dim> p;
    p.mass = 1.0;
    p.charge = 1.0;
    p.weight = 1.0;

    for (int i = 0; i < grid_size[0]; ++i) {
        for (int k = 0; k < 2; ++k) {
            p.position[0] = (i + 0.25 * k) * cell_size[0];
            p.v[0] = 1.0; p.v[1] = 0.0; p.v[2] = 0.0;
            pop.particles().push_back(p);
        }
    }

    pop.deposit();

    std::vector<Population<1>> pops = {pop};

    Field<dim> N{layout->allocate(Quantity::N), Quantity::N};
    VecField<dim> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    
    total_density(pops, N);
    bulk_velocity<dim>(pops, N, V);

    double mean_density = 0.0;
    for (auto ix = layout->primal_dom_start(Direction::X);
     ix <= layout->primal_dom_end(Direction::X); ++ix)
        mean_density += N(ix);
    mean_density /= grid_size[0];

    std::cout << "Mean density = " << mean_density << " (expected 2.0)\n";

    double mean_vx = 0.0;
    for (auto ix = layout->primal_dom_start(Direction::X);
     ix <= layout->primal_dom_end(Direction::X); ++ix)
        mean_vx += V.x(ix);
    mean_vx /= grid_size[0];

    std::cout << "Mean bulk Vx = " << mean_vx << " (expected 1.0)\n";

}
