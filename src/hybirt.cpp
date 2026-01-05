#include "vecfield.hpp"
#include "field.hpp"

#include "faraday.hpp"
#include "ampere.hpp"
#include "ohm.hpp"
#include "utils.hpp"
#include "gridlayout.hpp"
#include "boundary_condition.hpp"
#include "moments.hpp"
#include "pusher.hpp"
#include "diagnostics.hpp"
#include "population.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <memory>
#include <algorithm>


template<std::size_t dimension>
void average(Field<dimension> const& F1, Field<dimension> const& F2, Field<dimension>& Favg)
{
    // use std::transform to do an average of F1 and F2
    std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(), 
        [](double a, double b) {return 0.5*(a+b);}) ;
    // std::transform is used to clearly, safely, and efficiently apply a function element-wise to one or more containers without manual loops.
}


template<std::size_t dimension>
void average(VecField<dimension> const& V1, VecField<dimension> const& V2,
             VecField<dimension>& Vavg)
{
    average(V1.x, V2.x, Vavg.x);
    average(V1.y, V2.y, Vavg.y);
    average(V1.z, V2.z, Vavg.z);
}

double bx(double x)
{
    // Placeholder for a function that returns Bx based on x
    return 0.0; // Example value
}

double by(double x)
{
    // Placeholder for a function that returns By based on x
    return 1.0; // Example value
}


double bz(double x)
{
    // Placeholder for a function that returns Bz based on x
    return 0.0; // Example value
}

double density(double x)
{
    // Placeholder for a function that returns density based on x
    return 1.0; // Example value
}


void magnetic_init(VecField<1>& B, GridLayout<1> const& layout)
{
    // Initialize magnetic field B
    for (auto ix = layout.primal_dom_start(Direction::X); ix <= layout.primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout.coordinate(Direction::X, Quantity::Bx, ix);

        B.x(ix) = bx(x); // Bx
    }
    for (auto ix = layout.dual_dom_start(Direction::X); ix <= layout.dual_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout.coordinate(Direction::X, Quantity::By, ix);

        B.y(ix) = by(x); // By
        B.z(ix) = bz(x); // Bz, uniform magnetic field in z-direction
    }
}


// Compile and run this to check ... 

int main()
{
    double time                     = 0.;
    double final_time               = 10.0000;
    double dt                       = 0.001;
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 1;
    auto constexpr nppc                          = 100;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Enew{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Eavg{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bavg{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    VecField<dimension> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    Field<dimension> N{layout->allocate(Quantity::N), Quantity::N};

    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);

    std::vector<Population<1>> populations;
    populations.emplace_back("main", layout);
    for (auto& pop : populations)
        pop.load_particles(nppc, density);


    magnetic_init(B, *layout);
    boundary_condition->fill(B);

    Faraday<dimension> faraday{layout, dt};  
    Ampere<dimension> ampere{layout};
    Ohm<dimension> ohm{layout};
    Boris<dimension> push{layout, dt};


    ampere(B, J);
    boundary_condition->fill(J);
    for (auto& pop : populations)
    {
        pop.deposit();
        boundary_condition->fill(pop.flux());
        boundary_condition->fill(pop.density());
    }

    total_density(populations, N);
    bulk_velocity<dimension>(populations, N, V);
    
    ohm(B, J, N, V, E);
    boundary_condition->fill(E);

    diags_write_fields(B, E, V, N, time, HighFive::File::Truncate);
    diags_write_particles(populations, time, HighFive::File::Truncate);

    // ICN with 2 predictor steps
    while (time < final_time)
    {
        std::cout << "Time: " << time << " / " << final_time << "\n";
        // 1) First Prediction
        // --------------------

        // Predictor fields (using old moments)
        faraday(E, B, Bnew);
        boundary_condition->fill(Bnew);
    
        ampere(Bnew, J);
        boundary_condition->fill(J);
    
        ohm(Bnew, J, N, V, Enew);
        boundary_condition->fill(Enew);
    
        // Midpoint average fields
        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
    
        // Push particles with the midpoint fields (1st half time step)
        for (auto& pop : populations) {
            push(pop.particles(), Eavg, Bavg);
            boundary_condition->particles(pop.particles());
        }
    
        // Recompute the moments
        for (auto& pop : populations)
        {
            pop.deposit();
            boundary_condition->fill(pop.flux());
            boundary_condition->fill(pop.density());
        }
    
        total_density(populations, N);
        bulk_velocity<dimension>(populations, N, V);

        // 2) Second Predictor
        // --------------------

        faraday(Eavg, B, Bnew);
        boundary_condition->fill(Bnew);

        ampere(Bnew, J); 
        boundary_condition->fill(J);

        ohm(Bnew, J, N, V, Enew);
        boundary_condition->fill(Enew);

        // compute new midpoint
        average(E, Enew, Eavg);        
        average(B, Bnew, Bavg);

        // Push particles agains (2nd half time step)
        for (auto& pop : populations) {
            push(pop.particles(), Eavg, Bavg);
            boundary_condition->particles(pop.particles());
        }
    
        for (auto& pop : populations)
        {
            pop.deposit();
            boundary_condition->fill(pop.flux());
            boundary_condition->fill(pop.density());
        }

        total_density(populations, N);
        bulk_velocity<dimension>(populations, N, V);


        // Corrector fields 
        faraday(Eavg, B, Bnew);
        boundary_condition->fill(Bnew);
    
        ampere(Bnew, J);
        boundary_condition->fill(J);
    
        ohm(Bnew, J, N, V, Enew);
        boundary_condition->fill(Enew);
    
        // Advance the timestep
        E = Enew;
        B = Bnew;
    
        time += dt;
        diags_write_fields(B, E, V, N, time);
        std::cout << "**********************************\n";
        diags_write_particles(populations, time);
    }

    return 0;
}
