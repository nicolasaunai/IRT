#ifndef HYBIRT_POPULATION_HPP
#define HYBIRT_POPULATION_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "particle.hpp"

#include <random>
#include <optional>
#include <iostream>
#include <string>
#include <functional>


std::mt19937_64 getRNG(std::optional<std::size_t> const& seed)
{
    if (!seed.has_value())
    {
        std::random_device randSeed;
        std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                               randSeed(), randSeed(), randSeed(), randSeed()};
        return std::mt19937_64(seed_seq);
    }
    return std::mt19937_64(*seed);
}


void maxwellianVelocity(std::array<double, 3> const& V, std::array<double, 3> const& Vth,
                        std::mt19937_64& generator, std::array<double, 3>& partVelocity)
{
    std::normal_distribution<> maxwellX(V[0], Vth[0]);
    std::normal_distribution<> maxwellY(V[1], Vth[1]);
    std::normal_distribution<> maxwellZ(V[2], Vth[2]);

    partVelocity[0] = maxwellX(generator);
    partVelocity[1] = maxwellY(generator);
    partVelocity[2] = maxwellZ(generator);
}




template<std::size_t dimension>
class Population
{
public:
    Population(std::string name, std::shared_ptr<GridLayout<dimension>> grid)
        : m_name{name}
        , m_grid{grid}
        , m_flux{grid, {Quantity::Vx, Quantity::Vy, Quantity::Vz}}
        , m_density(m_grid->allocate(Quantity::N), {Quantity::N})
    {
        if (!grid)
            throw std::runtime_error("GridLayout is null");
    }


    void load_particles(int nppc, auto density)
    {
        static_assert(dimension == 1, "Population only implemented for 1D");
        auto randGen = getRNG(std::nullopt);
        std::array<double, 3> Vth{0.2, 0.2, 0.2}; // thermal velocity in each direction
        std::array<double, 3> V{0.0, 0.0, 0.0};   // bulk velocity

        for (auto iCell = m_grid->dual_dom_start(Direction::X);
             iCell <= m_grid->dual_dom_end(Direction::X); ++iCell)
        {
            auto const x      = m_grid->cell_coordinate(Direction::X, iCell);
            auto cell_density = density(x);

            auto cell_weight = cell_density / nppc;
            for (auto partIdx = 0; partIdx < nppc; ++partIdx)
            {
                Particle<1> particle;
                particle.position[0]
                    = x + 0.0 * m_grid->cell_size(Direction::X); // center of the cell
                maxwellianVelocity(V, Vth, randGen, particle.v);
                particle.weight = cell_weight;
                particle.mass   = 1.0; // hard-coded mass
                particle.charge = 1.0; // hard-coded charge

                m_particles.push_back(particle);
            }
        }
        std::cout << "Loaded " << m_particles.size() << " particles.\n";
    }

    void deposit()
    {
        static_assert(dimension == 1, "Population only implemented for 1D");
        for (auto& n : m_density)
        {
            n = 0.0; // Reset the field
        }

        for (auto& fx : m_flux.x)
            fx = 0.0;
        for (auto& fy : m_flux.y)
            fy = 0.0;
        for (auto& fz : m_flux.z)
            fz = 0.0;

        for (auto const& particle : m_particles)
        {
            double const dx    = m_grid->cell_size(Direction::X);
            double const x     = std::fmod(std::fmod(particle.position[0], dx * m_grid->nbr_cells(Direction::X)) + dx * m_grid->nbr_cells(Direction::X),
                                           dx * m_grid->nbr_cells(Direction::X));
            double       s     = x / dx;                              // in [0, Nx)
            int          i0    = static_cast<int>(std::floor(s));     // left node (primal)
            double       w1    = s - i0;                              // weight to the right
            double       w0    = 1.0 - w1;
            
            int left  = i0 + m_grid->dual_dom_start(Direction::X);
            int right = left + 1;
            
            // wrap the right index onto the next cell (periodic)
            if (right > m_grid->dual_dom_end(Direction::X)) {
                right = m_grid->dual_dom_start(Direction::X);
            }
            
            m_density(left)  += particle.weight * w0;
            m_density(right) += particle.weight * w1;
            
            m_flux.x(left)   += particle.weight * particle.v[0] * w0;
            m_flux.x(right)  += particle.weight * particle.v[0] * w1;
            m_flux.y(left)   += particle.weight * particle.v[1] * w0;
            m_flux.y(right)  += particle.weight * particle.v[1] * w1;
            m_flux.z(left)   += particle.weight * particle.v[2] * w0;
            m_flux.z(right)  += particle.weight * particle.v[2] * w1;
        }
    }

    auto& density() { return m_density; }
    auto const& density() const { return m_density; }

    auto& flux() { return m_flux; }
    auto const& flux() const { return m_flux; }

    auto& particles() { return m_particles; }
    auto const& particles() const { return m_particles; }

    auto name() const { return m_name; }

private:
    std::string m_name;
    std::shared_ptr<GridLayout<dimension>> m_grid;
    VecField<dimension> m_flux;
    Field<dimension> m_density;
    std::vector<Particle<dimension>> m_particles;
};

#endif
