#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>


template<std::size_t dimension>
class Pusher
{
protected:
    std::shared_ptr<GridLayout<dimension>> layout_;
    double dt_;

public:
    Pusher(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : layout_(layout)
        , dt_(dt)
    {
    }

    virtual void operator()(std::vector<Particle<dimension>>& particles,
                            VecField<dimension> const& E, VecField<dimension> const& B)
        = 0;

    virtual ~Pusher() {}
};



template<std::size_t dimension>
class Boris : public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt}
    {
    }

void operator()(std::vector<Particle<dimension>>& particles,
                VecField<dimension> const& E,
                VecField<dimension> const& B) override
{
    // --- CORRECTED LINE ---
    // Calculate domain width using functions that exist in GridLayout
    double const number_of_cells = (this->layout_->dual_dom_end(Direction::X) - this->layout_->dual_dom_start(Direction::X) + 1);
    double const domain_width = number_of_cells * this->layout_->cell_size(Direction::X);

    for (auto& particle : particles)
    {
        // Leapfrog: half-step position update
        particle.position[0] += 0.5 * particle.v[0] * this->dt_;
        
        // Compute grid index and remainder for interpolation
        double dx = this->layout_->cell_size(Direction::X);
        int iCell = static_cast<int>(particle.position[0] / dx);
        double rem = (particle.position[0] / dx) - iCell;

        // Interpolate electric and magnetic fields at particle position
        double Ex = this->interpolate(E.x, iCell, rem);
        double Ey = this->interpolate(E.y, iCell, rem);
        double Ez = this->interpolate(E.z, iCell, rem);

        double Bx = this->interpolate(B.x, iCell, rem);
        double By = this->interpolate(B.y, iCell, rem);
        double Bz = this->interpolate(B.z, iCell, rem);

        // Charge and mass
        double q = particle.charge;
        double m = particle.mass;

        // Half electric field acceleration: v_minus = v + (q E dt / 2m)
        double vx_minus = particle.v[0] + 0.5 * q * Ex * this->dt_ / m;
        double vy_minus = particle.v[1] + 0.5 * q * Ey * this->dt_ / m;
        double vz_minus = particle.v[2] + 0.5 * q * Ez * this->dt_ / m;

        // Boris rotation
        double tx = 0.5 * q * Bx * this->dt_ / m;
        double ty = 0.5 * q * By * this->dt_ / m;
        double tz = 0.5 * q * Bz * this->dt_ / m;

        double t2 = tx * tx + ty * ty + tz * tz;
        double sx = 2.0 * tx / (1.0 + t2);
        double sy = 2.0 * ty / (1.0 + t2);
        double sz = 2.0 * tz / (1.0 + t2);

        // v'
        double vpx = vx_minus + (vy_minus * tz - vz_minus * ty);
        double vpy = vy_minus + (vz_minus * tx - vx_minus * tz);
        double vpz = vz_minus + (vx_minus * ty - vy_minus * tx);

        // v+
        double vx_plus = vx_minus + (vpy * sz - vpz * sy);
        double vy_plus = vy_minus + (vpz * sx - vpx * sz);
        double vz_plus = vz_minus + (vpx * sy - vpy * sx);

        // Final half acceleration from E field
        particle.v[0] = vx_plus + 0.5 * q * Ex * this->dt_ / m;
        particle.v[1] = vy_plus + 0.5 * q * Ey * this->dt_ / m;
        particle.v[2] = vz_plus + 0.5 * q * Ez * this->dt_ / m;

        // Final half-step position update
        particle.position[0] += 0.5 * particle.v[0] * this->dt_;

        // APPLY PERIODIC PARTICLE BOUNDARY CONDITIONS
        if (particle.position[0] >= domain_width)
        {
            particle.position[0] -= domain_width;
        }
        else if (particle.position[0] < 0.0)
        {
            particle.position[0] += domain_width;
        }
    }
}

private:
    double interpolate(Field<dimension> const& field, int iCell, double reminder) const
    {
        if (this->layout_->centerings(field.quantity())[0] == this->layout_->dual)
        {
            if (reminder < 0.5)
                return field(iCell - 1) * (1.0 - reminder) + field(iCell) * reminder;
            else
                return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
        }
        return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
    }
};


#endif // PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
