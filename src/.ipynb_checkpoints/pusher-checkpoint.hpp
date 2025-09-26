#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP

#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>
#include <memory>
#include <array>
#include <cmath>

template<std::size_t dimension>
class Pusher
{
protected:
    std::shared_ptr<GridLayout<dimension>> layout_;
    double dt_;

public:
    Pusher(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : layout_(layout), dt_(dt) {}

    virtual void operator()(std::vector<Particle<dimension>>& particles,
                            VecField<dimension> const& E,
                            VecField<dimension> const& B) = 0;

    virtual ~Pusher() {}
};

template<std::size_t dimension>
class Boris : public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt} {}

    void operator()(std::vector<Particle<dimension>>& particles,
                    VecField<dimension> const& E,
                    VecField<dimension> const& B) override
    {
        for (auto& p : particles)
        {
            double dt = this->dt_;

            // 1. First half-step position update
            p.position[0] += 0.5 * p.v[0] * dt;

            // 2. Find cell index and remainder for interpolation
            double x_norm = p.position[0] / this->layout_->cell_size(Direction::X);
            int iCell = static_cast<int>(x_norm) + this->layout_->dual_dom_start(Direction::X);
            double reminder = x_norm - std::floor(x_norm);

            // 3. Interpolate E and B fields at particle position
            double Ex = interpolate(E.x, iCell, reminder);
            double Ey = interpolate(E.y, iCell, reminder);
            double Ez = interpolate(E.z, iCell, reminder);

            double Bx = interpolate(B.x, iCell, reminder);
            double By = interpolate(B.y, iCell, reminder);
            double Bz = interpolate(B.z, iCell, reminder);

            // 4. Half acceleration by electric field
            double qmdt2 = 0.5 * dt * (p.charge / p.mass);

            double vx_minus = p.v[0] + qmdt2 * Ex;
            double vy_minus = p.v[1] + qmdt2 * Ey;
            double vz_minus = p.v[2] + qmdt2 * Ez;

            // 5. Rotation due to magnetic field
            double tx = qmdt2 * Bx;
            double ty = qmdt2 * By;
            double tz = qmdt2 * Bz;

            double t2 = tx*tx + ty*ty + tz*tz;

            double sx = 2.0*tx / (1.0 + t2);
            double sy = 2.0*ty / (1.0 + t2);
            double sz = 2.0*tz / (1.0 + t2);

            // v' = v_minus + v_minus x t
            double vpx = vx_minus + (vy_minus*tz - vz_minus*ty);
            double vpy = vy_minus + (vz_minus*tx - vx_minus*tz);
            double vpz = vz_minus + (vx_minus*ty - vy_minus*tx);

            // v_plus = v_minus + v' x s
            double vx_plus = vx_minus + (vpy*sz - vpz*sy);
            double vy_plus = vy_minus + (vpz*sx - vpx*sz);
            double vz_plus = vz_minus + (vpx*sy - vpy*sx);

            // 6. Second half acceleration by electric field
            p.v[0] = vx_plus + qmdt2 * Ex;
            p.v[1] = vy_plus + qmdt2 * Ey;
            p.v[2] = vz_plus + qmdt2 * Ez;

            // 7. Second half-step position update
            p.position[0] += 0.5 * p.v[0] * dt;

        }
    }

private:
    // 1st-order linear interpolation
    double interpolate(Field<dimension> const& field, int iCell, double reminder) const
    {
        // Dual-centered fields are shifted by 0.5 cell
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

#endif
