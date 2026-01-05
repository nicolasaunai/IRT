#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>
#include <cmath>


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

    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                    VecField<dimension> const& B) override
    {
        for (auto& particle : particles)
        {
            // TODO implement the Boris pusher
            
            double x;
            double x_half;
            
            double tx;
            double ty;
            double tz;
            
            double sx;
            double sy;
            double sz;
            
            double vx;
            double vx_minus;
            double vx_prime;
            double vx_plus;
            
            double vy;
            double vy_minus;
            double vy_prime;
            double vy_plus;
            
            double vz;
            double vz_minus;
            double vz_prime;
            double vz_plus;


            x = particle.position[0];
            vx = particle.v[0];
            vy = particle.v[1];
            vz = particle.v[2];

            x_half = x + vx*this->dt_/2;

            double const iCell_float = x_half / this->layout_->cell_size(Direction::X)
                                   + this->layout_->dual_dom_start(Direction::X);
            int const iCell = static_cast<int>(iCell_float);
            double reminder = iCell_float - iCell;

            double const Lx = this->layout_->nbr_cells(Direction::X) *
                  this->layout_->cell_size(Direction::X);

            double Ex = interpolate(E.x, iCell, reminder);
            double Ey = interpolate(E.y, iCell, reminder);
            double Ez = interpolate(E.z, iCell, reminder);
            double Bx = interpolate(B.x, iCell, reminder);
            double By = interpolate(B.y, iCell, reminder);
            double Bz = interpolate(B.z, iCell, reminder);

            vx_minus = vx + particle.charge*this->dt_*Ex/(2*particle.mass);
            vy_minus = vy + particle.charge*this->dt_*Ey/(2*particle.mass);
            vz_minus = vz + particle.charge*this->dt_*Ez/(2*particle.mass);

            tx = particle.charge*this->dt_*Bx/(2*particle.mass);
            ty = particle.charge*this->dt_*By/(2*particle.mass);
            tz = particle.charge*this->dt_*Bz/(2*particle.mass);

            vx_prime = vx_minus + vy_minus*tz - vz_minus*ty;
            vy_prime = vy_minus + vz_minus*tx - vx_minus*tz;
            vz_prime = vz_minus + vx_minus*ty - vy_minus*tx;

            sx = 2*tx/(1+(tx*tx +ty*ty + tz*tz));
            sy = 2*ty/(1+(tx*tx +ty*ty + tz*tz));
            sz = 2*tz/(1+(tx*tx +ty*ty + tz*tz));

            vx_plus = vx_minus + vy_prime*sz - vz_prime*sy;
            vy_plus = vy_minus + vz_prime*sx - vx_prime*sz;
            vz_plus = vz_minus + vx_prime*sy - vy_prime*sx;

            vx = vx_plus + particle.charge*this->dt_*Ex/(2*particle.mass);
            vy = vy_plus + particle.charge*this->dt_*Ey/(2*particle.mass);
            vz = vz_plus + particle.charge*this->dt_*Ez/(2*particle.mass);

            x = x_half + vx*this->dt_/2;

            particle.position[0] = x;
            particle.v[0] = vx;
            particle.v[1] = vy;
            particle.v[2] = vz;
            
            
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


#endif
