#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>


template<std::size_t dimension>
class Pusher //every Pusher knows about the grid it’s working on and the size of the time step it advances particles by.
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

// Pusher<3> myPusher(layout, 0.01); would make a Pusher in 3D with timestep 0.01.


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

            auto dx = this->layout_->cell_size(Direction::X); 
            auto dt = this->dt_;

            particle.position[0] += particle.v[0]*dt/2;
                    
            auto iCell = static_cast<int>(particle.position[0] / dx); // (taken from layout mesh size) 
            auto remainder = (particle.position[0] / dx) - iCell;

            auto Ex = interpolate(E.x, iCell, remainder);
            auto Ey = interpolate(E.y, iCell, remainder);
            auto Ez = interpolate(E.z, iCell, remainder);

            auto Bx = interpolate(B.x, iCell, remainder);
            auto By = interpolate(B.y, iCell, remainder);
            auto Bz = interpolate(B.z, iCell, remainder);

            auto factor = (particle.charge*dt)/(2*particle.mass);
            
            auto vminus_x = particle.v[0] + factor * Ex;
            auto vminus_y = particle.v[1] + factor * Ey;
            auto vminus_z = particle.v[2] + factor * Ez;
            
            auto t_x = factor*Bx;
            auto t_y = factor*By;
            auto t_z = factor*Bz;

            auto vprime_x = vminus_x + (vminus_y * t_z - vminus_z * t_y);
            auto vprime_y = vminus_y - (vminus_x * t_z - vminus_z * t_x);
            auto vprime_z = vminus_z + (vminus_x * t_y - vminus_y * t_x);

            auto t_2 = t_x*t_x + t_y*t_y + t_z*t_z;
            
            auto s_x = 2*t_x/(1+t_2);
            auto s_y = 2*t_y/(1+t_2);
            auto s_z = 2*t_z/(1+t_2);

            auto vplus_x = vminus_x + (vprime_y * s_z - vprime_z * s_y); 
            auto vplus_y = vminus_y - (vprime_x * s_z - vprime_z * s_x);
            auto vplus_z = vminus_z + (vprime_x * s_y - vprime_y * s_x);
            
            particle.v[0] = vplus_x + factor*Ex;
            particle.v[1] = vplus_y + factor*Ey;
            particle.v[2] = vplus_z + factor*Ez;
            
            particle.position[0] += particle.v[0] * (dt/2);

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
