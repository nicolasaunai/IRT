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
            double const dx = this->layout_->cell_size(Direction::X);
            double const Lx = this->layout_->dom_size(Direction::X);
            double const dt = this->dt_;
 
            particle.position[0] += particle.v[0] * dt * 0.5;

            // Manual periodic wrapping 
            if (particle.position[0] < 0.0)
                particle.position[0] += Lx;
            else if (particle.position[0] >= Lx)
                particle.position[0] -= Lx;
            
            double x_cell = particle.position[0] / dx;
            int iCell = static_cast<int>(std::floor(x_cell)); // + this->layout_->dual_dom_start(Direction::X));
            double remainder = x_cell - iCell;
            
            // Clamp remainder safely to [0,1)
            if (remainder < 0.0) remainder = 0.0;
            else if (remainder >= 1.0) remainder = std::nextafter(1.0, 0.0);

            /*
            double const dx = this->layout_->cell_size(Direction::X);
            double const dt = this->dt_;
            
            particle.position[0] += particle.v[0] * dt * .5;
                
            int const iCell = static_cast<int>(particle.position[0]/dx); // + this->layout_->dual_dom_start(Direction::X)); 
            double remainder = (particle.position[0] / dx) - iCell;
            */
            
            double Ex = this->interpolate(E.x, iCell, remainder);
            double Ey = this->interpolate(E.y, iCell, remainder);
            double Ez = this->interpolate(E.z, iCell, remainder);
    
            double Bx = this->interpolate(B.x, iCell, remainder);
            double By = this->interpolate(B.y, iCell, remainder);
            double Bz = this->interpolate(B.z, iCell, remainder);
    
            double qmdt2 = (particle.charge / particle.mass) * 0.5 * dt;

            double vminus_x = particle.v[0] + qmdt2 * Ex;
            double vminus_y = particle.v[1] + qmdt2 * Ey;
            double vminus_z = particle.v[2] + qmdt2 * Ez;
    
            double tx = qmdt2 * Bx;
            double ty = qmdt2 * By;
            double tz = qmdt2 * Bz;
    
            double t2 = tx*tx + ty*ty + tz*tz;

            // painful vector manipulations

            double vprime_x = vminus_x + (vminus_y*tz - vminus_z*ty);
            double vprime_y = vminus_y + (vminus_z*tx - vminus_x*tz);
            double vprime_z = vminus_z + (vminus_x*ty - vminus_y*tx);
 
            double sx = 2.0 * tx / (1.0 + t2);
            double sy = 2.0 * ty / (1.0 + t2);
            double sz = 2.0 * tz / (1.0 + t2);

            double vplus_x = vminus_x + (vprime_y*sz - vprime_z*sy);
            double vplus_y = vminus_y + (vprime_z*sx - vprime_x*sz);
            double vplus_z = vminus_z + (vprime_x*sy - vprime_y*sx);
    
            particle.v[0] = vplus_x + qmdt2 * Ex;
            particle.v[1] = vplus_y + qmdt2 * Ey;
            particle.v[2] = vplus_z + qmdt2 * Ez;

            particle.position[0] += 0.5 * particle.v[0] * dt;
        }
    }

/*
private:
    double interpolate(Field<dimension> const& field, int iCell, double remainder) const
    {
        if (this->layout_->centerings(field.quantity())[0] == this->layout_->dual)
        {
            if (remainder < 0.5)
                return field(iCell - 1) * (1.0 - remainder) + field(iCell) * remainder;
            else
                return field(iCell) * (1.0 - remainder) + field(iCell + 1) * remainder;
        }
        return field(iCell) * (1.0 - remainder) + field(iCell + 1) * remainder;
    }
};
*/

private:
    double interpolate(Field<dimension> const& field, int iCell, double remainder) const
    {
        if (iCell < 0 || iCell + 1 >= field.size())
        {
            std::cout << "iCell " << iCell << " out of range [0, " <<  field.size() << "].\n";
        }
        if (remainder < 0. || remainder >= 1.)
        {
            std::cout << "Remainder " << remainder << " out of bound.\n";
        }

        const auto gsi = this->layout_->ghost_start(field.quantity(), Direction::X);
        const auto dsi = this->layout_->dom_start  (field.quantity(), Direction::X);
        const auto dei = this->layout_->dom_end    (field.quantity(), Direction::X);
        const auto gei = this->layout_->ghost_end  (field.quantity(), Direction::X);

        /*
        bool isDual = (this->layout_->centerings(field.quantity())[0] == this->layout_->dual);
        */

        int left_index = 0, right_index = 0;

        left_index  = iCell;
        right_index = iCell + 1;
/*      
        // Safety checks: ensure indices are within allocated range [gsi, gei]
        // If they fall outside, it's either a setup bug (ghost cells missing) or
        // a particle outside domain — handle explicitly.
#ifndef NDEBUG
        // In debug mode: fail fast so you catch the exact problem
        if (left_index < static_cast<int>(gsi) || right_index > static_cast<int>(gei)) {
            std::cerr << "Interpolate index OOB: left_index=" << left_index
                      << " right_index=" << right_index
                      << " valid=[" << gsi << "," << gei << "]"
                      << " remainder=" << remainder << " iCell=" << iCell << "\n";
            throw std::out_of_range("Interpolate index out of ghost range");
       }

#else
        // In release mode: clamp to edges of the allocated buffer to avoid crashes.
        // This avoids UB but may degrade correctness near boundaries.
        left_index  = std::max<int>(left_index,  static_cast<int>(gsi));
        right_index = std::min<int>(right_index, static_cast<int>(gei));

#endif
*/
        double f_left  = field(left_index);
        double f_right = field(right_index);

        return f_left * (1.0 - remainder) + f_right * remainder;
    }
};


#endif
