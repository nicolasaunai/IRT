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
class Boris: public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt}
    {
    }

    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E, VecField<dimension> const& B) override
    {
        
        
        for (auto& particle: particles)
        {
            double dt = this->dt_;
        
            
            // TODO implement the Boris pusher
            particle.position[0] += particle.v[0] * 0.5 * dt;
            double dx = this->layout_->cell_size(Direction::X) ;
            
            int iCell = static_cast<int>(particle.position[0]/dx);
            double remainder = particle.position[0]/dx - iCell;

            // Electric Fields on the Particle
            double E_x = interpolate(E.x, iCell, remainder);
            double E_y = interpolate(E.y, iCell, remainder);
            double E_z = interpolate(E.z, iCell, remainder);

            // Magnetic Fields on the Particle
            double B_x = interpolate(B.x, iCell, remainder);
            double B_y = interpolate(B.y, iCell, remainder);
            double B_z = interpolate(B.z, iCell, remainder);
            
            
            double qdtm = (particle.charge * dt) / (2 * particle.mass);

            // Half electric Acceleration
            double v_x_minus = particle.v[0] + (qdtm * E_x); 
            double v_y_minus = particle.v[1] + (qdtm * E_y); 
            double v_z_minus = particle.v[2] + (qdtm * E_z);

            // Time
            double t_x = qdtm * B_x;
            double t_y = qdtm * B_y;
            double t_z = qdtm * B_z;


            double v_x_prime = v_x_minus + ((v_y_minus * t_z) - (v_z_minus * t_y));
            double v_y_prime = v_y_minus + ((v_z_minus * t_x) - (v_x_minus * t_z));
            double v_z_prime = v_z_minus + ((v_x_minus * t_y) - (v_y_minus * t_x));

            // Magnetic Rotation
            double t_square = t_x * t_x + t_y * t_y + t_z * t_z;

            double s_x = (2 * t_x)/(1 + t_square);
            double s_y = (2 * t_y)/(1 + t_square);
            double s_z = (2 * t_z)/(1 + t_square);


            double v_x_plus = v_x_minus + ((v_y_prime * s_z) - (v_z_prime * s_y));
            double v_y_plus = v_y_minus + ((v_z_prime * s_x) - (v_x_prime * s_z));
            double v_z_plus = v_z_minus + ((v_x_prime * s_y) - (v_y_prime * s_x));

            // Half 
            particle.v[0] = v_x_plus + qdtm * E_x;
            particle.v[1] = v_y_plus + qdtm * E_y;
            particle.v[2] = v_z_plus + qdtm * E_z;

            // Final Position
            particle.position[0] += particle.v[0] * dt/2;


        };
    }

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


#endif
