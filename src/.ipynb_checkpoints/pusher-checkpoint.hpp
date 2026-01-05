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

    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                    VecField<dimension> const& B) override
    {
        for (auto& particle : particles)
        {
            // TODO implement the Boris pusher

            //Mid-step Position
            particle.position[0] += particle.v[0]*(this->dt_/2.0);


            //Fields on particle
            double x = particle.position[0];
            double dx = this->layout_->cell_size(Direction::X);  

            int iCell = (int)(x/dx);
            double reminder = x/dx - iCell;

            //fields on particle (this i get)
            std::array<double, 3> E_half;
            E_half[0] = interpolate(E.x, iCell, reminder);
            E_half[1] = interpolate(E.y, iCell, reminder);
            E_half[2] = interpolate(E.z, iCell, reminder);
         
            std::array<double, 3> B_half;
            B_half[0] = interpolate(B.x, iCell, reminder);
            B_half[1] = interpolate(B.y, iCell, reminder);
            B_half[2] = interpolate(B.z, iCell, reminder);

            //Half electric acceleration
            double q_dt_div_2m = (this->dt_*particle.charge)/(2.0*particle.mass);
            std::array<double, 3> v_minus;

            v_minus[0] = particle.v[0] + q_dt_div_2m*E_half[0];
            v_minus[1] = particle.v[1] + q_dt_div_2m*E_half[1];
            v_minus[2] = particle.v[2] + q_dt_div_2m*E_half[2];

            std::array<double, 3> t;
            t[0] = q_dt_div_2m*B_half[0];
            t[1] = q_dt_div_2m*B_half[1];
            t[2] = q_dt_div_2m*B_half[2];

            std::array<double, 3> v_prime;
            v_prime[0] = v_minus[0] + (v_minus[1]*t[2] - v_minus[2]*t[1]);
            v_prime[1] = v_minus[1] + (v_minus[2]*t[0] - v_minus[0]*t[2]);
            v_prime[2] = v_minus[2] + (v_minus[0]*t[1] - v_minus[1]*t[0]);

            //magnetic rotation
            std::array<double, 3> s;
            s[0]= (2.0*t[0])/(1+t[0]*t[0]);
            s[1]= (2.0*t[1])/(1+t[1]*t[1]);
            s[2]= (2.0*t[2])/(1+t[2]*t[2]);

            std::array<double, 3> v_plus;
            v_plus[0] = v_minus[0] + (v_prime[1]*s[2] - v_prime[2]*s[1]);
            v_plus[1] = v_minus[1] + (v_prime[2]*s[0] - v_prime[0]*s[2]);
            v_plus[2] = v_minus[2] + (v_prime[0]*s[1] - v_prime[1]*s[0]);

            //Half electric acceleration
            particle.v[0] = v_plus[0] + q_dt_div_2m*E_half[0];
            particle.v[1] = v_plus[1] + q_dt_div_2m*E_half[1];
            particle.v[2] = v_plus[2] + q_dt_div_2m*E_half[2];

            //Final position
            particle.position[0] += particle.v[0] * (this->dt_/2.0);
            
            
            
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
