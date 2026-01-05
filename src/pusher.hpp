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
            
             
            particle.position[0] += particle.v[0] * 0.5 * this->dt_;
            double dx = this->layout_->cell_size(Direction::X) ; 

            int iCell = (int)(particle.position[0]/dx) ; 
            
            double remainder = particle.position[0]/dx - iCell;

            double Ex = interpolate(E.x, iCell, remainder);
            double Ey = interpolate(E.y, iCell, remainder);
            double Ez = interpolate(E.z, iCell, remainder);
            
            double Bx = interpolate(B.x, iCell, remainder);
            double By = interpolate(B.y, iCell, remainder);
            double Bz = interpolate(B.z, iCell, remainder);
            
            
            double vminus_x = particle.v[0]+ 0.5 * particle.charge * Ex * this->dt_ / particle.mass; 
            double vminus_y = particle.v[1]+ 0.5 * particle.charge * Ey * this->dt_ / particle.mass ; 
            double vminus_z = particle.v[2] + 0.5 * particle.charge * Ez * this->dt_ / particle.mass ; 

            
            
            double tx = 0.5 * particle.charge  * Bx * this->dt_ / particle.mass; 
            double ty = 0.5 * particle.charge  * By * this->dt_ / particle.mass; 
            double tz = 0.5 * particle.charge  * Bz * this->dt_ / particle.mass; 

            double vprime_x = vminus_x + vminus_y * tz - vminus_z * ty ;
            double vprime_y = vminus_y + vminus_z * tx - vminus_x * tz ;
            double vprime_z = vminus_z + vminus_x * ty - vminus_y * tx ;
            
            double sx = 2 * tx / (1 + std::pow(tx,2) + std::pow(ty,2) + std::pow(tz,2)) ;
            double sy = 2 * ty / (1 + std::pow(tx,2) + std::pow(ty,2) + std::pow(tz,2))  ;
            double sz = 2 * tz / (1 + std::pow(tx,2) + std::pow(ty,2) + std::pow(tz,2))  ; 

            double vplus_x = vminus_x + vprime_y * sz - vprime_z * sy ; 
            double vplus_y = vminus_y + vprime_z * sx - vprime_x * sz ;
            double vplus_z = vminus_z + vprime_x * sy - vprime_y * sx ;

            particle.v[0] = vplus_x + 0.5 * particle.charge * Ex * this->dt_ / particle.mass;
            particle.v[1] = vplus_y + 0.5 * particle.charge * Ey * this->dt_ / particle.mass;
            particle.v[2] = vplus_z + 0.5 * particle.charge * Ez * this->dt_ / particle.mass;

            particle.position[0] += particle.v[0] * 0.5 * this->dt_;
            
        }
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
