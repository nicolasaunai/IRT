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
            const double dt = this->dt_;
            particle.position[0] += particle.v[0] * dt/2.0; 

            const double dx = this->layout_->cell_size((Direction::X));
            double x_over_dx = particle.position[0]/dx;
            int iCell = static_cast<int>(x_over_dx + this->layout_->dual_dom_start(Direction::X));
            double reminder = (x_over_dx + this->layout_->dual_dom_start(Direction::X)) - iCell;

            double Ex_interp = interpolate(E.x, iCell, reminder);
            double Ey_interp = interpolate(E.y, iCell, reminder);
            double Ez_interp = interpolate(E.z, iCell, reminder);
            double Bx_interp = interpolate(B.x, iCell, reminder);
            double By_interp = interpolate(B.y, iCell, reminder);
            double Bz_interp = interpolate(B.z, iCell, reminder);

            std::array<double, 3> v_minus;
            double qdt_over_2m = (particle.charge * dt) / (2.0 * particle.mass);
            v_minus[0] = particle.v[0] + Ex_interp * qdt_over_2m;
            v_minus[1] = particle.v[1] + Ey_interp * qdt_over_2m;
            v_minus[2] = particle.v[2] + Ez_interp * qdt_over_2m;

            std::array<double, 3> vector_t;
            vector_t[0] = Bx_interp * qdt_over_2m;
            vector_t[1] = By_interp * qdt_over_2m;
            vector_t[2] = Bz_interp * qdt_over_2m;

            std::array<double, 3> v_prime;
            v_prime[0] = v_minus[0] + ( v_minus[1] * vector_t[2] - v_minus[2] * vector_t[1] );
            v_prime[1] = v_minus[1] + ( v_minus[2] * vector_t[0] - v_minus[0] * vector_t[2] );
            v_prime[2] = v_minus[2] + ( v_minus[0] * vector_t[1] - v_minus[1] * vector_t[0] );

            std::array<double, 3> vector_s;
            double vector_t_squared = vector_t[0] * vector_t[0] + vector_t[1] * vector_t[1] + vector_t[2] * vector_t[2];
            double vector_s_prefactor = 2.0 / ( 1.0 + vector_t_squared );
            vector_s[0] = vector_s_prefactor * vector_t[0];
            vector_s[1] = vector_s_prefactor * vector_t[1];
            vector_s[2] = vector_s_prefactor * vector_t[2];

            std::array<double, 3> v_plus;
            v_plus[0] = v_minus[0] + ( v_prime[1] * vector_s[2] - v_prime[2] * vector_s[1] );
            v_plus[1] = v_minus[1] + ( v_prime[2] * vector_s[0] - v_prime[0] * vector_s[2] );
            v_plus[2] = v_minus[2] + ( v_prime[0] * vector_s[1] - v_prime[1] * vector_s[0] );

            particle.v[0] = v_plus[0] + Ex_interp * qdt_over_2m;
            particle.v[1] = v_plus[1] + Ey_interp * qdt_over_2m;
            particle.v[2] = v_plus[2] + Ez_interp * qdt_over_2m;

            particle.position[0] += particle.v[0] * dt/2.0; 
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
