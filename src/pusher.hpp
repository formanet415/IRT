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
            // Question: how does "this" work?
            double const dt = this->dt_;
            double const dx = this->layout_->cell_size(Direction::X);
            
            particle.position[0] += particle.v[0] * dt/2;

            auto iCell = static_cast<int>(particle.position[0]/dx);
            auto remainder = particle.position[0] - iCell;
            
            auto Ex = interpolate(E.x, iCell, remainder);
            auto Ey = interpolate(E.y, iCell, remainder);
            auto Ez = interpolate(E.z, iCell, remainder);

            auto Bx = interpolate(B.x, iCell, remainder);
            auto By = interpolate(B.y, iCell, remainder);
            auto Bz = interpolate(B.z, iCell, remainder);
            
            // Precompute constants
            auto qmdt2 = particle.charge / (2 * particle.mass) * dt;

            auto vminx = particle.v[0] + Ex * qmdt2;
            auto vminy = particle.v[1] + Ey * qmdt2;
            auto vminz = particle.v[2] + Ez * qmdt2;

            auto tx = Bx * qmdt2;
            auto ty = By * qmdt2;
            auto tz = Bz * qmdt2;

            // Boris rotation
            auto vprimex = vminx + (vminy * tz - vminz * ty);
            auto vprimey = vminy + (vminz * tx - vminx * tz);
            auto vprimez = vminz + (vminx * ty - vminy * tx);

            auto t2 = tx*tx + ty*ty + tz*tz;
            auto sx = 2 * tx / (1 + t2);
            auto sy = 2 * ty / (1 + t2);
            auto sz = 2 * tz / (1 + t2);

            auto vplusx = vminx + (vprimey * sz - vprimez * sy);
            auto vplusy = vminy + (vprimez * sx - vprimex * sz);
            auto vplusz = vminz + (vprimex * sy - vprimey * sx);

            particle.v[0] = vplusx + particle.charge*Ex/2/particle.mass*dt;
            particle.v[1] = vplusy + particle.charge*Ey/2/particle.mass*dt;
            particle.v[2] = vplusz + particle.charge*Ez/2/particle.mass*dt;

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
