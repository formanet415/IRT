#include <iostream>
#include <vector>
#include "highfive/highfive.hpp"
#include "field.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "ampere.hpp"
#include "faraday.hpp"
#include "boundary_condition.hpp"

void sin_bz()
{
    std::size_t constexpr dimension = 1;
    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);
    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    auto x_start = layout->coordinate(Direction::X, Quantity::Bx, layout->primal_dom_start(Direction::X));
    auto x_end   = layout->coordinate(Direction::X, Quantity::Bx, layout->primal_dom_end(Direction::X));
    auto pi    = 3.14159265358979323846;
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, ix);

        B.z(ix) = std::sin((x-x_start)/(x_end-x_start)*2*pi); // periodic, just like the boundary condition
    }
    // fill boundaries for B
    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);
    boundary_condition->fill(B);
    
    Ampere<dimension> ampere{layout};
    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    ampere(B, J);
    boundary_condition->fill(J);

    // write to file
    std::vector<double> x_data;
    std::vector<double> bz_data;
    std::vector<double> jy_data;
    for (auto ix = layout->primal_dom_start(Direction::X);
                 ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, ix);
        x_data.push_back(x);
        bz_data.push_back(B.z(ix));
        jy_data.push_back(J.y(ix));
    }

    {
        std::string filename = "sin_bz.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/x", x_data);
        file.createDataSet("/bz", bz_data);
        file.createDataSet("/jy", jy_data);
    }


    std::cout << "Wrote sin_bz.h5\n";
}

void sin_ey()
{
    // test faraday
    std::size_t constexpr dimension = 1;
    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);
    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    auto x_start = layout->coordinate(Direction::X, Quantity::Ey, layout->primal_dom_start(Direction::X));
    auto x_end   = layout->coordinate(Direction::X, Quantity::Ey, layout->primal_dom_end(Direction::X));
    auto pi    = 3.14159265358979323846;
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Ey, ix);

        E.y(ix) = std::sin((x-x_start)/(x_end-x_start)*2*pi); // periodic, just like the boundary condition
        B.z(ix) = 0.0;
    }
    // fill boundaries for E
    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);
    boundary_condition->fill(E);
    boundary_condition->fill(B);

    double dt = 0.05;
    Faraday<dimension> faraday{layout, dt};
    
    
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    faraday(E, B, Bnew);
    boundary_condition->fill(Bnew);

    // we set Ey = sin(x) and Faraday should give us dBz/dt = -dEy/dx = -cos(x)
    // predicted amplitude after dt = 0.05 is 0.05 and from the spatial derivative we get a factor 2pi/L = 2pi/100 = 0.0628


    // write to file
    std::vector<double> x_data;
    std::vector<double> ey_data;
    std::vector<double> bz_data;
    for (auto ix = layout->primal_dom_start(Direction::X);
                 ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Ey, ix);
        x_data.push_back(x);
        ey_data.push_back(E.y(ix));
        bz_data.push_back(Bnew.z(ix));
    }

    


    
}

int main()
{
    sin_bz();
    sin_ey();
}