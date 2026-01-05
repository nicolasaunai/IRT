#include <iostream>
#include <vector>
#include "highfive/highfive.hpp"
#include "field.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "ampere.hpp"
#include "faraday.hpp"
#include <cmath>

void test_ampere()
{   
    std::cout << "Running ampère test...";
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    for (auto ix = layout->ghost_start(Quantity::By, Direction::X); ix <= layout->ghost_end(Quantity::By, Direction::X); ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, static_cast<std::size_t>(ix));
        B.y(ix) = cos(x);
        B.z(ix) = sin(x);
    }

    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    
    Ampere<dimension> ampere{layout};
    ampere(B, J);

    std::vector<double> expected_Jy;
    std::vector<double> expected_Jz;
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Jy, ix);
        expected_Jy.push_back(- cos(x));
        expected_Jz.push_back(- sin(x));
        double delta_Jy = std::abs((J.y(ix) - expected_Jy.back()) / expected_Jy.back());
        double delta_Jz = std::abs((J.z(ix) - expected_Jz.back()) / expected_Jz.back());
        if (delta_Jy > 1e-2 || delta_Jz > 1e-2)
        {
            std::cout << "Ampere test failed at ix=" << ix << " x=" << x << " Jy=" << J.y(ix)
                      << " Jz=" << J.z(ix) << " delta_Jy=" << delta_Jy
                      << " delta_Jz=" << delta_Jz << "\n";
            return;
        }
    }
    {
        std::string filename = "ampere.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Jx", J.x.data());
        file.createDataSet("/Jy", J.y.data());
        file.createDataSet("/Jz", J.z.data());
        file.createDataSet("/expected_Jy", J.y.data());
        file.createDataSet("/expected_Jz", J.z.data());
    }
    std::cout << "Ampere test passed!\nOutput written to ampere.h5\n";
}

void test_faraday()
{
    std::cout << "Running faraday test...";
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    double dt                                    = 0.1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    for (auto ix = layout->ghost_start(Quantity::Bx, Direction::X); ix <= layout->ghost_end(Quantity::Bx,Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, static_cast<std::size_t>(ix));
        B.x(ix) = exp(-x);
        E.y(ix) = sin(x);
        E.z(ix) = cos(x);
    }
    for (auto ix = layout->ghost_start(Quantity::By, Direction::X); ix <= layout->ghost_end(Quantity::By, Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, static_cast<std::size_t>(ix));
        B.y(ix) = 0.0;
        B.z(ix) = 0.5; 
        
    }

    Faraday<dimension> faraday{layout, dt};
    faraday(B, Bnew, E);

    std::vector<double> expected_Bx;
    std::vector<double> expected_By;
    std::vector<double> expected_Bz;
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, ix);
        expected_Bx.push_back(exp(-x));
        double delta_Bx = std::abs(Bnew.x(ix) - exp(-x));
        if (delta_Bx > 1e-2)
        {
            std::cout << "Faraday test failed at ix=" << ix << " x=" << x << " Bx=" << Bnew.x(ix)
                      << " delta_Bx=" << delta_Bx << "\n";
            return;
        }
    }
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, ix);
        expected_By.push_back(- dt * sin(x));
        expected_Bz.push_back(0.5 - dt * cos(x));
        double delta_By    = std::abs((Bnew.y(ix) - expected_By.back())/expected_By.back());
        double delta_Bz    = std::abs((Bnew.z(ix) - expected_Bz.back())/expected_Bz.back());
        if (delta_By > 1e-3 || delta_Bz > 1e-3)
        {
            std::cout << "Faraday test failed at ix=" << ix << " x=" << x << " By=" << Bnew.y(ix)
                      << " Bz=" << Bnew.z(ix) << " delta_By=" << delta_By
                      << " delta_Bz=" << delta_Bz << "\n";
            return;
        }
    }
    {
        std::string filename = "faraday.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Bnewx", Bnew.x.data());
        file.createDataSet("/Bnewy", Bnew.y.data());
        file.createDataSet("/Bnewz", Bnew.z.data());
        file.createDataSet("/expected_Bx", expected_Bx);
        file.createDataSet("/expected_By", expected_By);
        file.createDataSet("/expected_Bz", expected_Bz);
    }
    std::cout << "Faraday test passed!\nOutput written to faraday.h5\n";
}

int main()
{
    test_ampere();
    test_faraday();
    return 0;
}