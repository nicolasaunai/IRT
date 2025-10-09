#include "ampere.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

void ampere_test()
{
    std::cout << "Running Ampere test...\n";
    
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};

    std::vector<double> x_primal_vals;
    std::vector<double> x_dual_vals;
    std::vector<double> Bx_vals, By_vals, Bz_vals;
    std::vector<double> Jx_vals, Jy_vals, Jz_vals;

    ////////////////
    
    for (auto ix = layout->ghost_start(Quantity::Bx, Direction::X); ix <= layout->ghost_end(Quantity::Bx, Direction::X); ++ix)
    {
        B.x(ix) = 0.0; //Bx is primal, though I am setting it at zero  
        Bx_vals.push_back(B.x(ix));
    }
    
    for (auto ix = layout->ghost_start(Quantity::By, Direction::X); ix <= layout->ghost_end(Quantity::By, Direction::X); ++ix)
    {
        B.y(ix) = std::cos(layout->coordinate(Direction::X, Quantity::By, ix)); //By is dual  
        By_vals.push_back(B.y(ix));
    }
    for (auto ix = layout->ghost_start(Quantity::Bz, Direction::X); ix <= layout->ghost_end(Quantity::Bz, Direction::X); ++ix)
    {
        B.z(ix) = std::sin(layout->coordinate(Direction::X, Quantity::Bz, ix)); //Bz is dual  
        Bz_vals.push_back(B.z(ix));
    }
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        double x_primal = layout->coordinate(Direction::X, Quantity::Bx, ix); // To plot primal things  
        x_primal_vals.push_back(x_primal);
    }
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        double x_dual = layout->coordinate(Direction::X, Quantity::By, ix); // To plot dual things
        x_dual_vals.push_back(x_dual);
    }
    
    /////////////////
    /*
    for (auto ix = layout->ghost_start(Direction::X); ix <= layout->ghost_end(Direction::X); ++ix)
    {
        double x_primal = layout->coordinate(Direction::X, Quantity::Bx, ix); // To plot primal things
        double x_dual = layout->coordinate(Direction::X, Quantity::By, ix); // To plot dual things
        
        B.x(ix) = 0.0; //Bx is primal, though I am setting it at zero
        B.y(ix) = std::cos(layout->coordinate(Direction::X, Quantity::By, ix)); //By is dual
        B.z(ix) = std::sin(layout->coordinate(Direction::X, Quantity::Bz, ix)); //Bz is dual

        x_primal_vals.push_back(x_primal);
        x_dual_vals.push_back(x_dual);
        Bx_vals.push_back(B.x(ix));
        By_vals.push_back(B.y(ix));
        Bz_vals.push_back(B.z(ix));
    }
    */
    ///////////////////

    Ampere<dimension> ampere(layout);
    ampere(B, J);
    
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        Jx_vals.push_back(J.x(ix)); //Dual (though it will be zero) 
    }    
  
    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        Jy_vals.push_back(J.y(ix)); //Primal
        Jz_vals.push_back(J.z(ix)); //Primal
    }    
    
    {
        std::string filename = "ampere_test.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/xprimal", x_primal_vals);
        file.createDataSet("/xdual", x_dual_vals);
        file.createDataSet("/Bx", Bx_vals);
        file.createDataSet("/By", By_vals);
        file.createDataSet("/Bz", Bz_vals);
        file.createDataSet("/Jx", Jx_vals);
        file.createDataSet("/Jy", Jy_vals);
        file.createDataSet("/Jz", Jz_vals);
    }
}


int main()
{
    ampere_test();
}
