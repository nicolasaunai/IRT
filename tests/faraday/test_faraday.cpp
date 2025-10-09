#include "faraday.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

void faraday_test()
{
    std::cout << "Running Faraday test...\n";
    
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    double dt = 0.001;

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    std::vector<double> x_primal_vals;
    std::vector<double> x_dual_vals;
    std::vector<double> Ex_vals, Ey_vals, Ez_vals;
    std::vector<double> Bx_vals, By_vals, Bz_vals;
    std::vector<double> Bxnew_vals, Bynew_vals, Bznew_vals;

    ////////////////
    //Only need to use the ghosts for E (creo) 
    
    for (auto ix = layout->ghost_start(Quantity::Ex, Direction::X); ix <= layout->ghost_end(Quantity::Ex, Direction::X); ++ix)
    {
        E.x(ix) = 0.0; //Ex is dual, though I am setting it at zero  
        Ex_vals.push_back(E.x(ix));
    }
    
    for (auto ix = layout->ghost_start(Quantity::Ey, Direction::X); ix <= layout->ghost_end(Quantity::Ey, Direction::X); ++ix)
    {
        E.y(ix) = std::sin(layout->coordinate(Direction::X, Quantity::Ey, ix)); //Ey is dual  
        Ey_vals.push_back(E.y(ix));
    }
    
    for (auto ix = layout->ghost_start(Quantity::Ez, Direction::X); ix <= layout->ghost_end(Quantity::Ez, Direction::X); ++ix)
    {
        E.z(ix) = std::cos(layout->coordinate(Direction::X, Quantity::Ez, ix)); //Ez is dual  
        Ez_vals.push_back(E.z(ix));
    }


    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        double x_primal = layout->coordinate(Direction::X, Quantity::Bx, ix); // To plot primal things  
        x_primal_vals.push_back(x_primal);

        B.x(ix) = 0.0; //Bx is dual, though I am setting it at zero
        Bx_vals.push_back(B.x(ix));
    }
    
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        double x_dual = layout->coordinate(Direction::X, Quantity::By, ix); // To plot dual things
        x_dual_vals.push_back(x_dual);
        
        B.y(ix) = 0.0; //By is dual, though I am setting it at zero
        By_vals.push_back(B.y(ix));
        
        B.z(ix) = 0.0; //Bz is dual, though I am setting it at zero
        Bz_vals.push_back(B.z(ix));       
    }

    Faraday<dimension> faraday(layout);
    faraday(B, E, Bnew, dt);

    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
    {
        Bxnew_vals.push_back(Bnew.x(ix)); //Primal (though it will be zero)
    }    
    
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X); ++ix)
    {
        Bynew_vals.push_back(Bnew.y(ix)); //Dual  
        Bznew_vals.push_back(Bnew.z(ix)); //Dual 
    }    
  
    
    
    {
        std::string filename = "faraday_test.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/xprimal", x_primal_vals);
        file.createDataSet("/xdual", x_dual_vals);
        file.createDataSet("/Bx", Bx_vals);
        file.createDataSet("/By", By_vals);
        file.createDataSet("/Bz", Bz_vals);
        file.createDataSet("/Ex", Ex_vals);
        file.createDataSet("/Ey", Ey_vals);
        file.createDataSet("/Ez", Ez_vals);
        file.createDataSet("/Bxnew", Bxnew_vals);
        file.createDataSet("/Bynew", Bynew_vals);
        file.createDataSet("/Bznew", Bznew_vals);
        file.createDataSet("/dt", dt);
    }
}


int main()
{
    faraday_test();
}
