// test_faraday.cpp
#include "faraday.hpp"
#include "gridlayout.hpp"
#include "vecfield.hpp"
#include <cmath>
#include <iostream>

#include "highfive/highfive.hpp"

int main() {
    constexpr std::size_t dim = 1;
    std::array<std::size_t, dim> grid_size = {100};
    std::array<double, dim> cell_size = {0.1};
    auto layout = std::make_shared<GridLayout<dim>>(grid_size, cell_size, 1);

    VecField<dim> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dim> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dim> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    const double k = 2.0;
    const double dt = 0.01;
    for (auto ix = layout->dual_dom_start(Direction::X);
         ix <= layout->dual_dom_end(Direction::X); ++ix) {
        double x = ix * cell_size[0];
        E.y(ix) = std::sin(k * x);
        E.z(ix) = std::cos(k * x);
        B.y(ix) = 0.0;
        B.z(ix) = 0.0;
    }

    Faraday<dim> faraday{layout, dt};
    faraday(E, B, Bnew);

    std::string filename = "test_faraday.h5";
    HighFive::File file(filename, HighFive::File::Truncate);
    file.createDataSet("/dBy", Bnew.y.data());
    file.createDataSet("/dBz", Bnew.z.data());

    
}
