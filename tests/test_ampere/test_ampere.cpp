// test_ampere.cpp
#include "ampere.hpp"
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

    VecField<dim> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dim> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};

    const double k = 2.0;
    for (auto ix = layout->dual_dom_start(Direction::X);
         ix <= layout->dual_dom_end(Direction::X); ++ix) {
        double x = ix * cell_size[0];
        B.y(ix) = std::sin(k * x);
        B.z(ix) = std::cos(k * x);
    }

    Ampere<dim> amp{layout};
    amp(B, J);  

    std::string filename = "test_ampere.h5";
    HighFive::File file(filename, HighFive::File::Truncate);
    file.createDataSet("/Jy", J.y.data());
    file.createDataSet("/Jz", J.z.data());
}
