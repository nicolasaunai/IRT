#include "faraday.hpp"
#include "gridlayout.hpp"
#include "vecfield.hpp"
#include "boundary_condition.hpp"
#include "diagnostics.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

double Ex_fun(double x) { return std::sin(2*M_PI*x); }
double Ey_fun(double x) { return std::cos(2*M_PI*x); }
double Ez_fun(double x) { return std::sin(4*M_PI*x); }

void initialize_E(VecField<1>& E, GridLayout<1> const& layout) {
    for (std::size_t ix = layout.primal_dom_start(Direction::X);
         ix <= layout.primal_dom_end(Direction::X); ++ix)
    {
        double x = layout.coordinate(Direction::X, Quantity::Ex, ix);
        E.x(ix) = Ex_fun(x);
    }
    for (std::size_t ix = layout.dual_dom_start(Direction::X);
         ix <= layout.dual_dom_end(Direction::X); ++ix)
    {
        double x = layout.coordinate(Direction::X, Quantity::Ey, ix);
        E.y(ix) = Ey_fun(x);
        E.z(ix) = Ez_fun(x);
    }
}

int main() {
    constexpr std::size_t dimension = 1;
    std::array<std::size_t, dimension> grid_size = {50};
    std::array<double, dimension> cell_size = {1.0 / 50};
    constexpr std::size_t nbr_ghosts = 1;
    double dt = 0.01;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> B(layout, {Quantity::Bx, Quantity::By, Quantity::Bz});
    VecField<dimension> Bnew(layout, {Quantity::Bx, Quantity::By, Quantity::Bz});
    VecField<dimension> E(layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez});

    initialize_E(E, *layout);

    for (std::size_t ix = layout->primal_dom_start(Direction::X);
         ix <= layout->primal_dom_end(Direction::X); ++ix) B.x(ix) = 0.0;
    for (std::size_t ix = layout->dual_dom_start(Direction::X);
         ix <= layout->dual_dom_end(Direction::X); ++ix) {
        B.y(ix) = 0.0;
        B.z(ix) = 0.0;
    }

    Faraday<dimension> faraday(layout, dt);

    faraday(E, B, Bnew);

    HighFive::File file("faraday_test.h5", HighFive::File::Truncate);

    std::vector<double> Bx_vec, By_vec, Bz_vec;
    for (std::size_t ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X); ++ix)
        Bx_vec.push_back(Bnew.x(ix));
    for (std::size_t ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X); ++ix) {
        By_vec.push_back(Bnew.y(ix));
        Bz_vec.push_back(Bnew.z(ix));
    }

    file.createDataSet("/Bx", Bx_vec);
    file.createDataSet("/By", By_vec);
    file.createDataSet("/Bz", Bz_vec);

    std::cout << "Faraday test complete, B saved to faraday_test.h5\n";

    return 0;
}
