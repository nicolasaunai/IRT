#include "population.hpp"
#include "gridlayout.hpp"
#include "field.hpp"
#include "vecfield.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>



template<class P>
double get_x(P const& p)
{
    if constexpr (requires { p.position[0]; }) return p.position[0];
    else if constexpr (requires { p.r[0]; })  return p.r[0];
    else {
        static_assert(sizeof(P) == 0, "Particle type has no position[0] or r[0]");
        return 0.0;
    }
}

template<class P>
void set_x(P& p, double x)
{
    if constexpr (requires { p.position[0] = x; }) p.position[0] = x;
    else if constexpr (requires { p.r[0] = x; })  p.r[0] = x;
    else static_assert(sizeof(P) == 0, "Particle type has no writable position[0] or r[0]");
}

template<class P>
double get_vx(P const& p)
{
    if constexpr (requires { p.v[0]; }) return p.v[0];
    else if constexpr (requires { p.velocity[0]; }) return p.velocity[0];
    else {
        static_assert(sizeof(P) == 0, "Particle type has no v[0] or velocity[0]");
        return 0.0;
    }
}

template<class P>
void set_vx(P& p, double vx)
{
    if constexpr (requires { p.v[0] = vx; }) p.v[0] = vx;
    else if constexpr (requires { p.velocity[0] = vx; }) p.velocity[0] = vx;
    else static_assert(sizeof(P) == 0, "Particle type has no writable v[0] or velocity[0]");
}

template<class P>
void set_vy(P& p, double vy)
{
    if constexpr (requires { p.v[1] = vy; }) p.v[1] = vy;
    else if constexpr (requires { p.velocity[1] = vy; }) p.velocity[1] = vy;
    else static_assert(sizeof(P) == 0, "Particle type has no writable v[1] or velocity[1]");
}

template<class P>
void set_vz(P& p, double vz)
{
    if constexpr (requires { p.v[2] = vz; }) p.v[2] = vz;
    else if constexpr (requires { p.velocity[2] = vz; }) p.velocity[2] = vz;
    else static_assert(sizeof(P) == 0, "Particle type has no writable v[2] or velocity[2]");
}

/* ---------- Utility: collect nonzero entries ---------- */
static std::vector<std::pair<std::size_t,double>> nonzeros_1d(auto const& field, double eps=1e-14)
{
    std::vector<std::pair<std::size_t,double>> nz;
    auto const n = field.data().size();
    for (std::size_t i = 0; i < n; ++i)
    {
        double val = field(i);
        if (std::abs(val) > eps) nz.emplace_back(i, val);
    }
    return nz;
}

int main()
{
    std::cout << "Testing Population::deposit() with a single deterministic particle...\n";

    constexpr std::size_t dim = 1;

    std::array<std::size_t, dim> grid_size = {16};
    std::array<double, dim>      cell_size = {0.2};
    constexpr std::size_t nbr_ghosts = 1;

    auto layout = std::make_shared<GridLayout<dim>>(grid_size, cell_size, nbr_ghosts);

    Population<dim> pop("pop_test", layout);

    /* --- Hand-craft exactly ONE particle --- */
    auto& parts = pop.particles();
    parts.clear();
    parts.resize(1);

    auto& p = parts[0];

    /* Put particle at x = (k + frac)*dx with frac = 0.30 (so weights are 0.70 and 0.30) */
    double const dx   = cell_size[0];
    std::size_t const k = 5;
    double const frac = 0.30;
    double const x    = (static_cast<double>(k) + frac) * dx;

    /* Choose weight and velocity so flux is easy to check */
    double const W  = 2.0;
    double const vx = 0.5;

    set_x(p, x);
    if constexpr (requires { p.weight = W; })  p.weight = W;
    else { std::cerr << "Particle has no .weight member\n"; return 2; }

    if constexpr (requires { p.charge = 1.0; }) p.charge = 1.0;  /* not needed for deposit check but should exist */

    set_vx(p, vx);
    set_vy(p, 0.0);
    set_vz(p, 0.0);

    /* --- Run deposit --- */
    pop.deposit();

    /* --- Check deposited density: exactly 2 nonzero entries, sum == W, values == W*(1-frac), W*frac --- */
    auto nzN = nonzeros_1d(pop.density());
    if (nzN.size() != 2)
    {
        std::cerr << "FAIL: density deposit should touch exactly 2 nodes, got " << nzN.size() << "\n";
        std::cerr << "Nonzeros were:\n";
        for (auto const& [i,v] : nzN) std::cerr << "  i=" << i << "  n=" << v << "\n";
        return 1;
    }

    std::sort(nzN.begin(), nzN.end(), [](auto const& a, auto const& b){ return a.first < b.first; });

    double const wl = 1.0 - frac;
    double const wr = frac;
    double const expected_leftN  = W * wl;
    double const expected_rightN = W * wr;

    double const got_leftN  = nzN[0].second;
    double const got_rightN = nzN[1].second;

    double const sumN = got_leftN + got_rightN;

    auto relerr = [](double a, double b){
        double den = std::max(1.0, std::max(std::abs(a), std::abs(b)));
        return std::abs(a-b)/den;
    };

    if (relerr(sumN, W) > 1e-12)
    {
        std::cerr << "FAIL: density conservation broken. sum(deposit)=" << sumN << " vs W=" << W << "\n";
        return 1;
    }

    if (relerr(got_leftN, expected_leftN) > 1e-12 || relerr(got_rightN, expected_rightN) > 1e-12)
    {
        std::cerr << "FAIL: wrong linear weights in density deposit.\n";
        std::cerr << "Expected (left,right)=(" << expected_leftN << "," << expected_rightN << ")\n";
        std::cerr << "Got      (left,right)=(" << got_leftN      << "," << got_rightN      << ")\n";
        return 1;
    }

    /* --- Check deposited flux in x, if available as pop.flux().x --- */
    auto nzFx = nonzeros_1d(pop.flux().x);
    if (nzFx.size() != 2)
    {
        std::cerr << "FAIL: flux-x deposit should touch exactly 2 nodes, got " << nzFx.size() << "\n";
        std::cerr << "Nonzeros were:\n";
        for (auto const& [i,v] : nzFx) std::cerr << "  i=" << i << "  Fx=" << v << "\n";
        return 1;
    }

    std::sort(nzFx.begin(), nzFx.end(), [](auto const& a, auto const& b){ return a.first < b.first; });

    double const expected_leftFx  = (W * vx) * wl;
    double const expected_rightFx = (W * vx) * wr;

    double const got_leftFx  = nzFx[0].second;
    double const got_rightFx = nzFx[1].second;

    if (relerr(got_leftFx, expected_leftFx) > 1e-12 || relerr(got_rightFx, expected_rightFx) > 1e-12)
    {
        std::cerr << "FAIL: wrong linear weights in flux deposit.\n";
        std::cerr << "Expected (left,right)=(" << expected_leftFx << "," << expected_rightFx << ")\n";
        std::cerr << "Got      (left,right)=(" << got_leftFx      << "," << got_rightFx      << ")\n";
        return 1;
    }

    std::cout << "PASS.\n";
    std::cout << "  particle x=" << get_x(p) << " (frac=" << frac << "), W=" << W << ", vx=" << get_vx(p) << "\n";
    std::cout << "  density nodes: i=" << nzN[0].first << " -> " << got_leftN
              << " , i=" << nzN[1].first << " -> " << got_rightN << "\n";
    std::cout << "  flux-x nodes:  i=" << nzFx[0].first << " -> " << got_leftFx
              << " , i=" << nzFx[1].first << " -> " << got_rightFx << "\n";
    return 0;
}

