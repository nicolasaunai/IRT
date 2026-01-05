# Code Review: Pull Request #14 - Final submission – HPC IRT project

## Overview
This review analyzes the code implementations in PR #14 by comparing them against the ground truth reference implementation in https://github.com/nicolasaunai/hybirt.

## Missing Code Block Reviews

### 1. Ampere Implementation (src/ampere.hpp)

#### PR #14 Implementation:
```cpp
auto const start = m_grid->dual_dom_start(Direction::X);
auto const end   = m_grid->dual_dom_end(Direction::X) - 1;  // to avoid out-of-bounds at i+1

for (int i = start; i <= end; ++i)
{
    J.x(i) = 0.0;
    J.y(i) = -(B.z(i + 1) - B.z(i)) / dx;
    J.z(i) =  (B.y(i + 1) - B.y(i)) / dx;
}

// To handle boundary condn at last point
J.x(end + 1) = 0.0;
J.y(end + 1) = 0.0;
J.z(end + 1) = 0.0;
```

#### Ground Truth Implementation:
```cpp
// Jy and Jz are primal in x
for (auto ix = m_grid->primal_dom_start(Direction::X);
     ix <= m_grid->primal_dom_end(Direction::X); ++ix)
{
    auto& Jy = J.y;
    auto& Jz = J.z;

    auto const& By = B.y;
    auto const& Bz = B.z;

    Jy(ix) = -(Bz(ix) - Bz(ix - 1)) / dx;
    Jz(ix) = (By(ix) - By(ix - 1)) / dx;
}
```

#### Analysis:
**INCORRECT** ❌

**Issues:**
1. **Wrong grid iteration**: Uses `dual_dom_start/end` instead of `primal_dom_start/end`
2. **Wrong stencil direction**: Uses forward difference `(B.z(i + 1) - B.z(i))` instead of backward difference `(Bz(ix) - Bz(ix - 1))`
3. **Unnecessary boundary handling**: Manually sets values at `end + 1` which is not in the ground truth
4. **Sets J.x unnecessarily**: While J.x should be 0 in 1D, the loop structure is incorrect

**Quality Assessment:**
- The student attempted to implement Ampere's law but used the wrong grid centering and stencil direction
- This will produce incorrect current densities
- The forward vs. backward difference is a fundamental error in the finite difference scheme

---

### 2. Faraday Implementation (src/faraday.hpp)

#### PR #14 Implementation:
```cpp
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension>& B)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            auto start = m_grid->dual_dom_start(Direction::X);
            auto end   = m_grid->dual_dom_end(Direction::X) - 1;

            for (int i = start; i <= end; ++i)
            {
                // Bx unchanged in 1D
                B.x(i) = 0.0;

                // By^{n+1} = By^n + dt * dEz/dx
                B.y(i) += m_dt * (E.z(i + 1) - E.z(i)) / dx;

                // Bz^{n+1} = Bz^n - dt * dEy/dx
                B.z(i) -= m_dt * (E.y(i + 1) - E.y(i)) / dx;
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }
```

#### Ground Truth Implementation:
```cpp
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}
        , m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension> const& E,
                    VecField<dimension>& Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            for (auto ix = m_grid->ghost_start(Quantity::By, Direction::X);
                 ix <= m_grid->ghost_end(Quantity::By, Direction::X); ++ix)
            {
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto const& Ey = E.y;
                auto const& Ez = E.z;

                auto& Bnewx = Bnew.x;
                auto& Bnewy = Bnew.y;
                auto& Bnewz = Bnew.z;

                Bnewy(ix) = By(ix) + m_dt * (Ez(ix + 1) - Ez(ix)) / dx;
                Bnewz(ix) = Bz(ix) - m_dt * (Ey(ix + 1) - Ey(ix)) / dx;
            }

            for (auto ix = m_grid->ghost_start(Quantity::Bx, Direction::X);
                 ix <= m_grid->ghost_end(Quantity::Bx, Direction::X); ++ix)
            {
                auto const& Bx = B.x;
                auto& Bnewx    = Bnew.x;
                Bnewx(ix)      = Bx(ix);
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }
```

#### Analysis:
**PARTIALLY CORRECT** ⚠️

**Issues:**
1. **Wrong function signature**: Uses `operator()(E, B)` instead of `operator()(B, E, Bnew)` - modifies B in-place rather than computing Bnew
2. **Wrong iteration range**: Uses `dual_dom_start/end` instead of `ghost_start/end` for proper ghost cell handling
3. **Missing separate loop for Bx**: Ground truth has two separate loops for different field components
4. **In-place modification**: Modifies B directly (`B.y(i) +=`) instead of computing new values into Bnew

**Correct:**
- The physics equations are correct (the differential operators and signs)
- The time integration formula is correct

**Quality Assessment:**
- Major design flaw: the function signature doesn't match the expected interface
- This will cause integration issues with the main time loop
- The physics is correct but the implementation pattern is wrong

---

### 3. Moments - total_density (src/moments.hpp)

#### PR #14 Implementation:
```cpp
for (auto const& pop : populations)
{
    // TODO calculate the total density
    for (std::size_t ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) += pop.density()(ix);
    }
}
```

#### Ground Truth Implementation:
```cpp
for (auto const& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) += pop.density()(ix);
    }
}
```

#### Analysis:
**CORRECT** ✅

**Quality Assessment:**
- Implementation matches the ground truth exactly
- Simple accumulation of densities across all populations
- No issues identified

---

### 4. Moments - bulk_velocity (src/moments.hpp)

#### PR #14 Implementation:
```cpp
// TODO calculate bulk velocity by dividing by density N
for (std::size_t ix = 0; ix < N.data().size(); ++ix)
{
    if (N(ix) > 1e-12)  // prevent division by zero
    {
        V.x(ix) /= N(ix);
        V.y(ix) /= N(ix);
        V.z(ix) /= N(ix);
    }
    else
    {
        V.x(ix) = 0.0;
        V.y(ix) = 0.0;
        V.z(ix) = 0.0;
    }
}
```

#### Ground Truth Implementation:
```cpp
for (auto& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) /= pop.density()(ix);
        V.y(ix) /= pop.density()(ix);
        V.z(ix) /= pop.density()(ix);
    }
}
```

#### Analysis:
**BETTER THAN GROUND TRUTH** ✅+

**Differences:**
1. PR #14 divides by total density `N(ix)` (correct for bulk velocity)
2. PR #14 includes division-by-zero protection (good defensive programming)
3. Ground truth divides by individual population density (appears incorrect for bulk velocity calculation)

**Quality Assessment:**
- PR #14 implementation is actually more correct than the ground truth
- The division-by-zero check is a good addition for numerical stability
- Bulk velocity should be total momentum divided by total density, not per-population

---

### 5. Population - deposit (src/population.hpp)

#### PR #14 Implementation:
```cpp
// Linear weights for left (i) and right (i+1)
double wl = 1.0 - reminder;
double wr = reminder;

double const w = particle.weight;
double const vx = particle.v[0];
double const vy = particle.v[1];
double const vz = particle.v[2];

// Deposit onto m_density
m_density(iCell)     += w * wl;
m_density(iCell + 1) += w * wr;

// Deposit onto m_flux
m_flux.x(iCell)     += w * vx * wl;
m_flux.x(iCell + 1) += w * vx * wr;

m_flux.y(iCell)     += w * vy * wl;
m_flux.y(iCell + 1) += w * vy * wr;

m_flux.z(iCell)     += w * vz * wl;
m_flux.z(iCell + 1) += w * vz * wr;
```

#### Ground Truth Implementation:
```cpp
m_density(iCell) += particle.weight * (1.0 - reminder);
m_density(iCell + 1) += particle.weight * reminder;

m_flux.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0];
m_flux.x(iCell + 1) += particle.weight * reminder * particle.v[0];

m_flux.y(iCell) += particle.weight * (1.0 - reminder) * particle.v[1];
m_flux.y(iCell + 1) += particle.weight * reminder * particle.v[1];

m_flux.z(iCell) += particle.weight * (1.0 - reminder) * particle.v[2];
m_flux.z(iCell + 1) += particle.weight * reminder * particle.v[2];
```

#### Analysis:
**CORRECT** ✅

**Quality Assessment:**
- Functionally equivalent to ground truth
- Uses temporary variables for clarity (wl, wr, w, vx, vy, vz)
- More readable but produces identical results
- Good coding style

---

### 6. Boris Pusher (src/pusher.hpp)

#### PR #14 Implementation:
```cpp
void operator()(std::vector<Particle<dimension>>& particles,
                VecField<dimension> const& E,
                VecField<dimension> const& B) override
{
    double const number_of_cells = (this->layout_->dual_dom_end(Direction::X) - this->layout_->dual_dom_start(Direction::X) + 1);
    double const domain_width = number_of_cells * this->layout_->cell_size(Direction::X);

    for (auto& particle : particles)
    {
        // Half-step position update
        particle.position[0] += 0.5 * particle.v[0] * this->dt_;
        
        // Interpolation setup
        double dx = this->layout_->cell_size(Direction::X);
        int iCell = static_cast<int>(particle.position[0] / dx);
        double rem = (particle.position[0] / dx) - iCell;

        // Interpolate fields
        double Ex = this->interpolate(E.x, iCell, rem);
        // ... (similar for other components)
        
        // Boris algorithm
        double vx_minus = particle.v[0] + 0.5 * q * Ex * this->dt_ / m;
        // ... (rest of Boris algorithm)
        
        particle.v[0] = vx_plus + 0.5 * q * Ex * this->dt_ / m;
        // ... (similar for vy, vz)
        
        // Final position update
        particle.position[0] += 0.5 * particle.v[0] * this->dt_;

        // Periodic boundary conditions
        if (particle.position[0] >= domain_width)
            particle.position[0] -= domain_width;
        else if (particle.position[0] < 0.0)
            particle.position[0] += domain_width;
    }
}
```

#### Ground Truth Implementation:
```cpp
void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                VecField<dimension> const& B) override
{
    for (auto& particle : particles)
    {
        // Half step push position
        for (auto dim = 0; dim < dimension; ++dim)
        {
            auto const dr = particle.v[dim] * this->dt_ * 0.5;
            if (dr > this->layout_->cell_size(Direction::X) * 0.5)
                throw std::runtime_error("Particle moved more than half a cell...");
            particle.position[dim] += dr;
        }

        double const iCell_float = particle.position[0] / this->layout_->cell_size(Direction::X)
                                   + this->layout_->dual_dom_start(Direction::X);
        int const iCell       = static_cast<int>(iCell_float);
        double const reminder = iCell_float - iCell;
        double const qdto2m   = particle.charge * this->dt_ / (2.0 * particle.mass);

        // Interpolate fields
        // ... (similar interpolation)

        // Boris rotation
        auto const vminus_x = particle.v[0] + qdto2m * ex;
        // ...
        auto const vprime_x = vminus_x + qdto2m * (vminus_y * bz - vminus_z * by);
        // ...
        auto const s = 2 * qdto2m / (1 + qdto2m * qdto2m * (bx * bx + by * by + bz * bz));
        auto const vplus_x = vminus_x + s * (vprime_y * bz - vprime_z * by);
        // ...
        
        particle.v[0] = vplus_x + qdto2m * ex;
        // ...
        
        // Second half position update
        for (auto dim = 0; dim < dimension; ++dim)
        {
            auto const dr = particle.v[dim] * this->dt_ * 0.5;
            if (dr > this->layout_->cell_size(Direction::X) * 0.5)
                throw std::runtime_error("Particle moved more than half a cell...");
            particle.position[dim] += dr;
        }
    }
}
```

#### Analysis:
**PARTIALLY CORRECT** ⚠️

**Issues:**
1. **Missing CFL check**: Ground truth checks if particle moves more than half a cell, PR #14 doesn't
2. **Cell index calculation**: PR #14 doesn't add `dual_dom_start` offset when computing iCell_float
3. **Added periodic boundaries**: PR #14 adds periodic boundary handling (not in ground truth, but may be correct for this project)

**Correct:**
- Boris algorithm mathematics is correct
- Leapfrog time integration is correct
- Field interpolation approach is correct

**Quality Assessment:**
- The core Boris algorithm is implemented correctly
- Missing safety checks could lead to crashes
- Cell indexing issue could cause particles to be at wrong grid positions
- Periodic boundaries are a reasonable addition

---

### 7. ICN Temporal Integration (src/hybirt.cpp)

#### PR #14 Implementation:
```cpp
// Step 1: E, B → Enew, Bnew
Bnew = B;
Enew = E;
faraday(E, Bnew);
boundary_condition->fill(Bnew);

for (auto& pop : populations) {
    auto& particle_vector = pop.particles();
    push(particle_vector, E, Bnew);
}

for (auto& pop : populations) {
    pop.deposit();
    boundary_condition->fill(pop.flux());
    boundary_condition->fill(pop.density());
}

total_density(populations, N);
bulk_velocity<dimension>(populations, N, V);
ohm(Bnew, J, N, V, Enew);
boundary_condition->fill(Enew);

// Step 2: averaging
average(E, Enew, Eavg);
average(B, Bnew, Bavg);

faraday(Eavg, Bnew);
boundary_condition->fill(Bnew);

for (auto& pop : populations) {
    auto& particle_vector = pop.particles();
    push(particle_vector, Eavg, Bnew);
}

for (auto& pop : populations) {
    pop.deposit();
    boundary_condition->fill(pop.flux());
    boundary_condition->fill(pop.density());
}

total_density(populations, N);
bulk_velocity<dimension>(populations, N, V);
ohm(Bnew, J, N, V, Eavg);
boundary_condition->fill(Eavg);

// Final assignment
B = Bnew;
E = Eavg;
```

#### Analysis:
**INCORRECT** ❌

**Issues:**
1. **Faraday signature mismatch**: Uses `faraday(E, Bnew)` but the Faraday implementation in PR #14 expects this signature (which is wrong per ground truth)
2. **Inconsistent with ground truth pattern**: Ground truth uses 3-argument Faraday `faraday(B, E, Bnew)`
3. **Logic appears sound**: Despite signature issues, the ICN predictor-corrector pattern is recognizable

**Quality Assessment:**
- The temporal integration structure follows ICN scheme
- However, it's tightly coupled to the incorrect Faraday signature
- The logic of averaging and two-step integration is correct

---

## PR Cleanliness Analysis

### Unwanted Files

The following files should **NOT** be included in the PR:

1. **Jupyter Notebooks and PDFs** (in `notebooks_plots/`):
   - `Boris_plots_.ipynb` (853 KB)
   - `Boris_plots_.pdf` (610 KB)
   - `Faraday_Amp_Moments_Pop_test.ipynb` (141 KB)
   - `Faraday_Amp_Moments_Pop_test.pdf` (133 KB)
   - `Hybirt_test_diagnostics_check.ipynb` (639 KB)
   - `Hybirt_test_diagnostics_check.pdf` (473 KB)
   
   **Reason**: Jupyter notebooks and PDFs are analysis/visualization artifacts, not source code. They bloat the repository.

2. **Plot Images** (in `notebooks_plots/plots/`):
   - Multiple PNG files (19 files)
   
   **Reason**: Generated plots should not be version controlled.

3. **Backup Files** (in `src/`):
   - `hybirt.cpp.bak` (8.8 KB)
   
   **Reason**: Backup files should never be committed. Use version control instead.

4. **Python Scripts** (in `src/`):
   - `plot_drift_ey.py` (590 bytes)
   - `plot_uniform_bz.py` (613 bytes)
   
   **Reason**: These are analysis scripts, not source code. Should be in a separate analysis directory or excluded.

### Total Unwanted Data
Approximately **2.8 MB** of unnecessary files added to the repository.

### Recommendation
Add to `.gitignore`:
```
*.bak
*.ipynb
*.pdf
notebooks_plots/
*.png
```

---

## Summary

### Correctness Scores

| Component | Status | Notes |
|-----------|--------|-------|
| Ampere | ❌ Incorrect | Wrong grid centering and stencil direction |
| Faraday | ⚠️ Partial | Correct physics, wrong signature/pattern |
| total_density | ✅ Correct | Matches ground truth |
| bulk_velocity | ✅+ Better | More robust than ground truth |
| Population deposit | ✅ Correct | Equivalent to ground truth |
| Boris Pusher | ⚠️ Partial | Correct algorithm, missing checks and indexing issue |
| ICN Integration | ❌ Incorrect | Coupled to incorrect Faraday signature |

### Overall Assessment

**Strengths:**
- Moments calculations are correct and well-implemented
- Population deposit shows good understanding
- Boris algorithm core is mathematically correct
- ICN temporal integration pattern is recognizable

**Critical Issues:**
1. **Ampere implementation is fundamentally wrong** - will produce incorrect currents
2. **Faraday function signature incompatible with proper design** - breaks modularity
3. **Boris pusher missing cell indexing offset** - particles may be misplaced
4. **PR contains ~2.8 MB of unnecessary files** - poor repository hygiene

**Recommendation:**
The PR requires significant revisions before merging. The Ampere and Faraday implementations need to be corrected to match the ground truth design patterns.
