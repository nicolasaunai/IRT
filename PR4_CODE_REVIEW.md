# PR #4 Code Review Report

This report analyzes the code implemented in PR #4 by comparing it against the ground truth implementation from https://github.com/nicolasaunai/hybirt.

## Overview of TODO Blocks in Master Branch

The following TODO statements were found in the master branch:
1. `src/faraday.hpp:15` - Implement the Faraday class
2. `src/ampere.hpp:26` - Implement Ampere calculation
3. `src/population.hpp:111` - Implement linear weighting deposit for density and flux
4. `src/pusher.hpp:49` - Implement the Boris pusher
5. `src/moments.hpp:20` - Calculate total density
6. `src/moments.hpp:43` - Calculate bulk velocity by dividing by density N
7. `src/hybirt.cpp:128` - Uncomment Faraday when implemented
8. `src/hybirt.cpp:156` - Implement ICN temporal integration

## Detailed Code Review by Missing Block

### 1. Ampere Implementation (`src/ampere.hpp`)

**PR #4 Implementation:**
```cpp
for (auto ix = m_grid->dual_dom_start(Direction::X);
     ix <= m_grid->dual_dom_end(Direction::X); ++ix)
{
    auto& Jx = J.x;
    Jx(ix) = 0.0;
}
// Jy, Jz are primal
for (auto ix = m_grid->primal_dom_start(Direction::X);
     ix <= m_grid->primal_dom_end(Direction::X); ++ix)
{
    auto const& Bx = B.x;
    auto const& By = B.y;
    auto const& Bz = B.z;
    
    auto& Jy = J.y;
    auto& Jz = J.z;
    
    Jy(ix) = - (Bz(ix)- Bz(ix - 1))/(dx);
    Jz(ix) = (By(ix)- By(ix - 1))/(dx);
}
```

**Ground Truth Implementation:**
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

**Analysis:**
- ✅ **Correctness**: The formulas for Jy and Jz are correct and match the ground truth
- ⚠️ **Quality Issues**:
  1. **Unnecessary code**: PR #4 includes a loop setting Jx to 0.0, which is not present in ground truth and appears unnecessary
  2. **Unused variables**: The variable `Bx` is declared but never used in PR #4
  3. **Variable order**: Ground truth declares references in a cleaner order
  4. **Spacing**: Missing space after minus sign in PR #4 (`Bz(ix)- Bz`)
- **Overall**: Functionally correct but contains unnecessary code and minor style inconsistencies

### 2. Faraday Implementation (`src/faraday.hpp`)

**PR #4 Implementation:**
```cpp
for (auto ix = m_grid->primal_dom_start(Direction::X);
     ix <= m_grid->primal_dom_end(Direction::X); ++ix)
{
    auto& Bx = Bnew.x;
    Bx(ix) = B.x(ix);
}
// By, Bz are dual
for (auto ix = m_grid->dual_dom_start(Direction::X);
     ix <= m_grid->dual_dom_end(Direction::X); ++ix)
{
    auto const& Ex = E.x;
    auto const& Ey = E.y;
    auto const& Ez = E.z;
    
    auto& By = Bnew.y;
    auto& Bz = Bnew.z;
    
    By(ix) = B.y(ix) + (Ez(ix + 1) - Ez(ix))*m_dt/dx; 
    Bz(ix) = B.z(ix) - (Ey(ix + 1) - Ey(ix))*m_dt/dx;
}
```

**Ground Truth Implementation:**
```cpp
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
```

**Analysis:**
- ✅ **Correctness**: The Faraday update formulas are mathematically correct
- ⚠️ **Quality Issues**:
  1. **Loop bounds**: PR #4 uses `dual_dom_start/end` and `primal_dom_start/end` instead of the more appropriate `ghost_start/end` used in ground truth
  2. **Loop order**: PR #4 updates Bx first, then By/Bz; ground truth does the opposite (By/Bz first, then Bx separately)
  3. **Unused variable**: `Ex` is declared but never used in PR #4
  4. **Operator spacing**: Missing spaces around operators in PR #4 (`*m_dt/dx` vs `m_dt * ... / dx`)
  5. **Function signature**: Ground truth has different parameter order: `(B, E, Bnew)` vs PR #4's `(E, B, Bnew)`
- **Overall**: Functionally correct but uses different loop bounds which may affect ghost cell handling

### 3. Population Deposit Implementation (`src/population.hpp`)

**PR #4 Implementation:**
```cpp
double const dx    = m_grid->cell_size(Direction::X);
double const x     = std::fmod(std::fmod(particle.position[0], dx * m_grid->nbr_cells(Direction::X)) + dx * m_grid->nbr_cells(Direction::X),
                               dx * m_grid->nbr_cells(Direction::X));
double       s     = x / dx;                              // in [0, Nx)
int          i0    = static_cast<int>(std::floor(s));     // left node (primal)
double       w1    = s - i0;                              // weight to the right
double       w0    = 1.0 - w1;

int left  = i0 + m_grid->dual_dom_start(Direction::X);
int right = left + 1;

// wrap the right index onto the next cell (periodic)
if (right > m_grid->dual_dom_end(Direction::X)) {
    right = m_grid->dual_dom_start(Direction::X);
}

m_density(left)  += particle.weight * w0;
m_density(right) += particle.weight * w1;

m_flux.x(left)   += particle.weight * particle.v[0] * w0;
m_flux.x(right)  += particle.weight * particle.v[0] * w1;
m_flux.y(left)   += particle.weight * particle.v[1] * w0;
m_flux.y(right)  += particle.weight * particle.v[1] * w1;
m_flux.z(left)   += particle.weight * particle.v[2] * w0;
m_flux.z(right)  += particle.weight * particle.v[2] * w1;
```

**Ground Truth Implementation:**
```cpp
double const iCell_float = particle.position[0] / m_grid->cell_size(Direction::X);
int const iCell_         = static_cast<int>(iCell_float);
double const reminder    = iCell_float - iCell_;
auto const iCell         = iCell_ + m_grid->dual_dom_start(Direction::X);


m_density(iCell) += particle.weight * (1.0 - reminder);
m_density(iCell + 1) += particle.weight * reminder;

m_flux.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0];
m_flux.x(iCell + 1) += particle.weight * reminder * particle.v[0];

m_flux.y(iCell) += particle.weight * (1.0 - reminder) * particle.v[1];
m_flux.y(iCell + 1) += particle.weight * reminder * particle.v[1];

m_flux.z(iCell) += particle.weight * (1.0 - reminder) * particle.v[2];
m_flux.z(iCell + 1) += particle.weight * reminder * particle.v[2];
```

**Analysis:**
- ✅ **Correctness**: Both implementations use linear weighting correctly
- ⚠️ **Quality Issues**:
  1. **Complexity**: PR #4 uses nested `std::fmod` for periodic boundary handling which is overly complex
  2. **Explicit wrapping**: PR #4 explicitly wraps the right index, while ground truth relies on boundary conditions
  3. **Cleaner approach**: Ground truth is simpler and more readable
  4. **Variable naming**: Ground truth uses `reminder` consistently; PR #4 uses `w0`/`w1` and separate `left`/`right` variables
- **Overall**: Functionally correct but unnecessarily complex; ground truth is cleaner

### 4. Boris Pusher Implementation (`src/pusher.hpp`)

**PR #4 Implementation:**
```cpp
double x;
double x_half;

double tx, ty, tz;
double sx, sy, sz;

double vx, vx_minus, vx_prime, vx_plus;
double vy, vy_minus, vy_prime, vy_plus;
double vz, vz_minus, vz_prime, vz_plus;

x = particle.position[0];
vx = particle.v[0];
vy = particle.v[1];
vz = particle.v[2];

x_half = x + vx*this->dt_/2;

double const iCell_float = x_half / this->layout_->cell_size(Direction::X)
                       + this->layout_->dual_dom_start(Direction::X);
int const iCell = static_cast<int>(iCell_float);
double reminder = iCell_float - iCell;

// ... field interpolation and Boris algorithm implementation
```

**Ground Truth Implementation:**
```cpp
// half step push position from t=n to t=n+1/2
for (auto dim = 0; dim < dimension; ++dim)
{
    auto const dr = particle.v[dim] * this->dt_ * 0.5;
    if (dr > this->layout_->cell_size(Direction::X) * 0.5)
    {
        throw std::runtime_error(
            "Particle moved more than half a cell size in one step in 1nd update");
    }
    particle.position[dim] += dr;
}

double const iCell_float = particle.position[0] / this->layout_->cell_size(Direction::X)
                           + this->layout_->dual_dom_start(Direction::X);
int const iCell       = static_cast<int>(iCell_float);
double const reminder = iCell_float - iCell;
double const qdto2m   = particle.charge * this->dt_ / (2.0 * particle.mass);

// ... Boris algorithm with cleaner variable naming
```

**Analysis:**
- ✅ **Correctness**: The Boris algorithm is correctly implemented in both versions
- ⚠️ **Quality Issues**:
  1. **Variable declarations**: PR #4 pre-declares many variables at the top, while ground truth uses them as needed
  2. **Safety checks**: Ground truth includes bounds checking to ensure particles don't move more than half a cell
  3. **Position update**: PR #4 uses intermediate `x_half` variable; ground truth directly updates `particle.position[dim]`
  4. **Code style**: Ground truth is more modern C++ with `const` and `auto`
  5. **Formula clarity**: Ground truth uses `qdto2m` constant for repeated charge*dt/(2*mass) calculations
  6. **Missing second half**: PR #4's implementation appears complete, but ground truth explicitly documents the two-step position update
- **Overall**: Functionally correct but ground truth has better safety checks and cleaner code style

### 5. Total Density Implementation (`src/moments.hpp`)

**PR #4 Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    N(ix) = 0;
}
for (auto const& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) += pop.density()(ix);
    }
}
```

**Ground Truth Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    N(ix) = 0;
}
for (auto const& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) += pop.density()(ix);
    }
}
```

**Analysis:**
- ✅ **Correctness**: Implementation is identical to ground truth
- ✅ **Quality**: Perfect match
- **Overall**: Excellent, no issues

### 6. Bulk Velocity Implementation (`src/moments.hpp`)

**PR #4 Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    V.x(ix) = 0;
    V.y(ix) = 0;
    V.z(ix) = 0;
}
for (auto& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) += pop.flux().x(ix);
        V.y(ix) += pop.flux().y(ix);
        V.z(ix) += pop.flux().z(ix);
    }
}
// TODO calculate bulk velocity by dividing by density N
constexpr double N_floor = 1e-12;
for (std::size_t ix = 0; ix < N.data().size(); ++ix) {
    double const n = std::max(N(ix), N_floor);
    V.x(ix) /= n; 
    V.y(ix) /= n; 
    V.z(ix) /= n;
}
```

**Ground Truth Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    V.x(ix) = 0;
    V.y(ix) = 0;
    V.z(ix) = 0;
}
for (auto& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) += pop.flux().x(ix);
        V.y(ix) += pop.flux().y(ix);
        V.z(ix) += pop.flux().z(ix);
    }
}
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

**Analysis:**
- ⚠️ **Correctness**: **DIFFERENT IMPLEMENTATION**
  - PR #4 divides by total density N with a floor value to avoid division by zero
  - Ground truth divides by each population's density individually
  - **PR #4's approach is arguably better** as it avoids division by zero with a safety floor
- ✅ **Quality**: PR #4 has better numerical stability with `N_floor`
- **Overall**: PR #4 implementation is actually superior in terms of numerical safety

### 7. Faraday Uncommented in hybirt.cpp

**PR #4 Implementation:**
```cpp
Faraday<dimension> faraday{layout, dt};  // TODO uncomment when Faraday is implemented
```

**Ground Truth Implementation:**
```cpp
Faraday<dimension> faraday{layout, dt};
```

**Analysis:**
- ✅ **Correctness**: Correctly uncommented
- ⚠️ **Quality**: TODO comment should be removed
- **Overall**: Functionally correct but comment should be cleaned up

### 8. ICN Temporal Integration Implementation (`src/hybirt.cpp`)

**PR #4 Implementation:**
```cpp
// Prediction 1 
faraday(E, B, Bnew);
boundary_condition->fill(Bnew);

ampere(Bnew, J);
boundary_condition->fill(J);

ohm(Bnew, J, N, V, Enew);
boundary_condition->fill(Enew);

average(E, Enew, Eavg);
average(B, Bnew, Bavg);

for (auto& pop : populations)
{
    push(pop.particles(), Eavg, Bavg);
    boundary_condition->particles(pop.particles());
    pop.deposit();
    boundary_condition->fill(pop.flux());
    boundary_condition->fill(pop.density());
}
total_density(populations, N);
bulk_velocity<dimension>(populations, N, V);

// Prediction 2
faraday(Eavg, B, Bnew);
// ... similar structure

// Correction 
faraday(Eavg, B, Bnew);
// ... final update
```

**Ground Truth Implementation:**
```cpp
// predictor 1
faraday(B, E, Bnew);
ampere(Bnew, J);
boundary_condition->fill(J);
ohm(Bnew, J, N, V, Enew);
average(E, Enew, Eavg);
average(B, Bnew, Bavg);
boundary_condition->fill(Eavg);
boundary_condition->fill(Bavg);

// save particles at t=n before updating
auto save{populations};
for (auto& pop : save)
{
    push(pop.particles(), Eavg, Bavg);
    boundary_condition->particles(pop.particles());
    pop.deposit();
    boundary_condition->fill(pop.flux());
    boundary_condition->fill(pop.density());
}
total_density(populations, N);
bulk_velocity<dimension>(populations, N, V);
save.clear(); // clear the save vector to free memory

// predictor 2
faraday(B, Eavg, Bnew);
// ... similar with populations (not save)

// corrector
faraday(B, Eavg, B);
ampere(B, J);
boundary_condition->fill(J);
ohm(B, J, N, V, E);
boundary_condition->fill(E);
```

**Analysis:**
- ⚠️ **Correctness**: **CRITICAL DIFFERENCES**
  1. **Faraday parameter order**: PR #4 uses `faraday(E, B, Bnew)` while ground truth uses `faraday(B, E, Bnew)` - this matches the signature difference noted in section 2
  2. **Particle save/restore**: Ground truth saves particles before predictor 1 and appears to use separate handling; PR #4 doesn't save/restore
  3. **Boundary filling**: PR #4 fills Bnew after Faraday in predictor steps; ground truth fills Eavg and Bavg after averaging
  4. **Corrector step**: Ground truth updates B in-place `faraday(B, Eavg, B)`; PR #4 updates Bnew then copies
- **Overall**: Significant algorithmic differences that may affect correctness - needs testing

## PR Cleanliness Review

### Files That Should Not Be in the Repository

The following files/directories were found that should NOT be in a production repository:

1. **Jupyter Notebook Checkpoints** (Multiple locations):
   - `./.ipynb_checkpoints/`
   - `./tests/test_population_deposit/.ipynb_checkpoints/`
   - `./tests/test_ampere/.ipynb_checkpoints/`
   - `./tests/test_faraday/.ipynb_checkpoints/`
   - `./tests/boris/.ipynb_checkpoints/`
   - `./src/.ipynb_checkpoints/`

2. **Checkpoint Files** (Should be in .gitignore):
   - `.ipynb_checkpoints/test_hybirt-checkpoint.ipynb`
   - `.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - `.ipynb_checkpoints/plot_ampere-checkpoint.ipynb`
   - `.ipynb_checkpoints/plot_faraday-checkpoint.ipynb`
   - Multiple checkpoint files in test directories

### Recommended Actions:

1. **Add to .gitignore**:
   ```
   .ipynb_checkpoints/
   *-checkpoint.*
   *.ipynb_checkpoints
   ```

2. **Remove from repository**:
   - All `.ipynb_checkpoints` directories and their contents
   - All `*-checkpoint.*` files

3. **Consider adding**:
   ```
   build/
   CMakeCache.txt
   CMakeFiles/
   *.h5
   ```

### Other Observations:

- No CMake build artifacts found (good)
- No temporary test output files visible
- No compiled binaries committed (good)

## Summary

### Correct Implementations:
- ✅ Total density calculation (moments.hpp)
- ✅ Ampere calculation (with minor style issues)

### Implementations with Minor Issues:
- ⚠️ Faraday (different loop bounds, parameter order)
- ⚠️ Population deposit (overly complex but functionally correct)
- ⚠️ Boris pusher (missing safety checks)

### Superior Implementation:
- 🌟 Bulk velocity (better numerical stability with floor value)

### Critical Issues:
- ❌ ICN temporal integration has significant algorithmic differences
- ❌ Faraday parameter order differs from ground truth
- ❌ Repository contains many unwanted checkpoint files

### Overall Assessment:

The PR #4 implementation shows good understanding of the physics and algorithms. Most implementations are functionally correct, though some differ in details from the ground truth. The bulk velocity implementation is actually superior to the ground truth due to better numerical handling. However:

1. The ICN temporal integration needs careful review as it differs significantly from ground truth
2. The Faraday function signature inconsistency needs to be resolved
3. Repository cleanliness is poor with many checkpoint files
4. Code style could be improved to match ground truth conventions

**Recommendation**: The code needs refinement before merging, particularly around the temporal integration scheme and cleanup of checkpoint files.
