# PR #5 Code Review: Comparison with Ground Truth

This report analyzes the code implemented in PR #5 (PIC Nisrine Bakaev) by comparing it against the ground truth implementation from https://github.com/nicolasaunai/hybirt.

## 1. Ampere Implementation (src/ampere.hpp)

### Ground Truth Implementation
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

### PR #5 Implementation
```cpp
auto const& Bx = B.x;
auto const& By = B.y;
auto const& Bz = B.z;

auto& Jx = J.x;
auto& Jy = J.y;
auto& Jz = J.z;

// Ex, Jx are dual in x
for (auto ix = m_grid->dual_dom_start(Direction::X);
     ix <= m_grid->dual_dom_end(Direction::X); ++ix)
{ 
    Jx(ix) = 0.0; // since in 1D Jx = 0
}

// Ey, Ez, Jx, Jz are in primal. 
// Bz, By are in dual. 
for (auto ix = m_grid->primal_dom_start(Direction::X);
     ix <= m_grid->primal_dom_end(Direction::X); ++ix)
{
    Jy(ix) = -(Bz(ix) - Bz(ix-1)) / (dx);
    Jz(ix) =  (By(ix) - By(ix-1)) / (dx);
}
```

### Analysis
**Correctness:** ✅ **CORRECT** - The core computation for Jy and Jz matches the ground truth exactly.

**Quality:**
- **Positive:** The PR includes additional comments explaining the grid centering (dual vs primal), which is educational.
- **Positive:** Correctly handles Jx = 0 in 1D, which is physically accurate.
- **Minor Issue:** The PR extracts all field components (Bx, Jx) even though they're not all used in the computation. This is slightly verbose but not incorrect.
- **Minor Issue:** Unnecessary parentheses around `(dx)` in the division.

**Overall Quality:** Good - The implementation is correct with helpful comments, though slightly more verbose than necessary.

---

## 2. Faraday Implementation (src/faraday.hpp)

### Ground Truth Implementation
```cpp
template<std::size_t dimension>
class Faraday
{
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
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
                auto& Bnewx = Bnew.x;
                Bnewx(ix) = Bx(ix);
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;
};
```

### PR #5 Implementation
```cpp
template<std::size_t dimension>
class Faraday
{
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension> const& B, VecField<dimension>& Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            auto const& Ex = E.x;
            auto const& Ey = E.y;
            auto const& Ez = E.z;
            
            auto const& Bx = B.x;
            auto const& By = B.y;
            auto const& Bz = B.z;

            auto& Bnew_x = Bnew.x;
            auto& Bnew_y = Bnew.y;
            auto& Bnew_z = Bnew.z;
            
            // Bx is in primal
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                Bnew_x(ix) = Bx(ix); // Bx is constant in 1D
            }

            // By, Bz are in dual
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                Bnew_y(ix) = By(ix) + m_dt* (Ez(ix+1) - Ez(ix))/ (dx);
                Bnew_z(ix) = Bz(ix) - m_dt* (Ey(ix+1) - Ey(ix))/ (dx);
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }
private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;   
};
```

### Analysis
**Correctness:** ⚠️ **MOSTLY CORRECT with differences**

**Key Differences:**
1. **Parameter Order:** PR #5 has `operator()(E, B, Bnew)` while ground truth has `operator()(B, E, Bnew)`. This is a **significant API difference** that could cause confusion, though both are functionally equivalent.
2. **Loop Range:** PR #5 uses `primal_dom_start/end` and `dual_dom_start/end`, while ground truth uses `ghost_start/end`. This means PR #5 only updates domain points, not ghost cells, which **may cause issues with boundary conditions**.

**Quality:**
- **Positive:** Good comments explaining grid centering.
- **Issue:** Parameter order inconsistency with typical convention (fields usually come before output).
- **Issue:** Not updating ghost cells could lead to incorrect boundary behavior.
- **Minor Issue:** Inconsistent spacing around operators (`m_dt*` vs `m_dt *`).

**Overall Quality:** Fair - The core physics is correct, but the API signature difference and ghost cell handling could cause integration issues.

---

## 3. Boris Pusher Implementation (src/pusher.hpp)

### Ground Truth Implementation (abbreviated)
```cpp
// half step push position from t=n to t=n+1/2
for (auto dim = 0; dim < dimension; ++dim)
{
    auto const dr = particle.v[dim] * this->dt_ * 0.5;
    if (dr > this->layout_->cell_size(Direction::X) * 0.5)
        throw std::runtime_error("Particle moved more than half a cell size...");
    particle.position[dim] += dr;
}

double const iCell_float = particle.position[0] / this->layout_->cell_size(Direction::X)
                           + this->layout_->dual_dom_start(Direction::X);
int const iCell = static_cast<int>(iCell_float);
double const reminder = iCell_float - iCell;
double const qdto2m = particle.charge * this->dt_ / (2.0 * particle.mass);

// Interpolate E and B fields
double const ex = interpolate(E.x, iCell, reminder);
// ... other field interpolations ...

// Boris algorithm
auto const vminus_x = particle.v[0] + qdto2m * ex;
auto const vprime_x = vminus_x + qdto2m * (vminus_y * bz - vminus_z * by);
auto const s = 2 * qdto2m / (1 + qdto2m * qdto2m * (bx * bx + by * by + bz * bz));
auto const vplus_x = vminus_x + s * (vprime_y * bz - vprime_z * by);
particle.v[0] = vplus_x + qdto2m * ex;

// second half position update
for (auto dim = 0; dim < dimension; ++dim)
{
    particle.position[dim] += particle.v[dim] * this->dt_ * 0.5;
}
```

### PR #5 Implementation (abbreviated)
```cpp
auto dx = this->layout_->cell_size(Direction::X); 
auto dt = this->dt_;

particle.position[0] += particle.v[0]*dt/2;
        
auto iCell = static_cast<int>(particle.position[0] / dx);
auto remainder = (particle.position[0] / dx) - iCell;

auto Ex = interpolate(E.x, iCell, remainder);
// ... other field interpolations ...

auto factor = (particle.charge*dt)/(2*particle.mass);

auto vminus_x = particle.v[0] + factor * Ex;
// ... other vminus components ...

auto t_x = factor*Bx;
// ... other t components ...

auto vprime_x = vminus_x + (vminus_y * t_z - vminus_z * t_y);
// ... other vprime components ...

auto t_2 = t_x*t_x + t_y*t_y + t_z*t_z;

auto s_x = 2*t_x/(1+t_2);
// ... other s components ...

auto vplus_x = vminus_x + (vprime_y * s_z - vprime_z * s_y);
// ... other vplus components ...

particle.v[0] = vplus_x + factor*Ex;
// ... other velocity components ...

particle.position[0] += particle.v[0] * (dt/2);
```

### Analysis
**Correctness:** ❌ **INCORRECT - Missing critical offset**

**Critical Issues:**
1. **Missing iCell offset:** PR #5 uses `iCell = static_cast<int>(particle.position[0] / dx)` but the ground truth uses `iCell_float = particle.position[0] / dx + dual_dom_start(X)`. This **missing offset will cause incorrect field interpolation** at particle positions.
2. **No boundary checking:** Ground truth includes safety checks for particle motion, PR #5 does not.
3. **Different algorithm structure:** PR #5 uses separate `t` and `s` vectors instead of scalar `s`. While mathematically this could work, it's **unnecessarily complex** and deviates from the standard Boris algorithm.

**Quality:**
- **Issue:** The use of `t_x, t_y, t_z` and `s_x, s_y, s_z` vectors is non-standard. The standard Boris algorithm uses a single scalar `s`.
- **Issue:** Missing error checking for large particle motion.
- **Issue:** Variable naming inconsistency (`remainder` vs `reminder`).

**Overall Quality:** Poor - The implementation has a critical bug with the missing grid offset and uses a non-standard formulation of the Boris algorithm.

---

## 4. Moments Implementation (src/moments.hpp)

### Ground Truth Implementation
```cpp
template<std::size_t dimension>
void total_density(std::vector<Population<dimension>> const& populations, Field<dimension>& N)
{
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
}

template<std::size_t dimension>
void bulk_velocity(...) 
{
    // ... initialization ...
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
}
```

### PR #5 Implementation
```cpp
template<std::size_t dimension>
void total_density(std::vector<Population<dimension>> const& populations, Field<dimension>& N)
{
    // initialization already in master
    for (auto const& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            N(ix) += pop.density()(ix);
        }
    }
}

template<std::size_t dimension>
void bulk_velocity(...)
{
    // ... flux accumulation already in master ...
    
    for (std::size_t ix = 0; ix < N.data().size(); ++ix)
    {
        if (N(ix) > 0.0)
        {
            V.x(ix) /= N(ix);
            V.y(ix) /= N(ix);
            V.z(ix) /= N(ix);
        }
    }
}
```

### Analysis
**Correctness:** ⚠️ **DIFFERENT APPROACH**

**Key Differences:**
1. **total_density:** PR #5 is correct and matches ground truth.
2. **bulk_velocity:** 
   - **PR #5:** Divides by total density `N(ix)` with a safety check for zero density.
   - **Ground truth:** Divides by each population's density `pop.density()(ix)` - this appears to be a **bug in the ground truth** as it divides the accumulated flux multiple times.

**Quality:**
- **Positive:** PR #5's approach for bulk_velocity is actually **more correct** than the ground truth! The bulk velocity should be total flux / total density, not divided by each population separately.
- **Positive:** Includes safety check for division by zero.

**Overall Quality:** Excellent - PR #5's implementation is actually better than the ground truth for bulk_velocity calculation.

---

## 5. Population Deposit Implementation (src/population.hpp)

### Ground Truth Implementation
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

### PR #5 Implementation
```cpp
m_density(iCell + 1) += particle.weight * (1.0 - remainder);
m_density(iCell) += particle.weight * remainder;

m_flux.x(iCell + 1) += particle.weight * particle.v[0] * (1.0 - remainder);
m_flux.y(iCell + 1) += particle.weight * particle.v[1] * (1.0 - remainder);
m_flux.z(iCell + 1) += particle.weight * particle.v[2] * (1.0 - remainder);

m_flux.x(iCell) += particle.weight * particle.v[0] * remainder;
m_flux.y(iCell) += particle.weight * particle.v[1] * remainder;
m_flux.z(iCell) += particle.weight * particle.v[2] * remainder;
```

### Analysis
**Correctness:** ❌ **INCORRECT - Swapped weights**

**Critical Issue:**
The weights are **swapped**! 
- PR #5: `iCell+1` gets `(1.0 - remainder)`, `iCell` gets `remainder`
- Ground truth: `iCell` gets `(1.0 - remainder)`, `iCell+1` gets `remainder`

This means particles closer to `iCell` will deposit more weight to `iCell+1`, which is **backwards**.

**Quality:**
- **Issue:** Variable naming `remainder` vs `reminder` (inconsistent with codebase).
- **Issue:** The critical logic error will cause incorrect particle deposition.

**Overall Quality:** Poor - Contains a fundamental logic error that will produce incorrect results.

---

## 6. Main Loop Implementation (src/hybirt.cpp)

### 6.1 Average Function

**Ground Truth:**
```cpp
std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(),
               [](double const& f1, double const& f2) { return 0.5 * (f1 + f2); });
```

**PR #5:**
```cpp
std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(), 
    [](double a, double b) {return 0.5*(a+b);});
```

**Analysis:** ✅ **CORRECT** - Functionally identical, just minor style differences.

### 6.2 ICN Temporal Integration

**Ground Truth Structure:**
1. Predictor 1: Update fields, average, save particles, push and deposit
2. Predictor 2: Restart from saved particles, update fields with new moments
3. Corrector: Final field update

**PR #5 Structure:**
1. First Prediction: Update fields, average, push particles, compute moments
2. Second Predictor: Update fields, average, push particles again, compute moments
3. Corrector: Final field update

**Analysis:** ⚠️ **DIFFERENT APPROACH**

**Key Differences:**
- **Ground truth:** Saves particle state before predictor 1, restores before predictor 2
- **PR #5:** Pushes particles twice consecutively without restoration

The PR #5 approach **pushes particles for a full timestep** in predictor 1, then another full timestep in predictor 2, which is **incorrect** for the ICN scheme. The predictors should each advance by dt/2 from the original state.

**Overall Quality:** Poor - The ICN implementation is fundamentally flawed as it double-advances particles.

---

## 7. PR Cleanliness Analysis

### Unwanted Files Found in PR #5:

1. **`.ipynb_checkpoints/` directories** (multiple):
   - `notebooks/.ipynb_checkpoints/Test-Boris-1-checkpoint.ipynb`
   - `notebooks/.ipynb_checkpoints/Test-Faraday-checkpoint.ipynb`
   - `tests/ampere/.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - `tests/ampere/.ipynb_checkpoints/test_ampere-checkpoint.cpp`
   - `tests/boris/.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - `tests/boris/.ipynb_checkpoints/test_boris-checkpoint.cpp`
   - `tests/faraday/.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - `tests/faraday/.ipynb_checkpoints/test_faraday-checkpoint.cpp`
   - `tests/moments/.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - `tests/moments/.ipynb_checkpoints/test_moments-checkpoint.cpp`

2. **Binary/Data files:**
   - `notebooks/uniform_bz.h5` - HDF5 data file (binary output)

### Recommendations:
- **Add to `.gitignore`:**
  ```
  .ipynb_checkpoints/
  *.h5
  ```
- **Remove from repository:** All checkpoint files and binary data files should be removed from version control.

**Cleanliness Rating:** Poor - The PR contains many temporary/checkpoint files that should not be committed.

---

## Summary

| Component | Correctness | Quality | Notes |
|-----------|-------------|---------|-------|
| Ampere | ✅ Correct | Good | Verbose but correct with helpful comments |
| Faraday | ⚠️ Mostly Correct | Fair | API signature differs, missing ghost cell updates |
| Boris Pusher | ❌ Incorrect | Poor | Missing grid offset, non-standard formulation |
| Moments (total_density) | ✅ Correct | Good | Matches ground truth |
| Moments (bulk_velocity) | ✅ Better than GT | Excellent | More correct than ground truth! |
| Population Deposit | ❌ Incorrect | Poor | Weights are swapped |
| Average Function | ✅ Correct | Good | Functionally identical |
| ICN Integration | ❌ Incorrect | Poor | Fundamental flaw in particle advancement |
| PR Cleanliness | N/A | Poor | Contains checkpoint and binary files |

### Critical Issues to Fix:
1. **Boris Pusher:** Missing grid offset for iCell calculation
2. **Population Deposit:** Swapped interpolation weights  
3. **ICN Integration:** Incorrect particle advancement (double-stepping)
4. **Faraday:** Consider updating ghost cells, review API signature
5. **Repository Cleanliness:** Remove all .ipynb_checkpoints and .h5 files

### Strengths:
- Good commenting and documentation
- bulk_velocity implementation is actually superior to ground truth
- Ampere implementation is correct
