# Code Review: PR #12 - HPC_Syrus_Michael

This document provides a detailed analysis of the code implementations in PR #12 compared to the ground truth implementation from the [hybirt repository](https://github.com/nicolasaunai/hybirt).

## Overview

PR #12 implements missing code blocks in the IRT repository. The missing code was identified by TODO comments in the master branch. This review compares each implementation against the ground truth to assess correctness and code quality.

---

## 1. Ampere's Law Implementation (`src/ampere.hpp`)

### Ground Truth Implementation
```cpp
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

### PR #12 Implementation
```cpp
for (int ix=1; ix< B.y.data().size() -1 ; ++ix)
{
    J.x(ix) = 0; 
    J.y(ix) = -(B.z(ix+1) - B.z(ix-1)) / (2.0*dx);
    J.z(ix) = (B.y(ix+1) - B.y(ix-1)) / (2.0*dx);
}

// Boundary Conditions
J.x(0) = J.x(B.y.data().size()-1) = 0;
J.y(0) = J.y(B.y.data().size()-1) = 0;
J.z(0) = J.z(B.z.data().size()-1) = 0;
```

### Analysis

**Correctness: ❌ INCORRECT**

The PR #12 implementation has several critical errors:

1. **Wrong discretization scheme**: Uses centered differences `(ix+1) - (ix-1) / (2.0*dx)` instead of backward differences `(ix) - (ix-1) / dx`
2. **Wrong loop bounds**: Uses `1 to size-1` instead of using the grid's `primal_dom_start/end` methods
3. **Unnecessary boundary condition code**: The grid layout should handle boundaries through its domain specifications
4. **Sets J.x to 0**: While J.x should be 0 in 1D, this is redundant if the field is properly initialized

**Quality Issues:**
- Poor spacing and formatting (e.g., `ix=1` should be `ix = 1`)
- Uses raw array indexing instead of grid layout methods
- Incorrect numerical discretization violates Maxwell's equations on a staggered grid

---

## 2. Faraday's Law Implementation (`src/faraday.hpp`)

### Ground Truth Implementation
```cpp
Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
    : m_grid{grid}, m_dt{dt}

void operator()(VecField<dimension> const& B, VecField<dimension> const& E,
                VecField<dimension>& Bnew)
{
    for (auto ix = m_grid->ghost_start(Quantity::By, Direction::X);
         ix <= m_grid->ghost_end(Quantity::By, Direction::X); ++ix)
    {
        Bnewy(ix) = By(ix) + m_dt * (Ez(ix + 1) - Ez(ix)) / dx;
        Bnewz(ix) = Bz(ix) - m_dt * (Ey(ix + 1) - Ey(ix)) / dx;
    }
    for (auto ix = m_grid->ghost_start(Quantity::Bx, Direction::X);
         ix <= m_grid->ghost_end(Quantity::Bx, Direction::X); ++ix)
    {
        Bnewx(ix) = Bx(ix);
    }
}
```

### PR #12 Implementation
```cpp
Faraday(std::shared_ptr<GridLayout<dimension>> grid)
    : m_grid{grid}

void operator()(VecField<dimension> const& B, VecField<dimension> const& E, 
                VecField<dimension>& Bnew, double dt)
{
    auto N = B.y.data().size();
    for (int ix=1; ix< N -1 ; ++ix)
    {
        Bnew.x(ix) = B.x(ix);
        Bnew.y(ix) = B.y(ix) - dt * (-(E.z(ix+1) - E.z(ix-1)) / (2.0 * dx));
        Bnew.z(ix) = B.z(ix) - dt * ((E.y(ix+1) - E.y(ix-1)) / (2.0 * dx));
    }

    // Boundary Conditions
    Bnew.x(0) = B.x(0);
    Bnew.y(0) = B.y(0);
    Bnew.z(0) = B.z(0);
    Bnew.x(N-1) = B.x(N-1);
    Bnew.y(N-1) = B.y(N-1);
    Bnew.z(N-1) = B.z(N-1);
}
```

### Analysis

**Correctness: ❌ INCORRECT**

Major issues:

1. **Wrong constructor**: `dt` should be stored as a member variable, not passed to the operator
2. **Wrong discretization**: Uses centered differences `(ix+1) - (ix-1) / (2.0*dx)` instead of forward differences `(ix+1) - (ix) / dx`
3. **Wrong loop bounds**: Should use `ghost_start/ghost_end` for proper staggered grid handling
4. **Redundant boundary code**: Properly designed grid layout handles this
5. **Incorrect sign and parentheses**: The double negative in `-(-(E.z...))` is confusing and suggests misunderstanding

**Quality Issues:**
- Inconsistent with the class design pattern (dt should be a member)
- Poor code organization (all components updated in one loop vs. separate loops)
- Formatting issues

---

## 3. Boris Pusher Implementation (`src/pusher.hpp`)

### Ground Truth Implementation
```cpp
// Position update (half step)
for (auto dim = 0; dim < dimension; ++dim)
{
    auto const dr = particle.v[dim] * this->dt_ * 0.5;
    if (dr > this->layout_->cell_size(Direction::X) * 0.5)
        throw std::runtime_error("Particle moved more than half a cell");
    particle.position[dim] += dr;
}

// Calculate iCell and reminder
double const iCell_float = particle.position[0] / this->layout_->cell_size(Direction::X)
                           + this->layout_->dual_dom_start(Direction::X);
int const iCell = static_cast<int>(iCell_float);
double const reminder = iCell_float - iCell;
double const qdto2m = particle.charge * this->dt_ / (2.0 * particle.mass);

// Interpolate fields and apply Boris algorithm
// ... (uses 's' formulation with proper vector algebra)
```

### PR #12 Implementation
```cpp
particle.position[0] += particle.v[0] * 0.5 * dt;
double dx = this->layout_->cell_size(Direction::X);

int iCell = static_cast<int>(particle.position[0]/dx);
double remainder = particle.position[0]/dx - iCell;

// Interpolate E and B
// ... (correctly done)

double qdtm = (particle.charge * dt) / (2 * particle.mass);

// Boris algorithm using t and s vectors
// ... (uses separate s_x, s_y, s_z instead of vectorized s)

particle.position[0] += particle.v[0] * dt/2;
```

### Analysis

**Correctness: ⚠️ MOSTLY CORRECT with issues**

Issues identified:

1. **Missing domain offset**: iCell calculation doesn't add `dual_dom_start(Direction::X)` like the ground truth
2. **No safety check**: Missing the particle displacement validation that prevents particles from moving more than half a cell
3. **Only updates position[0]**: Doesn't use a loop for all dimensions (though this is acceptable for 1D)
4. **Redundant dt declaration**: `double dt = this->dt_` is unnecessary

**Quality Issues:**
- Less robust than ground truth (missing safety checks)
- Formatting inconsistencies (spacing around operators)
- Uses expanded s_x, s_y, s_z instead of vectorized approach
- The variable naming is good (remainder vs reminder)

**Positive aspects:**
- The core Boris algorithm logic is correct
- Interpolation is done correctly
- The mathematical steps are in the right order

---

## 4. Particle Deposit Implementation (`src/population.hpp`)

### Ground Truth Implementation
```cpp
double const reminder = iCell_float - iCell_;

m_density(iCell) += particle.weight * (1.0 - reminder);
m_density(iCell + 1) += particle.weight * reminder;

m_flux.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0];
m_flux.x(iCell + 1) += particle.weight * reminder * particle.v[0];
// ... (similar for y and z)
```

### PR #12 Implementation
```cpp
double const remainder = iCell_float - iCell_;

double w_left = 1.0-remainder;
double w_right = remainder;

m_density(iCell) += particle.weight * w_left;
m_density(iCell+1) += particle.weight * w_right;

m_flux.x(iCell) += particle.weight * particle.v[0] * w_left;
m_flux.x(iCell+1) += particle.weight * particle.v[0] * w_right;
// ... (similar for y and z)
```

### Analysis

**Correctness: ✅ CORRECT**

The implementation is functionally correct. The linear weighting is properly implemented.

**Quality Issues:**
- **Unnecessary intermediate variables**: `w_left` and `w_right` add no value and reduce performance slightly
- **Minor formatting**: `1.0-remainder` should have spaces; `iCell+1` should be `iCell + 1`
- **Verbose**: Ground truth is more concise and equally readable

**Positive aspects:**
- Algorithm is correct
- Clear variable names (remainder vs reminder)
- Proper linear interpolation weights

---

## 5. Moments Calculation (`src/moments.hpp`)

### Ground Truth Implementation (`total_density`)
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

### PR #12 Implementation (`total_density`)
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    N(ix) = 0;
}
for (auto const& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
        N(ix) += pop.density()(ix);
}
```

### Analysis

**Correctness: ✅ CORRECT**

The total density implementation is correct and nearly identical to the ground truth.

**Quality Issues:**
- Missing braces around the inner for loop (minor style issue)

---

### Ground Truth Implementation (`bulk_velocity`)
```cpp
// Initialize V to zero
// Sum all fluxes
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

### PR #12 Implementation (`bulk_velocity`)
```cpp
// Initialize V to zero
// Sum all fluxes (correctly done)
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    if (N(ix) > 0)
    {
        V.x(ix) /= N(ix);
        V.y(ix) /= N(ix);
        V.z(ix) /= N(ix);
    }
    else
    {
        V.x(ix) = 0;
        V.y(ix) = 0;
        V.z(ix) = 0;
    }
}
```

### Analysis

**Correctness: ⚠️ DIFFERENT but arguably BETTER**

The PR implementation differs from ground truth:
- **PR divides by total density N(ix)** - This is the physically correct bulk velocity formula: V = Σ(n_i * v_i) / Σ(n_i)
- **Ground truth divides by each population's density** - This appears to be incorrect for multi-species plasmas

**Quality Issues:**
- The zero-check is good defensive programming
- However, if V is already zeroed, the else clause is redundant

**Verdict: The PR implementation is actually MORE CORRECT than the ground truth for this particular function.**

---

## 6. Average Function (`src/hybirt.cpp`)

### Ground Truth Implementation
```cpp
void average(Field<dimension> const& F1, Field<dimension> const& F2, Field<dimension>& Favg)
{
    std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(),
                   [](double const& f1, double const& f2) { return 0.5 * (f1 + f2); });
}
```

### PR #12 Implementation
```cpp
void average(Field<dimension> const& F1, Field<dimension> const& F2, Field<dimension>& Favg)
{
    std::transform(F1.data().begin(), F1.data().end(), F2.data().begin(), F1.data().begin(), 
                   [](double a, double b) { return 0.5 * (a + b); });
    // use std::transform to do an average of F1 and F2
}
```

### Analysis

**Correctness: ❌ INCORRECT**

Critical error:
1. **Wrong output iterator**: Uses `F1.data().begin()` instead of `Favg.begin()` - this modifies F1 instead of Favg!
2. **Wrong accessor method**: Uses `.data().begin()` instead of `.begin()` (may work but not idiomatic)

**Quality Issues:**
- This is a serious bug that would corrupt the input field
- Lambda parameter types are less const-correct

---

## 7. ICN Temporal Integration (`src/hybirt.cpp`)

### Ground Truth Implementation
```cpp
// Predictor 1
faraday(B, E, Bnew);
ampere(Bnew, J);
boundary_condition->fill(J);
ohm(Bnew, J, N, V, Enew);
average(E, Enew, Eavg);
average(B, Bnew, Bavg);
boundary_condition->fill(Eavg);
boundary_condition->fill(Bavg);

auto save{populations};  // Save particles at t=n
// Push particles and deposit
// ...

// Predictor 2
faraday(B, Eavg, Bnew);
// ... similar structure

// Corrector
faraday(B, Eavg, B);  // Note: updates B in-place
ampere(B, J);
// ...
```

### PR #12 Implementation
```cpp
// Prediction 1
faraday(B, E, Bnew, dt);
boundary_condition->fill(Bnew);
ampere(Bnew, J);
boundary_condition->fill(J);
ohm(Bnew, J, N, V, Enew);
boundary_condition->fill(Enew);
average(B, Bnew, Bavg);
average(E, Enew, Eavg);
// Push particles
// ...

// Prediction 2
faraday(B, Eavg, Bnew, dt);
// ... similar

// Correction
faraday(B, Eavg, Bnew, dt);  // Note: updates Bnew, not B
// ...
```

### Analysis

**Correctness: ⚠️ PARTIALLY CORRECT with issues**

Issues:
1. **Corrector step is wrong**: Updates `Bnew` instead of `B` in-place, breaking the algorithm
2. **Extra boundary fills**: Fills `Bnew` and `Enew` unnecessarily after faraday and ohm
3. **Missing particle save/restore**: Ground truth saves particles before predictor 1 to restart from t=n for predictor 2; PR doesn't do this
4. **Different call signature**: Passes `dt` to faraday (which is inconsistent with best practice but matches PR's Faraday implementation)

**Quality Issues:**
- The general structure is similar but key details are wrong
- The corrector step error is critical and will produce incorrect results

---

## 8. PR Cleanliness Assessment

### Files Added/Modified

The PR includes:
- Modified source files (appropriate)
- Modified `.gitignore` (appropriate - adds `*.ipynb_checkpoints`)
- **Added `tests/boris/Plotting.ipynb` (271 KB)** - ⚠️ **Should NOT be in the PR**

### Issues

1. **Jupyter Notebook in tests directory**: The file `tests/boris/Plotting.ipynb` is a large Jupyter notebook (271 KB) that contains plotting/analysis code. This should not be committed to the repository.
   - Notebooks often contain binary output data that bloats the repository
   - Personal analysis notebooks should be kept local or in a separate documentation/examples directory

2. **The `.gitignore` addition is incomplete**: While it adds `*.ipynb_checkpoints`, it doesn't prevent `.ipynb` files from being committed in tests directories.

### Recommendations

1. Remove `tests/boris/Plotting.ipynb` from the repository
2. Consider adding `*.ipynb` to `.gitignore` for test directories, or create a separate directory for analysis notebooks
3. No CMake temporary files or other build artifacts were found - this is good

---

## 9. Test Configuration Analysis

### Current State

Examining the test CMakeLists.txt files:

- **`tests/boris/CMakeLists.txt`**: ✅ Has `add_test(NAME test-boris COMMAND test-boris)` on line 10
- **`tests/waves/CMakeLists.txt`**: ❌ Missing `add_test()` command
- **`tests/ion_beam/CMakeLists.txt`**: ❌ Missing `add_test()` command

### Root CMakeLists.txt

The root `CMakeLists.txt` includes:
- Line 77: `add_subdirectory(tests/boris)` ✅
- Missing: `add_subdirectory(tests/waves)` ❌
- Missing: `add_subdirectory(tests/ion_beam)` ❌

### Required Fixes

1. Add `add_test()` commands to waves and ion_beam test CMakeLists
2. Add test subdirectories to root CMakeLists.txt
3. These changes should be in a new PR with "AUTOMATIC" in the title

---

## Summary

### Code Quality by Module

| Module | Correctness | Quality | Notes |
|--------|-------------|---------|-------|
| Ampere | ❌ Incorrect | Poor | Wrong discretization scheme, wrong bounds |
| Faraday | ❌ Incorrect | Poor | Wrong discretization, wrong design |
| Boris Pusher | ⚠️ Mostly OK | Fair | Missing safety checks and domain offset |
| Deposit | ✅ Correct | Good | Works but verbose |
| Total Density | ✅ Correct | Good | Nearly perfect |
| Bulk Velocity | ✅ Better than GT | Good | Actually more correct than ground truth |
| Average | ❌ Incorrect | Poor | Critical bug - modifies input |
| ICN Integration | ⚠️ Partial | Fair | Corrector step is wrong |

### Critical Issues

1. **Ampere and Faraday use wrong discretization**: These will produce incorrect physics
2. **Average function has a serious bug**: Modifies input instead of output
3. **ICN corrector step is wrong**: Won't converge properly
4. **Jupyter notebook should be removed**: Bloats repository unnecessarily

### Recommendations

**For the student:**
1. Study staggered grid discretization (Yee grid) for electromagnetic codes
2. Understand the difference between forward, backward, and centered differences
3. Review the ICN (Iterative Crank-Nicolson) algorithm carefully
4. Use proper grid layout methods instead of raw array indexing
5. Test the code against analytical solutions

**For the repository:**
1. The PR should not be merged as-is due to critical correctness issues
2. After fixes, enable all tests in CI/CD
3. Add unit tests for each operator (Ampere, Faraday, etc.)
