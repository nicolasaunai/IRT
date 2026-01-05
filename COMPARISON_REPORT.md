# Code Comparison Report: PR #15 vs hybirt Repository

## Executive Summary

This report compares the C++ and HPP code proposed in Pull Request #15 of the `nicolasaunai/IRT` repository with the corresponding files in the `nicolasaunai/hybirt` repository. The comparison focuses on implementation differences in five main source files and three test files.

**Date:** January 5, 2026  
**PR #15:** "HPC PROJECT - Bilosi Gaia" - Boris pusher implementation and validation  
**Repositories Compared:**
- **PR #15 Branch:** gaiabilosi:master (commit f18a52d)
- **Reference Repository:** nicolasaunai/hybirt (main/master branch)

---

## 1. ampere.hpp

### Overview
The Ampere class implements Ampère's law to calculate current density from magnetic field.

### Key Differences

#### PR #15 Implementation
- **Lines added:** 25 additions
- **Includes:** Added `#include "gridlayout.hpp"` (line 4)
- **Implementation approach:**
  - Sets `Jx = 0` everywhere on dual grid (lines 28-32)
  - Initializes boundary values for `Jy` and `Jz` at start to 0.0 (lines 39-40)
  - Uses backward difference for interior primal points (lines 42-48)
  - Formula: `J.y(ix) = -(B.z(ix) - B.z(ix-1))/dx` and `J.z(ix) = (B.y(ix) - B.y(ix-1))/dx`

#### hybirt Implementation
- **Implementation approach:**
  - No explicit `Jx` initialization
  - No boundary value initialization
  - Cleaner loop structure using references (lines 30-37)
  - Same backward difference formula but more concise
  - Operates directly on entire primal domain without special boundary handling

#### Analysis
**Similarities:**
- Both use backward finite differences on primal grid points
- Same mathematical formulas for `Jy` and `Jz`
- Both handle 1D case only

**Differences:**
1. **Boundary handling:** PR #15 explicitly sets boundary values to 0, hybirt does not
2. **Code organization:** hybirt uses local references (`auto&`) for cleaner syntax
3. **Jx handling:** PR #15 explicitly zeros out Jx on dual grid, hybirt omits this
4. **Include dependencies:** PR #15 adds explicit gridlayout.hpp include

**Quality Assessment:**
- hybirt version is more concise and cleaner
- PR #15's boundary initialization may be safer but could be overly conservative
- hybirt's use of references improves readability

---

## 2. faraday.hpp

### Overview
The Faraday class implements Faraday's law to update magnetic field from electric field.

### Key Differences

#### PR #15 Implementation
- **Lines added:** 51 additions (complete implementation from scratch)
- **Structure:**
  - Constructor with grid and dt parameters (lines 16-20)
  - Operator implementation with three-argument signature: `(B, E, Bnew)` (line 24)
  - Manual boundary handling: sets `Bnew.y(dend)` and `Bnew.z(dend)` to `B.y(dend)` and `B.z(dend)` (lines 47-48)
  - Loops only to `dend-1` to avoid out-of-bounds (line 51)
  - Forward difference: `Bnew.y(ix) = B.y(ix) + m_dt*(E.z(ix+1) - E.z(ix))/dx`
  - Separate loops for different field components

#### hybirt Implementation
- **Structure:**
  - Same constructor signature
  - Same operator signature
  - Cleaner loop structure using ghost domain (`ghost_start` to `ghost_end`)
  - Uses local references for field components (lines 34-42)
  - Forward difference with same formula
  - Handles all ghost points, no special boundary logic
  - Separate loops for By/Bz and Bx components

#### Analysis
**Similarities:**
- Same mathematical formulas (Faraday's law with forward differences)
- Same class structure and interface
- Both handle Bx separately (no evolution in 1D)

**Differences:**
1. **Domain handling:** 
   - PR #15 uses `dual_dom_start/end` with manual boundary exclusion
   - hybirt uses `ghost_start/ghost_end` for complete ghost coverage
2. **Code clarity:**
   - hybirt uses local references making code more readable
   - PR #15 has more explicit but verbose field access
3. **Boundary treatment:**
   - PR #15 explicitly copies boundary values and excludes last point
   - hybirt processes all ghost points uniformly

**Quality Assessment:**
- hybirt's ghost domain approach is more general and correct
- PR #15's manual boundary handling may miss edge cases
- hybirt version is cleaner and more maintainable

---

## 3. moments.hpp

### Overview
Functions to calculate plasma moments (total density and bulk velocity) from particle populations.

### Key Differences

#### total_density Function

**PR #15 Implementation (lines 20-24):**
```cpp
// TODO calculate the total density

// Accumulate density from each populaton
for(auto ix = 0; ix < N.data().size(); ++ix)
    N(ix) += pop.density()(ix);
```
- Simple accumulation loop
- Missing outer loop structure (indentation suggests it's inside the population loop)

**hybirt Implementation (lines 19-24):**
```cpp
for (auto const& pop : populations)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) += pop.density()(ix);
    }
}
```
- Explicit nested loop structure
- Clear iteration over all populations

#### bulk_velocity Function

**PR #15 Implementation (lines 47-57):**
```cpp
// TODO calculate bulk velocity by dividing by density N

// Normalize by total density to obtain mean velocity
for (auto ix = 0; ix < N.data().size(); ++ix) {
    if (N(ix) > 0.0) {
        // avoiding division by zero
        V.x(ix) /= N(ix);
        V.y(ix) /= N(ix);
        V.z(ix) /= N(ix);
    }
}
```
- Single loop with division by total density N
- Includes zero-check safety
- Correct algorithm: divides accumulated flux by total density

**hybirt Implementation (lines 47-56):**
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
- Nested loops over populations
- Divides by each population's density separately
- **INCORRECT:** This would divide multiple times by different densities
- No zero-check protection
- Contains commented debug line (line 44)

#### Analysis
**Similarities:**
- Both handle flux accumulation the same way (lines 34-42)
- Both initialize fields to zero

**Differences:**
1. **Correctness:**
   - **PR #15 bulk_velocity is CORRECT:** divides by total density N
   - **hybirt bulk_velocity is INCORRECT:** divides by individual population densities multiple times
2. **Safety:**
   - PR #15 includes division-by-zero protection
   - hybirt lacks this protection
3. **Code clarity:**
   - PR #15 has clearer comments
   - hybirt has debug code remnants

**Quality Assessment:**
- **PR #15 implementation is superior and mathematically correct**
- hybirt has a significant bug in bulk_velocity calculation
- PR #15 includes better safety checks

---

## 4. population.hpp

### Overview
The Population class manages particle data and deposits particle quantities onto grid.

### Key Differences

#### deposit() Function

**PR #15 Implementation (lines 104-125):**
```cpp
double const iCell_float = particle.position[0] / m_grid->cell_size(Direction::X);
int const iCell_         = static_cast<int>(iCell_float);
double const remainder    = iCell_float - iCell_;  // Fixed typo: reminder -> remainder
auto const iCell         = iCell_ + m_grid->dual_dom_start(Direction::X);

// TODO implement linear weighting deposit for the density and flux

// Linear weighting (1st order interpolation) to deposit particles onto nodes
// 2 nodes: iCell (left) iCell + 1 (right)
m_density(iCell) += particle.weight * (1.0 - remainder);
m_density(iCell + 1) += particle.weight * remainder;

m_flux.x(iCell) += particle.weight * particle.v[0] * (1.0 - remainder);
m_flux.x(iCell + 1) += particle.weight * particle.v[0] * remainder;

m_flux.y(iCell) += particle.weight * particle.v[1] * (1.0 - remainder);
m_flux.y(iCell + 1) += particle.weight * particle.v[1] * remainder;

m_flux.z(iCell) += particle.weight * particle.v[2] * (1.0 - remainder);
m_flux.z(iCell + 1) += particle.weight * particle.v[2] * remainder;
```
- Fixed variable name: `remainder` instead of `reminder`
- Separate deposit for each component
- Good comments explaining the algorithm

**hybirt Implementation (lines 104-122):**
```cpp
double const iCell_float = particle.position[0] / m_grid->cell_size(Direction::X);
int const iCell_         = static_cast<int>(iCell_float);
double const reminder    = iCell_float - iCell_;  // Typo: should be "remainder"
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
- Uses `reminder` (typo for "remainder")
- More compact notation: `weight * (1.0 - reminder) * v[i]`
- Same algorithm and logic

#### Analysis
**Similarities:**
- Identical algorithm for linear weighting
- Same structure and approach
- Both correctly implement first-order particle-to-grid interpolation

**Differences:**
1. **Variable naming:** PR #15 fixes typo (`remainder` vs `reminder`)
2. **Ordering:** PR #15 multiplies `weight * v * (1-r)`, hybirt uses `weight * (1-r) * v`
3. **Comments:** PR #15 has more explanatory comments
4. **Typo in comments:** PR #15 has "lest" instead of "left"

**Quality Assessment:**
- Functionally equivalent implementations
- PR #15 has better variable naming (fixes typo)
- Both could benefit from using local references for readability

---

## 5. pusher.hpp

### Overview
The Boris class implements the Boris particle pusher for integrating equations of motion.

### Key Differences

#### PR #15 Implementation (lines 48-108)
```cpp
// TODO implement the Boris pusher

double dt = this->dt_;

// Step 1: Half-step position update
particle.position[0] += 0.5 * particle.v[0] * dt;

// Step 2: Find iCell and remainder for interpolation
double x_norm = particle.position[0] / this->layout_->cell_size(Direction::X);
int iCell = static_cast<int>(x_norm) + this->layout_->dual_dom_start(Direction::X);
double remainder = x_norm - floor(x_norm);

// Step 3: Interpolate E and B fields at the particle position
double Ex = interpolate(E.x, iCell, remainder);
// ... (similar for Ey, Ez, Bx, By, Bz)

// Step 4: Half acceleration due to E field, defining v_minus
double qmdt2 = 0.5 * dt * (particle.charge / particle.mass);
double vx_minus = particle.v[0] + qmdt2 * Ex;
// ... (similar for vy_minus, vz_minus)

// Step 5: Rotation due to B field, defining t, s, v', v_plus
double tx = qmdt2 * Bx;
// ... (defines t, s vectors)
double t2 = tx * tx + ty * ty + tz * tz;
double sx = 2.0 * tx / (1.0 + t2);
// ... (similar for sy, sz)

// Cross products for rotation
double vpx = vx_minus + (vy_minus * tz - vz_minus * ty);
// ... (similar for vpy, vpz)
double vx_plus = vx_minus + (vpy * sz - vpz * sy);
// ... (similar for vy_plus, vz_plus)

// Step 6: Second half acceleration due to E
particle.v[0] = vx_plus + qmdt2 * Ex;
// ... (similar for v[1], v[2])

// Step 7: Second half-step position update
particle.position[0] += 0.5 * particle.v[0] * dt;
```
**Characteristics:**
- Detailed step-by-step comments
- Uses separate t and s vectors
- Position update split into two halves around velocity update
- Uses `floor()` for remainder calculation
- All vector components spelled out explicitly
- Uses `#include <cmath>` for floor function

#### hybirt Implementation (lines 48-111)
```cpp
// half step push position from t=n to t=n+1/2
// needed because fields are defined at t=n+1/2
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

// Interpolate E and B fields at particle position
// ... (similar interpolation)

// Calculate the half-step velocity
auto const vminus_x = particle.v[0] + qdto2m * ex;
// ... (similar for y, z)

auto const vprime_x = vminus_x + qdto2m * (vminus_y * bz - vminus_z * by);
// ... (similar cross product)

auto const s = 2 * qdto2m / (1 + qdto2m * qdto2m * (bx * bx + by * by + bz * bz));
auto const vplus_x = vminus_x + s * (vprime_y * bz - vprime_z * by);
// ... (similar cross product)

particle.v[0] = vplus_x + qdto2m * ex;
// ... (similar for v[1], v[2])

// velocity is now at t=n+1
// so we can compute the other half of the position update
for (auto dim = 0; dim < dimension; ++dim)
{
    auto const dr = particle.v[dim] * this->dt_ * 0.5;
    if (dr > this->layout_->cell_size(Direction::X) * 0.5)
    {
        throw std::runtime_error(
            "Particle moved more than half a cell size in one step in 2nd update");
    }
    particle.position[dim] += dr;
}
```
**Characteristics:**
- Generic dimension handling with loops
- Includes safety checks for particle displacement
- More compact mathematical notation
- Uses single `s` factor instead of separate s vector
- Variable naming: `reminder` instead of `remainder`
- More professional error handling

#### Analysis
**Similarities:**
- Both implement the same Boris algorithm
- Same mathematical steps: vminus → vprime → vplus
- Both use half-step position updates
- Identical interpolation approach

**Differences:**

1. **Algorithm formulation:**
   - **PR #15:** Uses separate t-vector and s-vector (more textbook)
   - **hybirt:** Uses combined qdto2m and single s scalar (more efficient)

2. **Code structure:**
   - **PR #15:** Hardcoded for 1D, explicit x/y/z components
   - **hybirt:** Generic loops over dimensions, more extensible

3. **Safety:**
   - **PR #15:** No displacement checks
   - **hybirt:** Validates particle doesn't move > 0.5 cell size

4. **Variable naming:**
   - **PR #15:** Uses `remainder`, clearer mathematical names (tx, ty, tz, sx, sy, sz)
   - **hybirt:** Uses `reminder` (typo), more compact names (qdto2m, s)
   - **Note:** hybirt also has a typo in error message: "1nd update" instead of "1st update"

5. **Position calculation:**
   - **PR #15:** `x_norm = pos/dx`, `iCell = int(x_norm) + start`, `remainder = x_norm - floor(x_norm)`
   - **hybirt:** `iCell_float = pos/dx + start`, `iCell = int(iCell_float)`, `reminder = iCell_float - iCell`

**Quality Assessment:**
- **hybirt is more robust:** includes safety checks, dimension-agnostic
- **PR #15 is more pedagogical:** clearer step-by-step comments, textbook formulation
- **hybirt is more production-ready:** better error handling, more maintainable
- Both are mathematically correct implementations of Boris algorithm

---

## 6. Test Files Comparison

### test_ampere.cpp

**PR #15 only:** New test file (151 lines)
- Tests Ampere's law with analytical functions: `By(x) = cos(x)`, `Bz(x) = sin(x)`
- Verifies numerical derivatives against analytical: `Jy = -cos(x)`, `Jz = -sin(x)`
- Uses tolerance proportional to `dx` (first-order accuracy check)
- Outputs results to HDF5 for visualization
- **Status:** Well-designed validation test

**hybirt:** No equivalent ampere test file found

### test_faraday.cpp

**PR #15 only:** New test file (176 lines)
- Tests Faraday's law with analytical E-fields: `Ey(x) = sin(x)`, `Ez(x) = cos(x)`
- Verifies incremental B-field updates: `ΔBy = -dt*sin(x)`, `ΔBz = -dt*cos(x)`
- Uses tolerance proportional to `dt*dx`
- Outputs to HDF5
- **Status:** Comprehensive validation test

**hybirt:** No equivalent faraday test file found

### test_population.cpp

**PR #15 only:** New test file (101 lines)
- Tests particle loading and deposition
- Validates mass conservation: sum(weights) = sum(deposited density)
- Tests moment calculations (total density, bulk velocity)
- Simple density profile (uniform density = 1.0)
- 5000 particles per cell
- **Status:** Good integration test

**hybirt:** No equivalent population test file found

### Additional Test Files in hybirt

**test_boris.cpp:**
- Comprehensive Boris pusher tests
- Two scenarios: uniform Bz field, and E×B drift test
- Validates energy conservation and drift velocities
- Outputs trajectory data to HDF5
- More extensive than PR #15 (which doesn't have a boris test)

**test_waves.cpp:**
- Tests wave propagation
- Not present in PR #15

**test_ionbeam.cpp:**
- Tests ion beam physics
- Not present in PR #15

---

## Summary of Key Findings

### 1. Implementation Quality

| File | PR #15 | hybirt | Winner |
|------|--------|--------|--------|
| ampere.hpp | More explicit boundary handling | Cleaner, more concise | hybirt (marginally) |
| faraday.hpp | Manual boundary handling | Ghost domain approach | hybirt |
| moments.hpp | **Correct algorithm** | **Incorrect bulk_velocity** | **PR #15** |
| population.hpp | Fixed typo, good comments | Typo in variable name | PR #15 (marginally) |
| pusher.hpp | Pedagogical, clear | Production-ready, robust | hybirt |

### 2. Major Differences

1. **moments.hpp:** PR #15 has the **correct** bulk_velocity implementation; hybirt has a **bug**
2. **pusher.hpp:** hybirt includes safety checks and dimension-agnostic code; PR #15 is more explicit
3. **faraday.hpp:** hybirt uses proper ghost domain handling; PR #15 uses manual boundaries
4. **Test coverage:** PR #15 adds new tests for ampere, faraday, and population; hybirt has boris test

### 3. Code Quality Observations

**PR #15 Strengths:**
- Correct moments calculation
- Good explanatory comments
- Fixes variable naming issues (remainder vs reminder)
- Comprehensive new test files

**PR #15 Weaknesses:**
- Less robust error handling in pusher
- Manual boundary handling in faraday could be error-prone
- Hardcoded for 1D only

**hybirt Strengths:**
- More production-ready code (error checks, extensibility)
- Cleaner code with references and better organization
- Ghost domain handling in faraday
- Dimension-agnostic loops where applicable

**hybirt Weaknesses:**
- **Critical bug in bulk_velocity calculation**
- Variable naming typos (reminder vs remainder)
- Some debug code remnants

### 4. Mathematical Correctness

- **Ampere:** Both correct, same algorithm
- **Faraday:** Both correct, hybirt has better implementation
- **Moments:** **PR #15 is correct, hybirt has a bug**
- **Population:** Both correct and equivalent
- **Boris Pusher:** Both correct, different styles

---

## Recommendations

### For PR #15 Review

1. **Accept moments.hpp** - It fixes a critical bug in hybirt
2. **Consider improvements to pusher.hpp:**
   - Add displacement safety checks from hybirt
   - Consider making it dimension-agnostic
3. **Consider improvements to faraday.hpp:**
   - Use ghost domain iteration instead of manual boundaries
4. **Test coverage:** Excellent addition of new tests
5. **Consider adding:** Boris pusher test similar to hybirt's test_boris.cpp

### For hybirt Repository

1. **Critical fix needed:** bulk_velocity in moments.hpp divides by wrong density
2. **Fix typos:** Use "remainder" consistently instead of "reminder"
3. **Consider adding:** Test files for ampere, faraday, and population

### For Code Merge Strategy

If merging from both sources:
- Take moments.hpp from PR #15 (correct algorithm)
- Take pusher.hpp structure from hybirt (add comments from PR #15)
- Take faraday.hpp approach from hybirt (ghost domains)
- Fix typos in all variable names
- Combine test suites from both

---

## Conclusion

PR #15 provides valuable implementations with **one critical correction** (moments.hpp) and excellent new test coverage. However, hybirt generally has more production-ready code structure with better error handling and extensibility. The ideal solution would combine:
- The correct moments calculation from PR #15
- The robust pusher implementation from hybirt
- The cleaner code organization from hybirt
- The comprehensive test suite from PR #15
- Fix all variable naming typos in both repositories

**Overall Assessment:** Both implementations have merit. PR #15 provides important bug fixes and test coverage, while hybirt provides more robust and maintainable code structure. Neither is strictly superior across all files.
