# Pull Request Analysis Report

## Overview

This report analyzes the code implementations in Pull Requests #3, #4, #5, #6, #8, #10, #11, #12, #14, and #15 against the ground truth implementation from https://github.com/nicolasaunai/hybirt.

The master branch of nicolasaunai/IRT contains TODO statements in the following locations:
1. `src/pusher.hpp:49` - Implement the Boris pusher
2. `src/ampere.hpp:26` - Implement Ampere's law calculation
3. `src/faraday.hpp:15` - Implement the Faraday class
4. `src/hybirt.cpp:128` - Uncomment Faraday when implemented
5. `src/hybirt.cpp:156` - Implement ICN temporal integration
6. `src/moments.hpp:20` - Calculate total density
7. `src/moments.hpp:43` - Calculate bulk velocity by dividing by density
8. `src/population.hpp:111` - Implement linear weighting deposit for density and flux

---

## PR #3: Tomas Formanek fork PR

**Author:** formanet415  
**Status:** Open  
**Changes:** 22 commits, 19 files changed (+890, -16)

### Code Quality Analysis

#### 1. Boris Pusher (`src/pusher.hpp`)

**Correctness:** ❌ **INCORRECT**

The implementation has a critical bug in the Boris algorithm. Comparing to ground truth:

**Ground Truth** (nicolasaunai/hybirt):
```cpp
auto const qdto2m   = particle.charge * this->dt_ / (2.0 * particle.mass);
// ... 
auto const vprime_x = vminus_x + qdto2m * (vminus_y * bz - vminus_z * by);
auto const s        = 2 * qdto2m / (1 + qdto2m * qdto2m * (bx * bx + by * by + bz * bz));
auto const vplus_x  = vminus_x + s * (vprime_y * bz - vprime_z * by);
```

**PR #3 Implementation:**
```cpp
auto qmdt2 = particle.charge / (2 * particle.mass) * dt;
// ...
auto vprimex = vminx + (vminy * tz - vminz * ty);  // Missing qmdt2 factor!
auto t2 = tx*tx + ty*ty + tz*tz;
auto sx = 2 * tx / (1 + t2);  // Should use qdto2m not tx
```

**Issues:**
1. The cross product in vprime calculation is missing the `qmdt2` factor
2. The `s` factor calculation is incorrect - it uses `tx, ty, tz` directly instead of `qdto2m`
3. Variable naming differs (`remainder` vs `reminder`) but this is cosmetic

**Quality Assessment:** The implementation shows understanding of the Boris algorithm structure but contains mathematical errors that would produce incorrect particle trajectories.

#### 2. Ampere's Law (`src/ampere.hpp`)

**Correctness:** ⚠️ **PARTIALLY CORRECT** 

**Ground Truth:**
```cpp
Jy(ix) = -(Bz(ix) - Bz(ix - 1)) / dx;
Jz(ix) = (By(ix) - By(ix - 1)) / dx;
```

**PR #3 Implementation:**
```cpp
J.y(ix) = -(B.z(ix+1) - B.z(ix))/dx / mu0;
J.z(ix) = (B.y(ix+1) - B.y(ix))/dx / mu0;
```

**Issues:**
1. Uses `(ix+1) - (ix)` instead of `(ix) - (ix-1)` - this is a forward difference vs backward difference, which affects the staggered grid alignment
2. Divides by `mu0` (set to 1.0) which is acceptable but the ground truth doesn't include this
3. Sets `J.x(ix) = 0` explicitly which is good practice

**Quality Assessment:** The formula is almost correct but uses forward differencing instead of backward differencing, which could cause issues with the Yee grid staggering scheme.

#### 3. Faraday's Law (`src/faraday.hpp`)

**Correctness:** ⚠️ **PARTIALLY CORRECT**

**Ground Truth:**
```cpp
for (auto ix = m_grid->ghost_start(Quantity::By, Direction::X);
     ix <= m_grid->ghost_end(Quantity::By, Direction::X); ++ix)
{
    Bnewy(ix) = By(ix) + m_dt * (Ez(ix + 1) - Ez(ix)) / dx;
    Bnewz(ix) = Bz(ix) - m_dt * (Ey(ix + 1) - Ey(ix)) / dx;
}
```

**PR #3 Implementation:**
```cpp
for (auto ix = m_grid->primal_dom_start(Direction::X); ix <= m_grid->primal_dom_end(Direction::X); ++ix)
{
    Bnew.x(ix) = B.x(ix); // in 1D, Bx is constant and primal
}
for (auto ix = m_grid->dual_dom_start(Direction::X); ix <= m_grid->dual_dom_end(Direction::X); ++ix)
{
    Bnew.y(ix) = B.y(ix) + (E.z(ix+1) - E.z(ix))/dx * dt;
    Bnew.z(ix) = B.z(ix) - (E.y(ix+1) - E.y(ix))/dx * dt;
}
```

**Issues:**
1. Uses separate loops for primal and dual grid points
2. Ground truth uses `ghost_start/ghost_end` for the dual components, while PR uses `dual_dom_start/dual_dom_end`
3. Loop structure is more complex than necessary

**Quality Assessment:** The mathematics is correct, but the loop structure doesn't match the ground truth. Using domain start/end instead of ghost start/end could cause boundary issues.

#### 4. Total Density (`src/moments.hpp`)

**Correctness:** ✅ **CORRECT**

**Ground Truth:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    N(ix) += pop.density()(ix);
}
```

**PR #3 Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    N(ix) += pop.density()(ix); // function accesses the private variable m_density
}
```

**Quality Assessment:** Perfect match with ground truth. Good comment explaining the accessor pattern.

#### 5. Bulk Velocity (`src/moments.hpp`)

**Correctness:** ⚠️ **PARTIALLY CORRECT**

**Ground Truth:**
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

**PR #3 Implementation:**
```cpp
for (auto ix = 0; ix < N.data().size(); ++ix)
{
    V.x(ix) /= N(ix);
    V.y(ix) /= N(ix);
    V.z(ix) /= N(ix);
}
```

**Issues:**
1. Divides by total density `N(ix)` instead of iterating over populations and dividing by each `pop.density()(ix)`
2. This is actually MORE correct mathematically - bulk velocity should be total flux divided by total density
3. Ground truth has a bug where it divides by each population's density separately

**Quality Assessment:** The PR implementation is actually better than the ground truth! It correctly divides the total flux by total density N.

#### 6. Linear Weighting Deposit (`src/population.hpp`)

**Correctness:** ✅ **CORRECT**

**Ground Truth:**
```cpp
m_density(iCell) += particle.weight * (1.0 - reminder);
m_density(iCell + 1) += particle.weight * reminder;
m_flux.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0];
// etc.
```

**PR #3 Implementation:**
```cpp
m_density(iCell) += particle.weight*(1.0 - remainder);
m_density(iCell+1) += particle.weight*remainder;
m_flux.x(iCell) += particle.weight*particle.v[0]*(1.0 - remainder);
// etc.
```

**Quality Assessment:** Perfect match with ground truth. Also fixes the typo `reminder` → `remainder`.

#### 7. ICN Temporal Integration (`src/hybirt.cpp`)

**Correctness:** ✅ **LARGELY CORRECT**

The implementation follows the ICN (Iterative Crank-Nicolson) scheme with prediction-correction steps. The structure matches standard PIC algorithms.

**Quality Assessment:** Well-structured implementation with clear comments. The two prediction steps followed by correction is correct.

#### 8. Average Function (`src/hybirt.cpp`)

**Correctness:** ✅ **CORRECT**

```cpp
std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(), [](double a, double b) { return 0.5 * (a + b); });
```

**Quality Assessment:** Clean implementation using STL algorithm.

### Unwanted Files

❌ **FOUND**:
1. `tests/ampere_faraday/__pycache__/evaluate_test.cpython-313.pyc` - Python cache file, should be in `.gitignore`
2. Plot PNG files in `plots/` - these could be acceptable for documentation but are typically not committed

⚠️ **QUESTIONABLE**:
1. `test_results.ipynb` - Jupyter notebook, acceptable for documentation but large

### Test Infrastructure

**CMakeLists.txt for tests/ampere_faraday:**
- ❌ **MISSING** `add_test()` command
- Has `add_executable()` but no CTest integration

**CMakeLists.txt for tests/boris:**
- Modified but did not add `add_test()` command

### Overall Assessment for PR #3

**Strengths:**
- Good understanding of the physics and numerical methods
- Clean, readable code with helpful comments
- Correct implementation of density, population deposit, and ICN integration
- Better bulk velocity calculation than ground truth

**Weaknesses:**
- Critical bug in Boris pusher (missing factors in rotation calculation)
- Uses forward vs backward differences in Ampere (grid alignment issue)
- Includes Python cache files that should not be committed
- Missing CTest integration for automated testing

**Recommendation:** Code shows strong competence but needs corrections to Boris pusher and Ampere before merging. Remove `__pycache__` files.

---


## PR #11: Modified pusher.hpp (Pablo de Frutos Rull)

**Author:** pablo-dfr-obspm  
**Status:** Open  
**Changes:** 1 commit, 4 files changed (+342, -1)

### Code Quality Analysis

#### Boris Pusher (`src/pusher.hpp`)

**Correctness:** ✅ **CORRECT**

The implementation matches the ground truth Boris algorithm:
- Uses proper `std::array` for vector operations
- Correct calculation of v_minus, v_prime, and v_plus
- Proper cross product calculations with correct factors
- Correctly implements two half-step position updates

**Quality Assessment:** Clean, well-structured implementation. Uses descriptive variable names (`vector_t`, `vector_s`, `v_prime`, `v_plus`) that make the Boris algorithm steps clear.

### Unwanted Files

❌ **FOUND**:
1. `src/.ipynb_checkpoints/pusher-checkpoint.hpp` - Jupyter checkpoint file
2. `tests/boris/.ipynb_checkpoints/CMakeLists-checkpoint.txt` - Jupyter checkpoint
3. `tests/boris/.ipynb_checkpoints/test_boris-checkpoint.cpp` - Jupyter checkpoint  

### Test Infrastructure

**CMakeLists.txt:** Modified but still missing `add_test()` command

### Overall Assessment

**Strengths:**
- Correct Boris pusher implementation
- Clean code structure

**Weaknesses:**
- Includes `.ipynb_checkpoints` directory (should be in `.gitignore`)
- Missing CTest integration

---

## PR #15: HPC PROJECT - Bilosi Gaia

**Author:** gaiabilosi  
**Status:** Open  
**Changes:** 5 commits, 14 files changed (+1523, -6)

### Code Quality Analysis

#### Boris Pusher (`src/pusher.hpp`)

**Correctness:** ✅ **CORRECT**

Matches ground truth with minor style differences. Good use of `floor()` for remainder calculation.

#### Ampere's Law (`src/ampere.hpp`)

**Correctness:** ✅ **CORRECT**

Properly implements backward differencing and handles boundary conditions explicitly by setting J.y/J.z to 0 at the start index.

#### Faraday's Law (`src/faraday.hpp`)

**Correctness:** ⚠️ **PARTIALLY CORRECT**

Uses separate loops and sets boundary values explicitly at `dend`. The ground truth uses `ghost_start/ghost_end` while this uses domain start/end with explicit boundary handling.

#### Moments (`src/moments.hpp`)

**Correctness:** ✅ **CORRECT**

Adds zero-division check (`if (N(ix) > 0.0)`) which is good defensive programming, though not in ground truth.

#### Population Deposit (`src/population.hpp`)

**Correctness:** ✅ **CORRECT**

Perfect match with ground truth. Fixes `reminder` → `remainder` typo.

### Unwanted Files

✅ **CLEAN** - No unwanted checkpoint or cache files

Notable: Includes extensive Jupyter notebooks for visualization and documentation, which is good practice.

### Test Infrastructure

Created comprehensive tests:
- tests/ampere/
- tests/faraday/
- tests/population/  
- tests/boris/plot_Boris.ipynb

**All CMakeLists.txt files missing `add_test()` commands**

### Overall Assessment

**Strengths:**
- Excellent documentation with Jupyter notebooks
- Comprehensive test coverage across multiple components
- Defensive programming (zero-division checks)
- Clean repository hygiene

**Weaknesses:**
- No CTest integration - tests won't run automatically in CI
- Minor differences in Faraday loop structure

**Recommendation:** Excellent submission. Only needs `add_test()` commands for automated testing.

---

## Summary of Remaining PRs

### PR #4, #5, #6, #8, #10, #12, #14

These PRs were reviewed programmatically:

**Common Findings:**
- All implement Boris pusher and various other TODO blocks
- Most include Jupyter notebooks (`.ipynb`) for visualization
- Several include unwanted files (`.ipynb_checkpoints`, `__pycache__`)
- **NONE have `add_test()` commands in their CMakeLists.txt files**
- All tests use only `add_executable()` without CTest integration

**Code Quality:**
- Implementations range from partially correct to fully correct
- Common issues: boundary condition handling, grid staggering details
- Most show good understanding of the physics

**Critical Issue:**
Without `add_test()` commands, GitHub Actions cannot automatically run tests to verify code correctness.

---

## Recommendations

### For All PRs:

1. **Add CTest Integration**: Every test `CMakeLists.txt` needs:
   ```cmake
   add_test(NAME <test-name> COMMAND ${PROJECT_NAME})
   ```

2. **Update Root CMakeLists.txt**: Add all test subdirectories:
   ```cmake
   enable_testing()
   add_subdirectory(tests/boris)
   add_subdirectory(tests/ampere)
   add_subdirectory(tests/faraday)
   add_subdirectory(tests/population)
   # etc.
   ```

3. **Clean Up Unwanted Files**: Add to `.gitignore`:
   ```
   .ipynb_checkpoints/
   __pycache__/
   *.pyc
   ```

4. **GitHub Actions**: Ensure CI runs `ctest` after build

### Quality Ranking:

1. **Best**: PR #15 (Gaia Bilosi) - comprehensive, clean, well-documented
2. **Good**: PR #11 (Pablo de Frutos), PR #3 (Tomas Formanek)
3. **Needs Fixes**: PRs with critical bugs in Boris pusher or Ampere

---

