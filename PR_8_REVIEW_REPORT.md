# PR #8 Code Review Report

## Overview
This report reviews PR #8 ("paulv") which implements missing code blocks in the IRT repository. The ground truth reference is taken from https://github.com/nicolasaunai/hybirt.

## 1. Boris Pusher Implementation (`src/pusher.hpp`)

### Ground Truth Implementation
The correct Boris pusher implementation from hybirt includes:
- Half-step position update to t=n+1/2
- Field interpolation at particle position
- Proper Boris algorithm with vminus, vprime, vplus calculations
- Correct formula for magnetic rotation: `s = 2*qdto2m / (1 + qdto2m^2 * (bx^2 + by^2 + bz^2))`
- Safety checks for particle movement
- Second half-step position update

### PR #8 Implementation Analysis

#### Correctness: ❌ **INCORRECT**

The PR #8 implementation has several critical errors:

1. **Incorrect Magnetic Rotation Formula**
   - **Ground Truth**: `s = 2*qdto2m / (1 + qdto2m^2 * (bx^2 + by^2 + bz^2))`
   - **PR #8**: `s[i] = (2.0*t[i])/(1+t[i]*t[i])` for each component separately
   - **Issue**: The s vector should use a single scalar denominator that includes all three magnetic field components. PR #8 incorrectly computes s per component, which is physically wrong.

2. **Missing Safety Checks**
   - Ground truth includes checks to ensure particles don't move more than half a cell in one step
   - PR #8 lacks these safety checks

3. **Position Update Logic**
   - Ground truth properly updates position in a loop over dimensions
   - PR #8 only updates x-direction (position[0])

4. **Field Interpolation**
   - Ground truth: Uses cleaner variable names and better structure
   - PR #8: Acceptable but less clean

#### Quality: **Poor**

- Missing critical physics correctness (magnetic rotation formula is wrong)
- Incomplete implementation (only 1D position updates)
- No safety validation
- Inconsistent with ground truth algorithm

---

## 2. Missing TODO Blocks Not Implemented in PR #8

The following TODO blocks were **NOT** addressed in PR #8:

### 2.1 Ampere Implementation (`src/ampere.hpp:26`)
- **Status**: Not implemented in PR #8
- **Expected**: Loop over primal grid points computing Jy and Jz from curl of B
- **Ground Truth**:
  ```cpp
  for (auto ix = m_grid->primal_dom_start(Direction::X);
       ix <= m_grid->primal_dom_end(Direction::X); ++ix)
  {
      Jy(ix) = -(Bz(ix) - Bz(ix - 1)) / dx;
      Jz(ix) = (By(ix) - By(ix - 1)) / dx;
  }
  ```

### 2.2 Total Density Calculation (`src/moments.hpp:20`)
- **Status**: Not implemented in PR #8
- **Expected**: Loop to accumulate density from all populations
- **Ground Truth**:
  ```cpp
  for (auto ix = 0; ix < N.data().size(); ++ix)
  {
      N(ix) += pop.density()(ix);
  }
  ```

### 2.3 Bulk Velocity Calculation (`src/moments.hpp:43`)
- **Status**: Not implemented in PR #8
- **Expected**: Divide flux by density to get bulk velocity
- **Ground Truth**:
  ```cpp
  for (auto ix = 0; ix < N.data().size(); ++ix)
  {
      V.x(ix) /= pop.density()(ix);
      V.y(ix) /= pop.density()(ix);
      V.z(ix) /= pop.density()(ix);
  }
  ```

### 2.4 Linear Weighting Deposit (`src/population.hpp:111`)
- **Status**: Not implemented in PR #8
- **Expected**: Particle-in-cell linear weighting for density and flux deposition
- **Ground Truth**:
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

### 2.5 ICN Temporal Integration (`src/hybirt.cpp:156`)
- **Status**: Not implemented in PR #8
- **Expected**: Iterative Crank-Nicolson scheme with predictor-corrector
- **Ground Truth**: Complex multi-step scheme including:
  - Predictor 1: Update B and E with Faraday and Ampere
  - Push particles with averaged fields
  - Predictor 2: Repeat with updated fields
  - Corrector: Final update

### 2.6 Average Function (`src/hybirt.cpp:28`)
- **Status**: Not implemented in PR #8
- **Expected**: Average two fields element-wise
- **Ground Truth**:
  ```cpp
  std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(),
                 [](double const& f1, double const& f2) { return 0.5 * (f1 + f2); });
  ```

---

## 3. PR Cleanliness Issues

### ❌ **Major Issue: Jupyter Notebook Checkpoint Files**

PR #8 includes numerous `.ipynb_checkpoints` directories and files that should **NOT** be in version control:

1. `.ipynb_checkpoints/CMakeLists-checkpoint.txt`
2. `src/.ipynb_checkpoints/field-checkpoint.hpp`
3. `src/.ipynb_checkpoints/gridlayout-checkpoint.hpp`
4. `src/.ipynb_checkpoints/hybirt-checkpoint.cpp`
5. `src/.ipynb_checkpoints/particle-checkpoint.hpp`
6. `src/.ipynb_checkpoints/pusher-checkpoint.hpp`
7. `src/.ipynb_checkpoints/utils-checkpoint.hpp`
8. `src/.ipynb_checkpoints/vecfield-checkpoint.hpp`
9. `tests/boris/.ipynb_checkpoints/CMakeLists-checkpoint.txt`
10. `tests/boris/.ipynb_checkpoints/test_boris-checkpoint.cpp`

**Total unwanted files**: 10 files across 12 changed files in the PR

### Recommendations:
- Add `.ipynb_checkpoints/` to `.gitignore`
- Remove all checkpoint files from the PR
- Only keep the actual modified source file: `src/pusher.hpp`

---

## 4. Test Infrastructure Review

### Tests CMakeLists.txt Analysis

#### ✅ `tests/boris/CMakeLists.txt`
- **add_test() present**: YES
- **Command**: `add_test(NAME test-boris COMMAND test-boris)`
- **Status**: Properly configured

#### ❌ `tests/waves/CMakeLists.txt`
- **add_test() present**: NO
- **Issue**: Missing `add_test()` command
- **Fix needed**: Add `add_test(NAME test-waves COMMAND test-waves)`

#### ❌ `tests/ion_beam/CMakeLists.txt`
- **add_test() present**: NO
- **Issue**: Missing `add_test()` command
- **Fix needed**: Add `add_test(NAME test-ion_beam COMMAND test-ion_beam)`

### Root CMakeLists.txt Analysis
- **tests/boris**: ✅ Included (line 77)
- **tests/waves**: ❌ Not included
- **tests/ion_beam**: ❌ Not included

---

## 5. Summary and Recommendations

### Critical Issues in PR #8
1. **Boris pusher implementation is INCORRECT** - Wrong magnetic rotation formula
2. **Multiple TODO blocks not implemented** (Ampere, moments, population, ICN integration)
3. **Poor code cleanliness** - 10 unwanted checkpoint files
4. **Incomplete implementation** - Only one TODO block attempted (and incorrectly)

### Required Actions
1. **Fix Boris pusher** with correct s-vector formula
2. **Remove all `.ipynb_checkpoints` files** from the PR
3. **Add `.ipynb_checkpoints/` to .gitignore**
4. **Consider implementing remaining TODO blocks** or clearly state PR scope
5. **Add missing add_test() commands** to waves and ion_beam test CMakeLists
6. **Add missing test subdirectories** to root CMakeLists.txt

### Overall Assessment
**Quality: Poor | Correctness: Failed | Cleanliness: Poor**

The PR requires significant rework before it can be merged. The single implemented feature (Boris pusher) contains a critical algorithmic error that would produce incorrect physics results.
