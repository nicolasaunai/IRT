# PR #18 Review Report

## Overview
This report reviews the code implementations in PR #18 "Pablo De Frutos Rull: final pull request with detailed notebooks" against the ground truth implementation from https://github.com/nicolasaunai/hybirt.

## Code Implementation Reviews

### 1. Ampere Implementation (`src/ampere.hpp`)

**Status**: ✅ **CORRECT**

**Analysis**:
The implementation in PR #18 correctly fills the TODO block for computing the current density J from the magnetic field B using Ampere's law.

**Comparison with Ground Truth**:
- Both implementations compute `Jy` and `Jz` as primal quantities in the x-direction
- The formulas are identical:
  - `Jy(ix) = -(Bz(ix) - Bz(ix-1))/dx`
  - `Jz(ix) = (By(ix) - By(ix-1))/dx`
- PR #18 includes an extra loop to set `Jx` to zero (dual quantity), which is correct but not present in ground truth

**Quality Assessment**:
- **Good**: The implementation is mathematically correct
- **Good**: Proper handling of grid quantities (primal/dual)
- **Minor difference**: PR #18 explicitly sets Jx=0, while ground truth doesn't (both are acceptable)
- **Code clarity**: Well-commented and structured

### 2. Faraday Implementation (`src/faraday.hpp`)

**Status**: ⚠️ **MOSTLY CORRECT WITH DIFFERENCES**

**Analysis**:
The implementation in PR #18 correctly implements Faraday's law to update the magnetic field from the electric field.

**Comparison with Ground Truth**:
- **Architecture difference**: Ground truth stores `dt` as a member variable and takes it in constructor, PR #18 passes `dt` as operator() parameter
- **Loop structure difference**: 
  - Ground truth uses `ghost_start/ghost_end` for all quantities
  - PR #18 uses `primal_dom_start/primal_dom_end` for Bx and `dual_dom_start/dual_dom_end` for By/Bz
- The core physics formulas are identical:
  - `Bnew_y(ix) = By(ix) + dt * (Ez(ix+1) - Ez(ix))/dx`
  - `Bnew_z(ix) = Bz(ix) - dt * (Ey(ix+1) - Ey(ix))/dx`
  - `Bnew_x(ix) = Bx(ix)` (constant)

**Quality Assessment**:
- **Good**: Mathematically correct formulas
- **Design choice**: Different API (dt as parameter vs member) - PR #18's approach is more flexible
- **Potential issue**: Loop bounds differ - ground truth uses ghost cells, PR #18 uses domain cells only
- **Recommendation**: The loop bounds difference could affect boundary handling

### 3. Boris Pusher Implementation (`src/pusher.hpp`)

**Status**: ⚠️ **FUNCTIONALLY CORRECT BUT SIMPLIFIED**

**Analysis**:
The implementation in PR #18 implements the Boris particle pusher algorithm.

**Comparison with Ground Truth**:
- **Major simplification**: PR #18 uses separate variables for intermediate steps (v_minus, vector_t, v_prime, vector_s, v_plus), while ground truth uses more compact notation
- **Position update difference**: 
  - Ground truth: Has safety checks for CFL condition (particle shouldn't move more than half a cell)
  - PR #18: No CFL safety checks
- **Variable naming**: PR #18 uses more descriptive names (e.g., `vector_t`, `vector_s`) vs ground truth compact style
- **Algorithm correctness**: The core Boris algorithm is correctly implemented in both

**Quality Assessment**:
- **Good**: Algorithm is mathematically correct
- **Missing feature**: No CFL violation checks (present in ground truth)
- **Code clarity**: More verbose but perhaps more readable
- **Potential issue**: Interpolation logic appears identical and correct

### 4. Population Deposit Implementation (`src/population.hpp`)

**Status**: ⚠️ **CORRECT WITH EXTRA CODE**

**Analysis**:
The implementation correctly performs linear weighting deposit for density and flux.

**Comparison with Ground Truth**:
- **Core algorithm**: Identical particle-to-grid deposit using linear (cloud-in-cell) weighting
- **Extra code in PR #18**: Contains a commented-out periodic boundary handling block that is unnecessary
- **Formula differences**: 
  - Ground truth: `m_flux.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0]`
  - PR #18: `m_flux.x(iLeft) += particle.v[0] * particle.weight * (1.0 - reminder)`
  - These are mathematically equivalent (commutativity)

**Quality Assessment**:
- **Good**: Correct implementation of linear weighting
- **Good**: Properly handles both density and flux deposition
- **Minor issue**: Contains commented-out unnecessary code
- **Good**: Self-aware comment acknowledging boundary conditions are handled elsewhere

### 5. Moments Implementation (`src/moments.hpp`)

**Status**: ❌ **INCORRECT**

**Analysis**:
The implementation has a critical bug in the bulk_velocity calculation.

**Comparison with Ground Truth**:
- **total_density**: Identical and correct
- **bulk_velocity ISSUE**: 
  - Ground truth divides by `pop.density()(ix)` for each population
  - PR #18 divides by total `N(ix)`
  - **Impact**: This is mathematically different and incorrect

**Detailed Issue**:
```cpp
// Ground truth (correct):
V.x(ix) /= pop.density()(ix);  // Divides by individual population density

// PR #18 (incorrect):
V.x(ix) = V.x(ix)/N(ix);  // Divides by total density N
```

**Quality Assessment**:
- **Critical bug**: The bulk velocity calculation is incorrect
- **Impact**: Results will be wrong for multi-population simulations
- **Note**: Ground truth itself has unusual logic (loops over populations but divides by each pop's density, which seems odd for "bulk" velocity)

### 6. ICN Temporal Integration (`src/hybirt.cpp`)

**Status**: ⚠️ **DIFFERENT IMPLEMENTATION**

**Analysis**:
The PR #18 implements an Implicit-Crank-Nicolson (ICN) temporal integration scheme, but differently from ground truth.

**Comparison with Ground Truth**:
- **Structure**: Both use predictor-corrector approach
- **Key difference in Predictor 1**:
  - Ground truth: Uses `save{populations}` to save particle state before first predictor
  - PR #18: Updates populations directly without saving
- **Boundary condition handling**:
  - Ground truth: More extensive use of boundary_condition->fill() after each step
  - PR #18: Less frequent boundary filling
- **Average function**: Identical implementation

**Quality Assessment**:
- **Functional difference**: The lack of particle state saving in PR #18 could lead to incorrect predictor-corrector behavior
- **Potential issue**: Predictor 2 should restart from time n, but PR #18's implementation continues from predictor 1 state
- **Design**: Ground truth's approach is more correct for a proper predictor-corrector scheme

### 7. Average Function (`src/hybirt.cpp`)

**Status**: ✅ **CORRECT**

**Analysis**:
Simple element-wise averaging using std::transform.

**Comparison**: Identical to ground truth.

## PR Cleanliness Issues

### ❌ Unnecessary Files Present

The PR contains several files that should **NOT** be committed to version control:

1. **Jupyter Notebook Checkpoint Files** (16 files):
   - `.ipynb_checkpoints/CMakeLists-checkpoint.txt`
   - All files in `src/.ipynb_checkpoints/` (15 files)
   - All files in `tests/ampere/.ipynb_checkpoints/` (2 files)
   - All files in `tests/boris/.ipynb_checkpoints/` (2 files)
   - All files in `tests/faraday/.ipynb_checkpoints/` (2 files)
   - All files in `tests/moments/.ipynb_checkpoints/` (2 files)
   - All files in `tests/waves/.ipynb_checkpoints/` (1 file)

2. **Why these shouldn't be committed**:
   - These are auto-generated by Jupyter Notebook
   - They are duplicates of actual source files
   - They clutter the repository
   - They can cause merge conflicts
   - Standard practice is to add `.ipynb_checkpoints/` to `.gitignore`

### ✅ Legitimate Files

The following files are legitimate additions:
- `tests/ampere/CMakeLists.txt` and `test_ampere.cpp`
- `tests/faraday/CMakeLists.txt` and `test_faraday.cpp`
- `tests/moments/CMakeLists.txt` and `test_moments.cpp`
- Modified source files with implemented TODO blocks

## Summary

### Code Quality Summary:
| File | Status | Correctness | Quality |
|------|--------|-------------|---------|
| ampere.hpp | ✅ Correct | Excellent | Minor differences |
| faraday.hpp | ⚠️ Different | Good | Different design choice |
| pusher.hpp | ⚠️ Simplified | Good | Missing safety checks |
| population.hpp | ⚠️ Extra code | Good | Has commented code |
| moments.hpp | ❌ Incorrect | **FAIL** | Critical bug |
| hybirt.cpp | ⚠️ Different | Questionable | Different algorithm |

### Critical Issues:
1. **moments.hpp**: Incorrect bulk_velocity calculation
2. **hybirt.cpp**: Missing particle state save/restore in predictor-corrector
3. **26 unnecessary checkpoint files** should be removed

### Recommendations:
1. Fix the bulk_velocity calculation in moments.hpp
2. Review and potentially fix the ICN integration to match ground truth
3. Add CFL checks in Boris pusher for robustness
4. Remove all `.ipynb_checkpoints/` files
5. Add `.ipynb_checkpoints/` to `.gitignore`
