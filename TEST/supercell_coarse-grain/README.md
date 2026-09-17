# Coarse-Graining Test Case

This directory contains test inputs for the coarse-graining feature.

## Files Required

1. **config** - Configuration parameters (see format below)
2. **alloy.upg** - Atomistic structure file (UPtight Geometry format)
3. **material.etb** - Empirical tight-binding parameters for materials used

## Config File Format

The `config` file contains the following parameters (one per line):

```
Config file format (one value per line):
  1:  structure file (.upg)
  2:  relativistic (.true./.false.)
  3:  Harrison scaling (.true./.false.)
  4:  c-axis (3 floats)
  5:  solver (LK / JD / LO)
  6:  nVB  (standard mode)
  7:  nCB  (standard mode)
  8:  lambda_vb  (eV)
  9:  lambda_cb  (eV)
  10: n_blocks   (for CG & ICG & ICGN)
  11: cg_emin    (eV)
  12: cg_emax    (eV)
  13: imbalance  (METIS)
  14: icg_core_emin  (eV)
  15: icg_core_emax  (eV)
  16: icg_top_buffer   (eV)
  17: icg_bottom_buffer (eV)
  18: icg_epsilon    (threshold factor)
  19: sub_tolerance (CG block solver tolerance)
  20: icgn_selfenergy_order (0,1,2,...)
  21: icgn_E0     (eV, 0.0 = auto = core window midpoint)
  22: n_up   (n smallest positive eigenvalues for AAD, 0 to skip)
  23: n_down (n largest  negative eigenvalues for AAD, 0 to skip)
  24: check_neumann_convergence (.true./.false.)
  25: power_iteration_max_iterations
  26: power_iteration_tolerance
```

## Running the Test

### Compile

From `src/lib_uptight/`:
```bash
make test_supercell
```

### Run

From `TEST/supercell_coarse-grain/`:
```bash
../../src/lib_uptight/test_supercell
```

The program will:
- Read configuration from `config`
- Build the structure and Hamiltonian
- Solve at Gamma point (k = 0, 0, 0) only using individual modes in the order: standard full diagonalization (standard), coarse-graining (cg), improved coarse-graining (icg), improved coarse-graining + Neumann corrections for self-energy (icgn)
- Write results to `eigenvalues_<mode-label>.dat` where `<mode-label>` is either "standard", "cg", "icg", "icgn"

## Santity checks:
1. Setting `n_blocks` to 1 in line 10 should make results of other modes exactly the same as that of the "standard" mode, no matter the other inputs
2. Setting a very large energy window in "cg", "icg" and "icgn" modes that covers the whole spectrum of the Hamiltonian should make the results of these modes exactly the same as that of the "standard" mode, no matter the other inputs, no matter the number of blocks

