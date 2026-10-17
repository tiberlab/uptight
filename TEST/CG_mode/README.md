# Coarse-Graining Test Case (`CG_mode`)

Standalone driver that compares **full** diagonalization with the two
supported coarse-graining partition strategies in the current codebase.

## Modes

| Mode | Label | `partition_type` | Key inputs |
|------|--------|------------------|------------|
| 1 | STANDARD | — | none (full H) |
| 2 | COUPLING | 0 | `num_blocks` (METIS/graph) |
| 3 | NESTED | 2 | `nested_axes`, `nested_grids`, `layer_tolerance` |

Flow for modes 2/3 (same as production path in `uptight.F90`):

```
UPT_configure_coarse_graining
  → UPT_set_coarse_graining_partition
  → [nested only] UPT_set_coarse_graining_nested
  → [nested only] UPT_set_coarse_graining_nested_tiers
  → upt_hamiltonian   # sparse_ham → icgn_prepare → swap_in
  → UPT_lapack / UPT_jd / UPT_lanczos   # solve + icgn_finalize_eigenvectors
```

Solvers only see the operator currently in `upt%ham` (original or reduced).

## Files

| File | Role |
|------|------|
| `input` | Parameters (format below) |
| `superlattice.upg` | Structure |
| `GaAs.etb` / `AlAs.etb` | Material parameters |
| `../../src/lib_uptight/test_CG_mode.f90` | Driver source |

## `input` format (one value per line)

```
 1: structure file (.upg)
 2: relativistic (.true./.false.)
 3: Harrison scaling (.true./.false.)
 4: c-axis (3 floats)
 5: solver (LK / JD / LO)
 6: nVB
 7: nCB
 8: lambda_vb (eV)
 9: lambda_cb (eV)
10: partition_type   (0 = coupling, 2 = nested; both modes are always run)
11: num_blocks       (coupling only; 0 = auto; ignored for nested)
12: metis_imbalance
13: core_emin (eV)
14: core_emax (eV)
15: top_buffer (eV)
16: bottom_buffer (eV)
17: outer_coupling_threshold (epsilon)
18: add_core_outers  (0 or 1)
19: neumann_order    (<0 = off)
20: expansion_energy E0 (0.0 = auto = core midpoint)
21: nested_axes  (3 ints in {1,2,3}, permutation)
22: nested_grids (3 positive ints)
23: layer_tolerance (Angstrom; <=0 = auto)
24: n_up   (AAD vs standard: n smallest positive evals)
25: n_down (AAD vs standard: n largest negative evals)
```

## Build & run

From `src/lib_uptight/`:

```bash
make uptight.a
# link/compile the driver against the uptight objects (project Makefile target)
make test_CG_mode   # or equivalent local target
```

From `TEST/CG_mode/`:

```bash
../../src/lib_uptight/test_CG_mode
```

Outputs:

- `eigenvalues_standard.dat`
- `eigenvalues_coupling.dat`
- `eigenvalues_nested.dat`

Each file lists sorted eigenvalues plus AAD vs the standard reference when
`n_up` / `n_down` > 0.

## Sanity checks

1. **Coupling with `num_blocks = 1`** and a wide core window covering the full
   spectrum → eigenvalues should match STANDARD (within solver tolerance).
2. **Nested with `nested_grids = (1,1,1)`** and the same wide window → same.
3. Narrower windows reduce dimension; AAD on near-gap states should stay small
   if the window still covers the bands of interest.
