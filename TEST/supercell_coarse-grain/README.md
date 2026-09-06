# Coarse-Graining Test Case

This directory contains test inputs for the coarse-graining feature.

## Files Required

1. **config** - Configuration parameters (see format below)
2. **alloy.upg** - Atomistic structure file (UPtight Geometry format)
3. **material.etb** - Empirical tight-binding parameters for materials used

## Config File Format

The `config` file contains the following parameters (one per line):

```
Line 1:  Structure filename (e.g., alloy.upg)
Line 2:  Relativistic flag (.true. or .false.)
Line 3:  Harrison scaling (.true. or .false.)
Line 4:  C-axis direction (3 real numbers, e.g., 0.0 0.0 1.0)
Line 5:  Solver choice (LK=LAPACK, JD=Jacobi-Davidson, LO=Lanczos)
Line 6:  Number of valence bands (used only if n_blocks=1)
Line 7:  Number of conduction bands (used only if n_blocks=1)
Line 8:  Valence band energy guess (eV)
Line 9:  Conduction band energy guess (eV)
Line 10: Number of blocks (1=no coarse-graining, >1=enable coarse-graining)
Line 11: Energy window minimum Emin (eV)
Line 12: Energy window maximum Emax (eV)
Line 13: METIS imbalance tolerance (e.g., 0.03 for 3%)
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
- Solve at Gamma point (k = 0, 0, 0) only
- Write results to `eigenvalues.dat`

## Modes of Operation

### Standard Mode (n_blocks = 1)

When `n_blocks = 1`, the program runs in standard mode:
- Solves for `nVB` valence bands and `nCB` conduction bands
- Uses standard tight-binding without coarse-graining
- Energy window parameters (Emin, Emax) are ignored

### Coarse-Graining Mode (n_blocks > 1)

When `n_blocks > 1`, coarse-graining is enabled:
- Partitions atoms into `n_blocks` blocks using METIS
- Retains only eigenstates with energies in [Emin, Emax]
- Solves reduced Hamiltonian
- Returns ALL eigenvalues in the energy window
- nVB and nCB parameters are ignored

## Output

The program creates `eigenvalues.dat` with:
- Header with metadata (number of bands, solve time, coarse-graining status)
- Two columns: index and energy (eV)
- Energies sorted in descending order

Example output:
```
# Eigenvalues (eV) - sorted descending
# Total number of bands: 247
# Solve time (s):   12.345678
# Coarse-graining: ENABLED
# Dimension: 5040 -> 247
# Rank reduction:  95.10%
#
# Index    Energy(eV)
     1      12.84523100
     2      12.81204522
     3      12.77831945
     ...
```

## Performance Comparison

To compare standard vs coarse-grained calculations:

1. Run with `n_blocks = 1` (standard)
2. Run with `n_blocks = 16` (coarse-grained)
3. Compare:
   - Solve time
   - Eigenvalues (should match within tolerance)
   - Memory usage

## Recommended Settings

Based on Liu et al. (2022) paper:

| System Size | n_blocks | Emin (eV) | Emax (eV) | Expected Speedup |
|-------------|----------|-----------|-----------|------------------|
| Small (<1000 atoms) | 1 | - | - | Baseline |
| Medium (1000-5000) | 8-16 | -0.5 | 13.0 | 10-20× |
| Large (>5000) | 16-64 | -0.5 | 13.0 | 20-40× |

**Energy window selection**:
- For GaAs/AlAs systems: -0.5 to 13.0 eV captures VB and low CB
- For wider gap materials: adjust accordingly
- Window must cover all bands of interest

**Imbalance tolerance**:
- 0.03 (3%) is recommended
- Lower values: better balance, slower partitioning
- Higher values: faster partitioning, may increase cut edges

## Troubleshooting

**"retained rank is smaller than requested states"**
- Energy window [Emin, Emax] is too narrow
- Solution: Widen the window

**"energy window retained no states"**
- Window doesn't overlap any eigenvalues
- Solution: Check band edge energies and adjust window

**"METIS is unavailable or failed"**
- METIS not installed or not found
- Solution: Install METIS and reconfigure uptight with `--with-metis`
- Or: Use `n_blocks = 1` (no coarse-graining)

**Solver fails with coarse-graining**
- Try different solver (LK for small systems, JD for large)
- Check if reduced dimension is sufficient for solver parameters

## Validation

To validate against the Liu et al. (2022) paper results:

1. Create GaAs/AlAs quantum well structure (12 unit cells)
2. Run with different block counts: 1, 2, 4, 8, 16, 32, 64
3. Energy window: -0.5 to 13.0 eV
4. Compare:
   - Eigenvalue errors (should be < 1% for 16-32 blocks)
   - Speedup (should reach 30-40× at optimal block count)
   - Rank reduction (should be > 60%)
