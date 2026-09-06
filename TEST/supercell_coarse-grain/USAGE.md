# Test Supercell - Coarse-Graining Test Program

## Quick Start

### 1. Compile

```bash
cd /home/luan/Documents/GitHub/uptight/src/lib_uptight
make test_supercell
```

This creates the executable `test_supercell` in `src/lib_uptight/`.

### 2. Prepare Input Files

In `TEST/supercell_coarse-grain/`, you need:

1. **alloy.upg** - Your structure file (UPtight Geometry format)
2. **material.etb** files - ETB parameters for your materials
3. **config** - Configuration file (provided as template)

### 3. Edit Config File

Edit `config` with your parameters. Default template:

```
alloy.upg          # Line 1:  Structure file
.true.             # Line 2:  Relativistic (.true./.false.)
.true.             # Line 3:  Harrison scaling (.true./.false.)
0.0 0.0 1.0        # Line 4:  C-axis direction
LK                 # Line 5:  Solver (LK/JD/LO)
10                 # Line 6:  Number of VB (if n_blocks=1)
10                 # Line 7:  Number of CB (if n_blocks=1)
-0.5               # Line 8:  VB energy guess (eV)
2.0                # Line 9:  CB energy guess (eV)
1                  # Line 10: Number of blocks (1=standard, >1=coarse-grain)
-0.5               # Line 11: Emin (eV)
13.0               # Line 12: Emax (eV)
0.03               # Line 13: METIS imbalance tolerance
```

### 4. Run

```bash
cd /home/luan/Documents/GitHub/uptight/TEST/supercell_coarse-grain
../../src/lib_uptight/test_supercell
```

Output will be written to `eigenvalues.dat`.

## Understanding the Output

### Console Output

The program prints:
- Configuration summary
- Structure building progress
- Hamiltonian construction
- Coarse-graining statistics (if enabled)
- Solve time and number of bands

### eigenvalues.dat File

Format:
```
# Eigenvalues (eV) - sorted descending
# Total number of bands: 247
# Solve time (s):   1.234567
# Coarse-graining: ENABLED
# Dimension: 5040 -> 247
# Rank reduction:  95.10%
#
# Index    Energy(eV)
     1      12.84523100
     2      12.81204522
     ...
```

## Modes

### Standard Mode (n_blocks = 1)

```bash
# In config file, set line 10 to 1
1
```

- Uses standard tight-binding (no coarse-graining)
- Solves for nVB valence + nCB conduction bands
- Returns exactly nVB + nCB eigenvalues
- Emin/Emax are ignored

### Coarse-Graining Mode (n_blocks > 1)

```bash
# In config file, set line 10 to desired blocks (e.g., 16)
16
```

- Enables coarse-graining with METIS partitioning
- Solves for ALL eigenvalues in energy window [Emin, Emax]
- nVB and nCB are ignored
- Returns all retained eigenvalues
- Requires METIS installed

## Performance Comparison

Use the provided script:

```bash
cd /home/luan/Documents/GitHub/uptight/TEST/supercell_coarse-grain
./run_comparison.sh
```

This automatically:
1. Runs standard calculation (n_blocks=1)
2. Runs coarse-grained calculations with 2, 4, 8, 16, 32, 64 blocks
3. Compares eigenvalues
4. Reports speedup and accuracy
5. Creates output files: `eigenvalues_standard.dat`, `eigenvalues_cg_*.dat`

## Choosing Parameters

### Number of Blocks

| System Size (atoms) | Recommended n_blocks | Expected Speedup |
|---------------------|---------------------|------------------|
| < 500              | 1 (standard)        | Baseline         |
| 500 - 2000         | 4 - 8              | 5-10×            |
| 2000 - 5000        | 8 - 16             | 10-20×           |
| 5000 - 10000       | 16 - 32            | 20-40×           |
| > 10000            | 32 - 64            | 30-50×           |

**Note**: Optimal block count depends on system. Too few blocks = limited speedup. Too many blocks = overhead increases.

### Energy Window [Emin, Emax]

**For GaAs/AlAs systems**:
- Conservative: -0.5 to 13.0 eV (captures VB + low CB)
- Narrow: -0.5 to 9.0 eV (faster, may miss some CB states)
- Wide: -4.0 to 16.0 eV (slower, includes more states)

**For other materials**:
- Set Emin below valence band maximum
- Set Emax above highest conduction band of interest
- Window must cover ALL bands you want to compute

**Rule of thumb**: Window should extend ~1 eV beyond bands of interest.

### Solver Choice

| Solver | Best For | Notes |
|--------|----------|-------|
| LK (LAPACK) | Small systems (< 2000 orbitals) | Dense solver, exact but slow for large systems |
| JD (Jacobi-Davidson) | Large systems | Iterative, GPU support, good for sparse systems |
| LO (Lanczos) | Medium to large | Iterative, good convergence |

**With coarse-graining**: Reduced system is smaller, so LK becomes viable for larger original systems.

## Validation

To validate implementation against Liu et al. (2022):

1. **Create test structure**: GaAs/AlAs quantum well (12 unit cells well, 512 total unit cells)

2. **Run baseline**:
   ```bash
   # config line 10: n_blocks = 1
   ```

3. **Run coarse-grained**:
   ```bash
   # Test blocks: 1, 2, 4, 8, 16, 32, 64, 128, 256
   # Energy window: -0.5 to 13.0 eV
   ```

4. **Expected results** (from paper):
   - **Accuracy**: < 0.1% error for 16-64 blocks
   - **Speedup**: 30-40× at 8-64 blocks
   - **Rank reduction**: > 60%

5. **Plot results**:
   - Eigenvalue error vs number of blocks (Figure 3 in paper)
   - Speedup vs number of blocks (Figure 4 in paper)

## Troubleshooting

### Compilation Errors

```bash
# If "module not found" errors:
cd /home/luan/Documents/GitHub/uptight/src/lib_uptight
make clean
make

# Check dependency order in Makefile
```

### "METIS is unavailable or failed"

Two options:

**Option 1**: Install METIS (for full functionality)
```bash
# See main COARSE_GRAINING_GUIDE.md for METIS installation
```

**Option 2**: Use standard mode
```bash
# In config, set n_blocks = 1
```

### "retained rank is smaller than requested states"

Energy window [Emin, Emax] doesn't retain enough states.

**Solution**: Widen the window
```bash
# In config, increase Emax or decrease Emin
# Example: change from (-0.5, 9.0) to (-0.5, 13.0)
```

### "energy window retained no states"

Energy window doesn't overlap any eigenvalues.

**Solution**: Adjust window to cover band edges
```bash
# Check your material's band gap
# For GaAs: VB ~0 eV, CB ~1.5 eV
# Set Emin < 0, Emax > 2.0
```

### Solve Fails or Takes Too Long

**Solutions**:
1. Use faster solver (JD instead of LK)
2. Reduce energy window
3. Use more blocks (smaller reduced dimension)
4. Check verbose output for detailed error messages

### Results Don't Match Expected

**Checklist**:
1. Verify structure file is correct
2. Check ETB parameters loaded correctly
3. Ensure energy window covers bands of interest
4. Try exactness test: n_blocks=1 with full window should match standard
5. Check console output for warnings

## Example Workflow

Complete example for GaAs/AlAs quantum well:

```bash
# 1. Compile
cd /home/luan/Documents/GitHub/uptight/src/lib_uptight
make test_supercell

# 2. Setup test
cd ../../TEST/supercell_coarse-grain
# Copy your alloy.upg and .etb files here

# 3. Edit config for standard run
# Set line 10 to: 1
vi config

# 4. Run standard
../../src/lib_uptight/test_supercell > log_standard.txt
mv eigenvalues.dat eigenvalues_standard.dat

# 5. Edit config for coarse-graining
# Set line 10 to: 16
vi config

# 6. Run coarse-grained
../../src/lib_uptight/test_supercell > log_cg.txt
mv eigenvalues.dat eigenvalues_cg.dat

# 7. Compare results
diff <(head -20 eigenvalues_standard.dat) <(head -20 eigenvalues_cg.dat)

# 8. Extract timing
grep "Solve time" eigenvalues_standard.dat
grep "Solve time" eigenvalues_cg.dat

# 9. Calculate speedup
python3 -c "
std = 12.345  # Replace with actual time from step 8
cg = 0.678    # Replace with actual time from step 8
print(f'Speedup: {std/cg:.1f}x')
"
```

## Advanced Usage

### Custom Energy Windows

Test different windows to find optimal accuracy/speed tradeoff:

```bash
# Tight window (fast)
Emin = 0.0, Emax = 3.0  # Only near band edges

# Medium window (recommended)
Emin = -0.5, Emax = 13.0  # Standard for GaAs/AlAs

# Wide window (accurate)
Emin = -4.0, Emax = 16.0  # Captures more physics
```

### Block Size Sweep

Systematically test different block counts:

```bash
for NBLOCKS in 1 2 4 8 16 32 64; do
    # Update config line 10
    sed -i "10s/.*/$NBLOCKS/" config
    ../../src/lib_uptight/test_supercell > log_${NBLOCKS}.txt
    mv eigenvalues.dat eigenvalues_${NBLOCKS}.dat
done
```

### Batch Processing

Process multiple structures:

```bash
for STRUCT in structure1.upg structure2.upg structure3.upg; do
    # Update config line 1
    sed -i "1s/.*/$STRUCT/" config
    ../../src/lib_uptight/test_supercell > log_${STRUCT}.txt
done
```

## Support

For issues or questions:
1. Check README.md in this directory
2. See main COARSE_GRAINING_GUIDE.md in project root
3. Review Liu et al. (2022) paper methodology
4. Check uptight documentation

## References

- Liu, T.-X., Mao, L., Pistol, M.-E., & Pryor, C. (2022). Coarse-grained tight-binding models. *J. Phys.: Condens. Matter*, 34, 125901.
- METIS: https://github.com/KarypisLab/METIS
