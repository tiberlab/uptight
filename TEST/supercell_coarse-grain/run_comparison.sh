#!/bin/bash
# Script to compare standard vs coarse-grained calculations

TEST_SUPERCELL="../../src/lib_uptight/test_supercell"

if [ ! -f "$TEST_SUPERCELL" ]; then
    echo "ERROR: test_supercell not found. Please compile first:"
    echo "  cd ../../src/lib_uptight"
    echo "  make test_supercell"
    exit 1
fi

if [ ! -f "alloy.upg" ]; then
    echo "ERROR: alloy.upg not found in current directory"
    echo "Please provide structure file"
    exit 1
fi

echo "========================================"
echo "Coarse-Graining Performance Comparison"
echo "========================================"
echo ""

# Backup original config
cp config config.backup

# Test 1: Standard calculation (n_blocks = 1)
echo "Test 1: Standard calculation (no coarse-graining)"
echo "----------------------------------------"
sed -i '10s/.*/1/' config  # Set n_blocks = 1

$TEST_SUPERCELL > output_standard.log 2>&1

if [ -f "eigenvalues.dat" ]; then
    mv eigenvalues.dat eigenvalues_standard.dat
    NBANDS_STD=$(grep "Total number of bands" eigenvalues_standard.dat | awk '{print $6}')
    TIME_STD=$(grep "Solve time" eigenvalues_standard.dat | awk '{print $5}')
    echo "  Bands: $NBANDS_STD"
    echo "  Time:  ${TIME_STD}s"
else
    echo "  ERROR: Standard calculation failed"
    cat output_standard.log
    exit 1
fi
echo ""

# Test 2-7: Coarse-grained calculations with different block counts
for NBLOCKS in 2 4 8 16 32 64; do
    echo "Test: Coarse-graining with $NBLOCKS blocks"
    echo "----------------------------------------"
    
    # Set n_blocks in config (line 10)
    sed -i "10s/.*/$NBLOCKS/" config
    
    $TEST_SUPERCELL > output_cg_${NBLOCKS}.log 2>&1
    
    if [ -f "eigenvalues.dat" ]; then
        mv eigenvalues.dat eigenvalues_cg_${NBLOCKS}.dat
        
        NBANDS=$(grep "Total number of bands" eigenvalues_cg_${NBLOCKS}.dat | awk '{print $6}')
        TIME=$(grep "Solve time" eigenvalues_cg_${NBLOCKS}.dat | awk '{print $5}')
        REDUCTION=$(grep "Rank reduction" eigenvalues_cg_${NBLOCKS}.dat | awk '{print $4}')
        
        # Calculate speedup
        SPEEDUP=$(echo "scale=2; $TIME_STD / $TIME" | bc)
        
        echo "  Bands:     $NBANDS"
        echo "  Time:      ${TIME}s"
        echo "  Speedup:   ${SPEEDUP}x"
        echo "  Reduction: $REDUCTION"
        
        # Compare first 10 eigenvalues if both have them
        echo "  Comparing eigenvalues..."
        python3 - << EOF 2>/dev/null || echo "  (Python not available for comparison)"
import sys
try:
    # Read eigenvalues
    std = []
    with open('eigenvalues_standard.dat', 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                std.append(float(parts[1]))
    
    cg = []
    with open('eigenvalues_cg_${NBLOCKS}.dat', 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                cg.append(float(parts[1]))
    
    # Compare first 10 or min available
    n_compare = min(10, len(std), len(cg))
    max_error = 0.0
    for i in range(n_compare):
        error = abs(std[i] - cg[i])
        if error > max_error:
            max_error = error
    
    print(f"  Max error (first {n_compare} bands): {max_error:.6f} eV")
    
    if max_error < 0.01:
        print("  ✓ Excellent agreement (< 0.01 eV)")
    elif max_error < 0.1:
        print("  ✓ Good agreement (< 0.1 eV)")
    else:
        print("  ⚠ Warning: Large errors")
        
except Exception as e:
    print(f"  Error comparing: {e}", file=sys.stderr)
EOF
        
    else
        echo "  ERROR: Calculation failed"
        cat output_cg_${NBLOCKS}.log | tail -20
    fi
    echo ""
done

# Restore original config
mv config.backup config

echo "========================================"
echo "Summary written to:"
echo "  eigenvalues_standard.dat"
echo "  eigenvalues_cg_*.dat"
echo "  output_*.log"
echo "========================================"
echo ""
echo "To plot results, use provided plot script or:"
echo "  gnuplot plot_comparison.gp"
