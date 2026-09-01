# Coarse-grained tight-binding solver layer

Uptight can optionally reduce a Hamiltonian after `sparse_ham` has assembled
it.  Atoms are partitioned with weighted METIS graph partitioning, each
partition is diagonalized in a dense local basis, and only local eigenvectors
inside one energy window are retained.  The final solver receives
`Q^H H Q`; returned eigenvectors are lifted with `Q` into the original
atomic-orbital basis.

The feature is disabled by default.  It retains the original CSR Hamiltonian
for existing inspection and matrix-element APIs, and supports a single MPI
rank (including the single-GPU JD path).  Configure with
`upt_set_coarse_graining`, then build the Hamiltonian normally.  METIS is an
optional configure-time dependency; enabling the feature in a build without
METIS produces an actionable error.

Implementation details and validation are kept next to the implementation in
`src/lib_uptight/coarse_grain.f90`.  The algorithm follows Liu et al.,
"Coarse-grained tight-binding models", J. Phys.: Condens. Matter 34 (2022)
125901.
