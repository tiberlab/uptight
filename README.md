uptight
=======

Uptight is a code for empirical tight binding (ETB) calculations. Its original
implementation dates back to 2002 by Jerome Gleize. Subsequently, it has
been further developed in the Opto- and Nanoelectronics group at the
Department of Electronic Engineering of Tor Vergata University of Rome,
under the guidance of Prof. Aldo Di Carlo. Main contributions are due to
Alessandro Pecchia, and Anh-Luan Phan.

The ETB implementation is based on the following papers:
* Slater J C, and Koster G F: "Simplified LCAO Method for the Periodic
  Potential Problem”, Phys. Rev. 94, no. 6, 1498–1524, 1954.
* Jancu J-M, Scholz R, Beltram F and Bassani F: "Empirical spds*
  tight-binding calculation for cubic semiconductors: general method
  and material parameters", Phys. Rev. B 57, 6493, 1998
* Tan Y, Povolotskyi M, Kubis T, Boykin T B and Klimeck G: "Transferable
  tight-binding model for strained group IV and III-V materials and
  heterostructures", Phys. Rev. B 94, 045311, 2016

Uptight has been interfaced to TiberCAD multiscale device simulation software
and used mainly through the latter. Several paramterizations are available for
III-V and III-nitride materials, Si, Ge, and in limited amount for 2D materials.

If you use uptight in your work, you may consider citing the following work,
in addition to one of the mentioned papers before, where appropriate:

* Anh-Luan Phan et al.: "Empirical tight-binding method for large-supercell
  simulations of disordered semiconductor alloys", Phys. Scr. 99, 075903, 2024

License
-------

Uptight is distributed as is, under the conditions of the Lesser GNU General Public
License (LGPL) version 3.0.


Installation
------------

### Requirements

- **Operating System:** The following instructions were tested on GNU/Linux (Debian-based distros). MacOS has not been tested but should work similarly. Windows are not officially supported, however using Windows Subsytem for Linux (WSL) may help.
- **Fortran Compiler:** GNU Fortran (`gfortran`) or Intel Fortran (`ifort`)
- **Build Tools:** `make`
- **Mathematical Libraries:**
  - BLAS and LAPACK (standard linear algebra libraries)
  - Optionally, ARPACK (for eigenvalue problems)
- **Optional:** MPI libraries (OpenMPI or MPICH) for parallel computations

### Installing dependencies

On Ubuntu / Debian:

```bash
sudo apt update
sudo apt install gfortran make libblas-dev liblapack-dev libarpack-dev libopenmpi-dev openmpi-bin
```

MPI packages are optional unless you want parallel support.

### Building Uptight

Unlike typical Autotools projects, Uptight **does not** use a root-level `make` or `make install`. Instead, the main build happens inside the `src/lib_uptight/` directory.

**Step 0: Download the source code**

First, clone the Uptight repository from your terminal:

```bash
git clone https://github.com/tiberlab/uptight.git
```

Then, change into the project root directory:

```bash
cd uptight
```

**Step 1: Configure**

From the root directory, run:

```bash
./src/configure
```

This generates necessary configuration files.

**Step 2: Build the core library**

Change into the library directory:

```bash
cd src/lib_uptight
```

To compile **only** the core static library `uptight.a` (recommended):

```bash
make clean
make uptight.a
```

This creates `src/lib_uptight/uptight.a`, which contains all core Uptight routines.

**Step 3 (optional): Build example tool**

To also build the example executable `PARAMETERIZER`:

```bash
make PARAMETERIZER
```

This links against `uptight.a` and produces `src/lib_uptight/PARAMETERIZER`.

### Notes on usage and portability

- The static library `uptight.a` can be linked into your own Fortran programs.
- The `PARAMETERIZER` executable depends on system libraries (LAPACK, BLAS, Fortran runtime).
- You can copy `PARAMETERIZER` to another Linux machine, but it **requires** compatible versions of these libraries installed there.
- The executable is **not** fully standalone; missing libraries will cause it to fail.
- To run your own programs using Uptight, link with:

```bash
gfortran -o my_program my_program.f90 src/lib_uptight/uptight.a -llapack -lblas
```

### Summary

- **Install dependencies first.**
- **Run `./src/configure` from the root.**
- **Build the library with `make uptight.a` inside `src/lib_uptight/`.**
- **Optionally build `PARAMETERIZER` with `make PARAMETERIZER`.**

Usage
-----

Uptight primarily provides a Fortran static library (`uptight.a`) for empirical tight binding calculations.

You can:

- Link `uptight.a` into your own Fortran programs.
- Use the example executable `PARAMETERIZER` (if built) to perform parameterization tasks.
- Integrate Uptight with TiberCAD multiscale simulation software.

Refer to the example input files and test cases in the `TEST/` directory for guidance on preparing inputs and running calculations.
