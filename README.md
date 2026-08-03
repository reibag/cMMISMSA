An implementation of Classical Molecular Mechanics Implicit Solvation Model Surface Area (cMMISMSA)
---------------------------------------------------------------------------------------------------

## Based on the following papers:

* A general screened Coulomb potential based implicit solvent model: Calculation of secondary structure of small peptides (https://doi.org/10.1002/qua.1210)

* MM-ISMSA: An Ultrafast and Accurate Scoring Function for Protein–Protein Docking (https://doi.org/10.1021/ct300497z)


<pre>
    Rev 5 - February 2021
      - Fixed bugs in PDB atom typing and protein typing based on dictionary
      - Fixed bugs in PDB procesing (multichain systems)
      - Automatic charges for protein (AMBER) and ligand (PEOE) in PDBs
      - Unitary tests

    Rev 4 - March 2018
      - Added full vdW/qq/solv residue decomposition in output for COMBINE-like methods
      - Added support for XTC GROMACS trajectories (still require AMBER topology)
      - Excluded atoms mask selects individual atoms and not full residues to allow
        better analysis of protein interfaces and Ala-scanning
      - LAPACK can be disabled on copilation time to enable a more portable code

    Rev 3 - January 2016
      - MultiPDB support (experimental) for GROMACS exports

    Rev 2 - November 2015
      - Finished QH entropy approximation. Fixed bugs.
      - Multiple files trajectories (experimental, limited naming scheme)
      - Exclusion regions for protein regions analysis
      - Fixed bugs regarding residue decomposition

    Rev 1 - January-April 2014
     C reimplementation with fancy stuff:
      - No dry topology needed. Waters/ions stripped on-the-fly
      - "Advanced" atom masks to select ligand
      - Automatic "box info" detection
      - PDB with residue contributions in b-factor columns
      - Much better memory management (no more SegFaults on < 2Gb boxes).
        But OVERLAPS still limited to 300
      - Corrected Hydrogen bond model (no more ionic problems). 3 layers
      - Simpler output (parseable CSV files)
      - Debug info control
      - Quasiharmonic entropy approximation thanks to the indexation of
        the trajectory
      - Faster/Slower depending on the case than Fortran MMISMSA :/

</pre>

## Build instructions (main branch compatibility + entropy)

The project now supports entropy-enabled builds directly from the Makefile in `src/`.

### 1) Build with entropy (default)

From `src/`:

```bash
make clean
make
```

This default build enables entropy (`-D_ENTROPY`) and links against LAPACKE/OpenBLAS.

### 2) Build without entropy

From `src/`:

```bash
make noentropy
```

Use this mode when LAPACKE/OpenBLAS are not available.

### 3) Linux example

Install required development packages (package names may vary by distro), then build:

```bash
cd src
make clean
make
```

If headers/libraries are installed in a non-standard prefix:

```bash
make LAPACK_PREFIX=/opt/lapack
```

### 4) MSYS2 / MinGW64 example

Install dependencies in MSYS2 and build inside a MinGW64 shell:

```bash
pacman -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-openblas mingw-w64-x86_64-lapack
cd src
make clean
make
```

The Makefile also accepts an explicit prefix, e.g.:

```bash
make LAPACK_PREFIX=/mingw64
```

### 5) Windows + make from PowerShell/cmd

If you invoke `make` from PowerShell/cmd (instead of MSYS2 shell), pass the MinGW64 prefix explicitly when needed:

```powershell
cd src
make clean
make LAPACK_PREFIX=C:/msys64/mingw64
```

### Notes

* Build target: `MM`
* Entropy target: `make entropy`
* Debug target: `make debug`
* Clean target is portable for both Unix-like shells and Windows cmd
