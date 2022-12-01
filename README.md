# WHAM Code!

## Usage

This code has a basic set of command line arguments

  -h, --help        show this help message and exit

  -Nw NW            Number of Windows

  -rlow RLOW        Low cutoff

  -rhigh RHIGH      High cutoff

  -nbin NBIN        Histogram bins

  -k K              Force constant

  -plot PLOT        If 1, plots histograms in matplotlib

  -subfile SUBFILE  File name for colvar in subdirectory

  -unit UNIT        [0] Angstrom [1] bohr (output always angstrom, kcal/mol)

Note that rlow and rhigh should be selected in your given unis (i.e. if unit=bohr, then make sure rlow and rhigh are selected in bohr). Everything gets converted to angstrom based units in the end. I.e. if you select bohr, you should give your force constant in hartree/bohr^2 and distance files in bohr, and then it will get coverted to kcal/mol/angstrom^2 and angstrom in the end when the pmf is calculated.

Note: the code presently uses a harmonic function as the bias function, of the form U=0.5*k*x^2 (note the inclusion of the 0.5 factor).

## Required Files:

1) Subdirectories (0-Nw-1) each with a biased file that contains just the distances of the collective variable.

2) "wham_metadata.info" A file that can contain whatever, as long as the second column is the location of the biasing potential for each window (in angstroms).

## Notes

Currently does only 1D problems, and assumes that the collective variable is not periodic. 




