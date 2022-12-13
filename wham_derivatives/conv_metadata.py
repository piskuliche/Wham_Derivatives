#!/usr/bin/env python

import os, shutil
import numpy as np

def Conv(Iargs):
    """Converts the Grossfield WHAM-code-like metadata file units

    This converts units for the Grossfield WHAM-code metafile. This metafile is in the format
    of one line per window, each line having filename, bin location, and k.

    Args:
        Iargs (argparse): Input arguments from argparse
    
    Raises:
        ValueError: Incorrect input arguments were provided

    """
    bohr_to_angstrom = 0.529177
    hartree_to_kcal = 627.509
    
    col, x, k = Read_Meta()
    if Iargs.xunits not in ['bohr','angstrom']:
        raise ValueError("Error: xunits must be bohr or angstrom")
    if Iargs.eunits not in ['hartree, kcal']:
        raise ValueError("Error: units must be in hartree or kcal")
    
    if Iargs.override_k != 0.0:
        k = [Iargs.override_k]*len(x)

    if Iargs.xunits == 'bohr':
        x = x*bohr_to_angstrom
    if Iargs.eunits == 'hartree':
        k = k*hartree_to_kcal/(bohr_to_angstrom**2.)
    
    with open("wham_metadata.info",'w') as f:
        for i in range(len(col)):
            line = "%s %10.5f %10.5f\n" % (col[i], x[i], k[i])
            f.write(line)

    return

def Read_Meta():
    try:
        col, x, k = np.genfromtxt("wham_metadata.info",usecols=(0,1,2),unpack=True)
    except:
        raise OSError("Error: Couldn't read wham_metadata.info")

    shutil.copyfile("wham_metadata.info",'wham_metadata.info.old')

    return col, x, k

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-xunits', default='bohr',type=str,
                         help='Distance units type [bohr, angstrom]')
    parser.add_argument('-kunits', default='hartree',type=str,
                         help='Energy units type [hartree, kcal]')
    parser.add_argument('-override_k', default=0.0, type=float,
                         help='Value to override force constant with')
    Iargs = parser.parse_args()

    Conv(Iargs)