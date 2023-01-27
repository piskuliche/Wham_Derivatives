#!/usr/bin/env python
import numpy as np
import argparse
from wham_derivatives import wham_class

def Main(rlow=1.5, rhigh=8.0, nbins=200, deriv=0, enerfile="flucts.inp", metafile="wham_metadata.info", maxiter=1000000):
    """This is an example of how these functions can be called to calculate WHAM

    Args:
        rlow (float): Lower cutoff for the WHAM calculation [default=1.5]
        rhigh (float): Upper cutoff for the WHAM calculation [default=8.0]
        nbins (int): Number of windows for the WHAM calculation [defaul=200]
        deriv (int): [0] Don't calculate derivative [1] Calculate derivative [default=0]
        enerfile (str): Name of file with information about columns and energies [default='flucts.inp']
        metafile (str): Name of file with information about umbrella potential locations [default='wham_metadata.info']

    Raises: 
        OSError: Metadata file was non-existent or in wrong format.

    """
    nwindows=0
    xc, k=[], []
    try:
        xc,k=np.genfromtxt(metafile,usecols=(1,2),unpack=True)
        nwindows = len(xc)
    except:
        exit("Error: Trouble grabbing windows from metafile")
    
    whammed = wham_class.Wham(xc, k[0], rlow, rhigh, nwindows, nbins, eweight=True)
    whammed.Do_WHAM()
    whammed.Plot_PMF()
    whammed.Do_WHAM_D("TotEng", maxiter=maxiter)


if __name__ == "__main__":


    from wham_class import Wham

    parser=argparse.ArgumentParser()
    parser.add_argument('-rlow', default=1.5,type=float,
                         help='Low cutoff in angstroms')
    parser.add_argument('-rhigh', default=8,type=float,
                         help='High cutoff in angstroms')
    parser.add_argument('-nbins', default=100, type=int,
                         help='Histogram bins')
    parser.add_argument('-deriv',default=0, type=int,
                         help="[0] Turn off derivative calculation (default) [1] Turn on derivative calculation")
    parser.add_argument('-enerfile', default="flucts.inp", type=str,
                         help="File name with energy information, default flucts.inp")
    parser.add_argument('-metafile', default='wham_metadata.info', type=str,
                         help='File name with directory, loc, and k info.')
    Iargs = parser.parse_args()

    Main(rlow=Iargs.rlow, rhigh=Iargs.rhigh, nbins=Iargs.nbins, deriv=Iargs.deriv, enerfile=Iargs.enerfile, metafile=Iargs.metafile)