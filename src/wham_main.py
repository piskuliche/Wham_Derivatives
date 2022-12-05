#!/usr/bin/env python
import numpy as np
import argparse
from wham_class import Wham

def Main(Iargs):
    nwindows=0
    xc, k=[], []
    try:
        xc,k=np.genfromtxt(Iargs.metafile,usecols=(1,2),unpack=True)
        nwindows = len(xc)
    except:
        exit("Error: Trouble grabbing windows from metafile")
    
    whammed = Wham(xc, k[0], Iargs.rlow, Iargs.rhigh, nwindows, Iargs.nbins, eweight=True)
    whammed.Do_WHAM()
    #whammed.Plot_PMF()
    whammed.Do_WHAM_D("TotEng",maxiter=100000)

    return

if __name__ == "__main__":
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

    Main(Iargs)