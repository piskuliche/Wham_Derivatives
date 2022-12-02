#!/usr/bin/env python
import argparse, os, pickle
import numpy as np

def Main(Iargs):
    """Does the pre-processing.

    This does basic preprocessing of files. All unit conversions, skippage,
    etc needs to be done ahead of time. 

    Args:
        Iargs (argparse): Input arguments

    """

    Write_WHAM_Files(Iargs)

    return

def Write_WHAM_Files(Iargs):
    """Writes files in a very simple, WHAM-readable format. 

    This takes files, and writes them in a very WHAM-readable format, 
    mainly by converting them to pckl files in a centralized directory.

    Args:
        Iargs (argparse): Input arguments
    """
    if not os.path.exists("wham_pckl"):
        os.makedirs('wham_pckl')
    
    nwindows=0
    try:
        k=np.genfromtxt(Iargs.metafile,usecols=1,unpack=True)
        nwindows = len(k)
    except:
        exit("Error: Trouble grabbing windows from metafile")
    
    etypes = []
    if Iargs.deriv == 1:
        try:
            etypes=np.genfromtxt(Iargs.enerfile, usecols=0, dtype=str)
            print("Grabbing energies")
            print("Selected energy types:",*etypes)
        except:
            exit("Error: Trouble grabbing energies")

    svdata = []
    for window in range(nwindows):
        if window%10 == 0: print(window)
        data, en = [], {}
        
        data = np.genfromtxt(str(window)+"/"+Iargs.subfile,
                             usecols=Iargs.subcol,unpack=True)
        pickle.dump(data,open('wham_pckl/window_%d.pckl'%window,'wb'))
        if Iargs.deriv == 1:
            for key in etypes:
                en[key]=np.genfromtxt(str(window)+"/"+key+"_init.out",usecols=0)
            pickle.dump(en,open("wham_pckl/ener_%d.pckl"%window,'wb'))


if __name__ == "__main__":

    parser=argparse.ArgumentParser()
    parser.add_argument('-subfile', default="lif.distance", type=str,
                         help="File name for colvar in subdirectory")
    parser.add_argument('-subcol',default=1, type=int, 
                         help="Column of subfile for colvar")
    parser.add_argument('-deriv',default=0, type=int,
                         help="[0] Turn off derivative calculation (default) [1] Turn on derivative calculation")
    parser.add_argument('-enerfile', default="flucts.inp", type=str,
                         help="File name with energy information, default flucts.inp")
    parser.add_argument('-metafile', default='wham_metadata.info', type=str,
                         help='File name with directory, loc, and k info.')
    Iargs = parser.parse_args()

    Main(Iargs)