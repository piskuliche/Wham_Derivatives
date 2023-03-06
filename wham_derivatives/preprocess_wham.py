#!/usr/bin/env python
import argparse, os, pickle
import numpy as np

def Write_WHAM_Files(metafile="wham_metadata.info", subfile="final.colvars",subcol=1,deriv=0, enerfile="flucts.inp" ):
    """Writes files in a very simple, WHAM-readable format. 

    This takes files, and writes them in a very WHAM-readable format, 
    mainly by converting them to pckl files in a centralized directory.

    Iargs.deriv sets whether this also reads energies (and writes those to pickle files).

    Args:
        metafile (str): Name of metadata file, [default=wham_metadata.info]
        subfile (str): Name of subdirectory file, [default=final.colvars]
        subcol (int): Column of subfile to read, [default=1]
        deriv (int): Whether to read energies, [default=0]
        enerfile (str): Name of energy file, [default=fluct.inp]
    
    Raises:
        OSError: Input files provided were incorrect, or in the wrong format.

    """
    if not os.path.exists("wham_pckl"):
        os.makedirs('wham_pckl')
    
    nwindows=0
    try:
        k=np.genfromtxt(metafile,usecols=1,unpack=True)
        nwindows = len(k)
    except:
        raise OSError("Error: Trouble grabbing windows from metafile")
    
    etypes = []
    if deriv == 1:
        try:
            etypes=np.genfromtxt(enerfile, usecols=0, dtype=str)
            print("Grabbing energies")
            print("Selected energy types:",*etypes)
        except:
            raise OSError("Error: Trouble grabbing energies")

    svdata = []
    # Loops over windows here to write pickle files
    for window in range(nwindows):
        if window%10 == 0: print(window)
        data, en = [], {}
        
        data = np.genfromtxt(str(window)+"/"+subfile,
                             usecols=subcol,unpack=True)
                             
        pickle.dump(data,open('wham_pckl/window_%d.pckl'%window,'wb'))

        if deriv == 1:
            for key in etypes:
                en[key]=np.genfromtxt(str(window)+"/"+key+"_init.out",usecols=0)
            pickle.dump(en,open("wham_pckl/ener_%d.pckl"%window,'wb'))


if __name__ == "__main__":

    parser=argparse.ArgumentParser()
    parser.add_argument('-subfile', default="final.colvars", type=str,
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

    Write_WHAM_Files(metafile=Iargs.metafile, subfile=Iargs.subfile, subcol=Iargs.subcol, deriv=Iargs.deriv, enerfile=Iargs.enerfile)