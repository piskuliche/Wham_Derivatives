#!/usr/bin/env python
import os, shutil
import numpy as np

def Conv_Files(Iargs):
    """Converts distance units for subdirectory files.

    Args:
        Iargs (argparse): Input arguments.

    """
    from conv_metadata import Read_Meta

    col, x, k = Read_Meta()

    bohr_to_angstrom = 0.529177
    
    nwindows = len(x)

    if Iargs.xunits not in ['bohr','angstrom']:
        exit("Error: xunits must be bohr or angstrom")

    for i in range(nwindows):
        t, xvals = [], []
        try:
            t, xvals = np.genfromtxt("%d/%s" % (i,Iargs.subfile), usecols=(0,1), unpack=True)
        except:
            raise OSError("Error: Couldn't find subfile")
        new_x = xvals[Iargs.start:Iargs.stop:Iargs.skip]
        t = t[Iargs.start:Iargs.stop:Iargs.skip]
        if Iargs.xunits == 'bohr':
            new_x = new_x*bohr_to_angstrom
        print(np.shape(t),np.shape(new_x))
        np.savetxt("%d/final.colvars"%i, np.c_[t,new_x])
    


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-xunits', default='bohr',type=str,
                         help='Distance units type [bohr, angstrom]')
    parser.add_argument('-subfile', default='hartree',type=str,
                         help='Energy units type [hartree, kcal]')
    parser.add_argument('-subcol', default=1,type=int,
                         help='Column from which file to read.')                     
    parser.add_argument('-start', default=0, type=int,
                         help='Starting index for data')
    parser.add_argument('-stop', default=-1, type=int,
                         help='Stopping index for data')
    parser.add_argument('-skip', default=1, type=int,
                         help='Skipping index for data')                                          

    Iargs = parser.parse_args()

    Conv_Files(Iargs)