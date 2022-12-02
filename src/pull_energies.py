#!/usr/bin/env python
import numpy as np
import sys

def Pull_LAMMPS(Iargs, n):
    """Pulls energy data from LAMMPS

    Pulls energy for all columns in the LAMMPS log file and writes them
    as a file for each one.

    Args:
        n (int): Number of windows

    """
    data={}
    with open("%d/%s"%(n, Iargs.logfile), 'r') as f:
        lines=f.readlines()
        flag=0
        keys=[]
        for line in lines:
            if "Loop" in line:
                flag=0
                print("stop")
            if flag == 1 and "colvars" not in line:
                for key in keys:
                    data[key].append(float(line.strip().split()[loc[key]]))
            if "Step Time" in line:
                flag=1
                data={}
                loc={}
                keys=line.strip().split()
                count = 0
                for key in keys:
                    data[key]=[]
                    loc[key]=count
                    count+=1
                print("start")

    for key in data:
        data[key].pop()
        np.savetxt("%d/%s_init.out"%(n, key), np.c_[data[key]])
        if key == "Volume":
            L = np.array(data[key])**(1./3.)
            np.savetxt("%d/L.dat"%(n), np.c_[L])
    
def Pull_CP2K(Iargs, n):
    """File for pulling energies from CP2K

    CP2K Energies are pulled using numpy. The data is pulled using numpy,
    and is only pulled for columns from enerfile. Also has the option
    to do unit conversions as needed when saving the energies to their
    files. To do so, use third column of flucts.inp file

    Args:
        n (int): Number of windows

    """
    hartree_to_kcal = 627.509
    try:
        cols, conv = np.genfromtxt(Iargs.enerfile, usecols=(1, 2), unpack=True, dtype=None)
        keys = np.genfromtxt(Iargs.enerfile, usecols=0, unpack=True, dtype=str)
    except:
        raise OSError("Error: Trouble reading %s"%Iargs.enerfile)
    print(keys)
    data=np.genfromtxt("%d/%s"%(n, Iargs.logfile), unpack=True)
    for i, key in enumerate(keys):
        energy = data[cols[i]]
        if conv[i] == 1:
            energy = energy * hartree_to_kcal
        np.savetxt("%d/%s_init.out"%(n, key), np.c_[energy])

def Pull_Energies(Iargs):
    try:
        k=np.genfromtxt(Iargs.metafile,usecols=1,unpack=True)
        nwindows = len(k)
    except:
        exit("Error: Trouble grabbing windows from metafile")
    
    nwindows = len(k)

    if Iargs.program == "LAMMPS":
        for n in range(nwindows):
            Pull_LAMMPS(Iargs, n)
    elif Iargs.program == "CP2K":
        for n in range(nwindows):
            Pull_CP2K(Iargs, n)
    else:
        raise ValueError("Error: Improper program selection")




if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-program', default="LAMMPS", type=str,
                         help="Which program [CP2K or LAMMPS]")
    parser.add_argument('-logfile', default="log.production", type=str,
                         help="File name for log file")
    parser.add_argument('-enerfile', default="flucts.inp", type=str,
                         help="File name with energy information, default flucts.inp")
    parser.add_argument('-metafile', default='wham_metadata.info', type=str,
                         help='File name with directory, loc, and k info.')
    Iargs = parser.parse_args()

    Pull_Energies(Iargs)