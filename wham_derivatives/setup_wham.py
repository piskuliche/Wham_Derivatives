#!/usr/bin/env python

import os, shutil, sys
import numpy as np


def Generate_Forceconsts():
    """Writes force constants out to a file

    This function takes user input, and turns it into a file which includes all the force constants.

    """
    f = open('force.consts','w')
    num_consts = int(input('How many force consts?\n'))
    prev_bins = 0
    for i in range(num_consts):
        num_bins = int(input('How many bins for const %s?\n' % str(i+1)))
        start = float(input('What is the starting distance (in A)?'))
        end = float(input('What is the ending distance (in A)?'))
        sep = end - start
        kval = float(input('What is the force constant (in hartree/bohr^2 aka cp2k internal)?\n'))
        for j in range(num_bins):
            ro = start+j*sep/float(num_bins)
            f.write('%s %s %s\n' % (prev_bins+1+j, kval, ro))
        prev_bins += num_bins

    f.close()

def Setup_CP2K():
    """Takes the force.consts file and setups cp2k simulation

    This writes the file wham_metadata.info, run_windows.sh, the collective.inc file in
    each windows directory.

    """
    K, ro = np.genfromtxt('force.consts', usecols=(1,2), unpack=True)
    bins = len(K)
    
    angpau=0.529177249
    kcalmolphartree=627.509
    conv = kcalmolphartree/angpau**2
    meta = open('wham_metadata.info', 'w')
    sub = open('run_windows.sh','w')

    sub.write('#!/bin/bash\n')
    sub.write('#SBATCH --job-name=wham_bin\n')
    sub.write('#SBATCH --partition=sixhour\n')
    sub.write('#SBATCH --output=wham_bin-%A_%a.out\n')
    sub.write('#SBATCH --nodes=1\n')
    sub.write('#SBATCH --ntasks-per-node=10\n')
    sub.write('#SBATCH --constraint=intel\n')
    sub.write('#SBATCH --mem=50G\n')
    sub.write('#SBATCH --time=06:00:00\n')
    sub.write('#SBATCH --array 0-%s\n\n\n\n' % (bins))

    sub.write('module load cp2k/6.1/popt\n\n')


    sub.write('cd $SLURM_ARRAY_TASK_ID\n')
    sub.write('mpirun -np 10 cp2k.popt inp_const.cp2k > run1.new\n')
    sub.write("sed 's/\([ \t]\+[^ \t]*\)\{3\}$//' out.colvar > lif.distance\n")
    sub.write('cd ../\n')

    for i in range(0,bins):
        if not os.path.isdir(str(i)):
            os.mkdir(str(i))
        f = open(str(i)+'/collective.inc','w')
        f.write('\t&COLLECTIVE\n')
        f.write('\t COLVAR 1\n')
        f.write('\t INTERMOLECULAR TRUE\n')
        f.write('\t TARGET [angstrom] %s\n' % ro[i])
        f.write('\t\t&RESTRAINT\n')
        f.write('\t\t K %s\n' % K[i])
        f.write('\t\t&END RESTRAINT\n')
        f.write('\t&END COLLECTIVE\n')
        f.close()
        shutil.copyfile('inp_const.cp2k',str(i)+'/inp_const.cp2k')
        #shutil.copyfile('conv.py',str(i)+'/conv.py')
        # Outputs in angstrom, kcal/(mol*angstrom)
        meta.write('%s/output.distance %s %s\n' % (str(i), ro[i] , K[i]*conv))
        
    meta.close()
    sub.close

def Setup_LAMMPS():
    """Takes the force.consts file and setups LAMMPS simulation

    This writes the file wham_metadata.info, run_windows.sh, the collective.inc file in
    each windows directory.

    """
    K, ro = np.genfromtxt('force.consts', usecols=(1,2), unpack=True)
    bins = len(K)
    angpau=0.529177249
    kcalmolphartree=627.509
    conv = kcalmolphartree/angpau**2
    meta = open('wham_metadata.info', 'w')
    sub = open('run_windows.sh','w')
    sub.write('#!/bin/bash\n')
    sub.write('#SBATCH --job-name=wham_bin\n')
    sub.write('#SBATCH --partition=sixhour,shontz,laird,thompson\n')
    sub.write('#SBATCH --constraint=intel\n')
    sub.write('#SBATCH --output=output.log\n')
    sub.write('#SBATCH --nodes=1\n')
    sub.write('#SBATCH --ntasks-per-node=20\n')
    sub.write('#SBATCH --time=6:00:00\n')
    sub.write('#SBATCH --array 0-%s\n\n\n' % (bins-1))
    sub.write('module purge\n')
    sub.write('module load lammps/lmps2021\n\n')


    sub.write('cd $SLURM_ARRAY_TASK_ID\n')
    sub.write('mpirun lmp_mpi <in.ip \n')
    sub.write('pull_energies.py log.production\n')
    sub.write('rm log.production log.lammps *.BAK P*out Temp*out Time*out\n')
    sub.write('cd ../\n')

    for i in range(0,bins):
        if not os.path.isdir(str(i)):
            os.mkdir(str(i))
        f = open(str(i)+'/lif.bias','w')
        """
        with open("lif.colvar",'r') as g:
            lines = g.readlines()
            for line in lines:
                line = line
                if "RRR" in line:
                    f.write("  centers %s %s\n" % (ro[i],ro[i]))
                elif "FFF" in line:
                    f.write("  forceConstant %s\n" % (K[i]*conv*2.0))
                else:
                    f.write("%s" % line)
        """
        f.write("fix cv all restrain bond 1186 1187 %s %s %s %s\n" % (K[i]*conv, K[i]*conv, ro[i], ro[i]))
        f.write("fix_modify cv energy yes\n")
        f.close()
        shutil.copyfile('in.ip',str(i)+'/in.ip')
        shutil.copyfile('lmps.include.coeffs',str(i)+'/lmps.include.coeffs')
        
        # Outputs in angstrom, kcal/(mol*angstrom)
        meta.write('%s/lif.distance %s %s\n' % (str(i), ro[i], K[i]*conv/2.0))

    meta.close()
    sub.close

if __name__ == "__main__":
    # runs and chooses the program.
    if len(sys.argv) == 1:
        print("Usage: python setup_wham.py -type [1,2 or 3] ")
        sys.exit("For more information run program with '-h' option")
    if "-h" in sys.argv:
        print("Usage: python setup_wham.py -type [name] ")
        print("This code sets up the wham calculation and sets the force constants for each window.")
        print("type: 1) generate force constants, 2) cp2k setup, 3) lammps setup")
        print("Note run 1 then 2 or 1 then 3 not 1 2 3.")
        sys.exit("Exiting")
    if "-type" in sys.argv:
        index = sys.argv.index("-type")+1
        type = int(sys.argv[index])
    else:
        sys.exit("type must be specified")

    if type == 1:
        Generate_Forceconsts()
    elif type == 2:
        Setup_CP2K()
    elif type == 3:
        Setup_LAMMPS()
    else:
        print("Invalid options chosen - try again")

